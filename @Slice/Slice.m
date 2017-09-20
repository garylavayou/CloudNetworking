%%  Network Slice
% Support network resource allocation scheme.
%%
classdef Slice < VirtualNetwork
    % Specify the properties that can only be modified by Physcial Network directly
    properties (SetAccess = {?Slice,?CloudNetwork,?SliceFlowEventDispatcher})
        Variables;      % [x,z]
        x_path;             % allocated bandwidth of paths
        z_npf;              % allocated resource on nodes
        alpha_f;
        weight;             % weight for the slice user's utility.
        I_path_function;        % used by <singleSliceOptimization>.
    end
    
    properties (Access = protected)
        As_res;             % coefficient matrix of the processing constraints
        x0;                 % start point
        link_load;          % temporary results of link load.
        node_load;          % temporary results of node load.
        flow_rate;          % temporary results of flow rate.
    end
    
    properties(Dependent, SetAccess = protected, GetAccess={?Slice, ?CloudNetwork})
        num_vars;           % number of optimization variables in the problem
        num_varz;           % number of variables in vector z_npf
        num_lcon_res;       % number of linear constraints
        % local_path_id;
    end
        
    methods
        function this = Slice(slice_data)
            if nargin == 0
                slice_data = [];
            end
            this@VirtualNetwork(slice_data);
            if nargin == 0
                return;
            end
            %%%
            % Link price
            if isfield(slice_data, 'link_price')
                this.VirtualLinks.Price = slice_data.link_price;
            else
                this.VirtualLinks.Price  = zeros(this.NumberVirtualLinks,1);
            end
            %%%
            % Data center node price
            if isfield(slice_data, 'node_price')
                this.VirtualDataCenters.Price = slice_data.node_price;
            else
                this.VirtualDataCenters.Price = zeros(this.NumberDataCenters,1);
            end
            
            if isfield(slice_data,'Weight')
                this.weight = slice_data.Weight;
            end
            
            if isfield(slice_data, 'alpha_f')
                this.alpha_f = slice_data.alpha_f;
            end
            this.getAs_res;
        end
    end
    
    methods
        function n = get.num_vars(this)
            n = (this.NumberVNFs*this.NumberDataCenters+1)*this.NumberPaths;
        end
        
        function n = get.num_varz(this)
            n = this.NumberVNFs*this.NumberDataCenters*this.NumberPaths;
        end
        
        function n = get.num_lcon_res(this)
            n = size(this.As_res,1);
        end
        
        function setPathBandwidth(this, x_path)
            if nargin == 1
                x_path = this.Variables.x;
            end
            p = 1;
            for i = 1:this.NumberFlows
                pathlist = this.FlowTable{i,'Paths'}.paths;
                for l = 1:length(pathlist)
                    pathlist{l}.bandwidth = x_path(p);
                    p = p + 1;
                end
            end
        end
                
        function r = getFlowRate(this, path_vars)
            if nargin == 1
                r = this.I_flow_path * this.Variables.x;
            else
                r = this.I_flow_path * path_vars;
            end
        end
        
        function ye = getLinkLoad(this, path_vars)
            % retrive the link load of the slice, given the path variables.
            if nargin == 1
                ye = this.I_edge_path * this.Variables.x;
            else
                ye = this.I_edge_path * path_vars;
            end
        end
        
        %         function vn = getNodeLoad(this, node_vars)
        %             if nargin == 1
        %                 node_vars = this.Variables.z;
        %             end
        %             %%
        %             % |node_vars| is index by |(node,path,function)|.
        %             % node_load = sum(f, node_vars(:,:,f).*I_node_path).
        %             NN = this.NumberVirtualNodes;
        %             NP = this.NumberPaths;
        %             NV = this.NumberVNFs;
        %             vn = zeros(NN,1);
        %             np = NN * NP;
        %             z_index = 1:np;
        %             for i = 1:NV
        %                 node_vars_fi = reshape(node_vars(z_index), NN, NP);
        %                 vn = vn + sum(this.I_node_path.*node_vars_fi,2);
        %                 z_index = z_index + np;
        %             end
        %%%
        % retrive the node load of the slice, given the node variables. The node variables
        % represent the resource allocation of data center nodes.
        %
        % |v_n|: data center's resource consumption.
        function v_n = getNodeLoad(this, node_vars)
            if nargin == 1
                node_vars = this.Variables.z;
            end
            %%
            % |node_vars| is index by |(node,path,function)|.
            % node_load = sum(f, node_vars(:,:,f).*I_node_path).
            NC = this.NumberDataCenters;
            NP = this.NumberPaths;
            v_n = zeros(NC,1);
            np = NC * NP;
            z_index = 1:np;
            for i = 1:this.NumberVNFs
                node_vars_fi = reshape(node_vars(z_index), NC, NP);
                v_n = v_n + sum(this.I_node_path.*node_vars_fi,2);
                z_index = z_index + np;
            end
            
            %% Alternative way to compute the node load.
            %             col_index = (1:NC:((NP-1)*NC+1))';
            %             col_index = repmat(col_index, 1, NV);
            %             for c = 2:NV
            %                 col_index(:,c) = col_index(:,c-1) + NC*NP;
            %             end
            %             col_index = col_index(:);
            %             As = zeros(NC, NP*NV);
            %             for row_index = 1:NC
            %                 As(row_index, col_index) = repmat(this.I_node_path(row_index,:),1, NV);
            %                 col_index = col_index + 1;
            %             end
            %             vn = As*node_vars;
        end
                
        function c = link_unit_cost(this)
            % the virtual links's unit cost
            c = this.Parent.getLinkField('UnitCost', this.VirtualLinks.PhysicalLink); 
        end
        
        function c = node_unit_cost(this)
            % the virtual data center nodes's unit cost
            c = this.Parent.getNodeField('UnitCost', this.getDCNI);   
        end
              
        function sc = getSliceCost(this, node_load, link_load, model)
            sc = this.getResourceCost(this, node_load, link_load, model);
        end
        
        function r = getRevenue(this)
            if isempty(this.flow_rate)
                r = this.weight*sum(fcnUtility(this.FlowTable.Rate));
            else
                r = this.weight*sum(fcnUtility(this.flow_rate));
            end
        end
        
        function b = checkFeasible(this, vars, options)
            if nargin <= 1 || isempty(vars)
                vars = [this.Variables.x; this.Variables.z];
            end
            if nargin <= 2 || ~isfield(options, 'Display')
                options.Display = 'off';
            end
            if nargin <= 2 || ~isfield(options, 'Tolerance')
                options.Tolerance = 10^-3;
            end
            
            ub = this.As_res*vars;
            if isempty(find(ub > options.Tolerance, 1))
                b = true;
            else
                b = false;
                if ~strcmp(options.Display, 'off')
                    if issparse(ub)
                        cprintf('Comments', 'sparse vector\n');
                    end
                    warning('checkFeasible: Maximal violation is %.4f.', full(max(ub)));
                end
            end
        end
        
        profit = optimalFlowRate( this, options );
        [utility, node_load, link_load] = priceOptimalFlowRate(this, x0, options);
        [payment, grad, pseudo_hess] = fcnLinkPricing(this, link_price, link_load);
        [payment, grad, pseudo_hess] = fcnNodePricing(this, node_price, node_load);        
    end
    
    methods (Static)
        %% TODO: move to Network and split it according to type
        slice_template = loadSliceTemplate(index);

        % Objective function and gradient
        [profit, grad]= fcnProfit(var_x, S, options);
        [profit, grad] = fcnSocialWelfare(x_vars, S, options);

        % Hessian matrix
        hess = fcnHessian(var_x, ~, S, options);

    end
    
    methods (Access = private)
        %         function sc = getResourceCost(this, node_load, link_load, model)
        %             if nargin <= 1 || isempty(node_load)
        %                 node_load = this.VirtualNodes.Load;
        %             end
        %             if nargin <= 2 || isempty(link_load)
        %                 link_load = this.VirtualLinks.Load;
        %             end
        %             if nargin <=3
        %                 warning('model is set as Approximate.');
        %                 model = 'Approximate';
        %             end
        %
        %             pn = this.Parent;
        %             link_uc = pn.getLinkField('UnitCost', this.VirtualLinks.PhysicalLink); % the virtual links's unit cost
        %             node_uc = pn.getNodeField('UnitCost', this.VirtualNodes.PhysicalNode); % the virtual nodes's unit cost
        %             epsilon = pn.unitStaticNodeCost;
        %
        %             if strcmp(model, 'Approximate')
        %                 sc = dot(link_uc, link_load) + dot(node_uc, node_load) ...
        %                     + pn.phis_n*sum(node_load)+pn.phis_l*sum(link_load)+...
        %                     pn.static_factor*epsilon/pn.NumberSlices;
        %             elseif strcmp(model, 'Accurate')
        %                 sc = dot(link_uc, link_load) + dot(node_uc, node_load);
        %             else
        %                 error('error: invalid model %s', model);
        %             end
        %         end
                %%%
        % getSliceCost  When compute the static cost, the Capacity of all physical nodes
        % and links is included.
        % |epsilon/pn.NumberSlices| is a constant. To keep consistence with other methods,
        % this part should not be ignored. See also getNetworkCost and getStaticCost.
        % The calculation is not absolutely precise, since it cannot be decide that the
        % static cost should be arributed to which slices.
        %
        % When calculate network cost as a single slice, this method equals to
        % _getNetworkCost_ .
        function rc = getResourceCost(this, node_load, link_load)
            if nargin <= 1 || isempty(node_load)
                node_load = this.VirtualDataCenters.Load;
            end
            if nargin <= 2 || isempty(link_load)
                link_load = this.VirtualLinks.Load;
            end
            
            %% A temporary slice should be assigned the parent network
            % so that |this.Parent| is valid.
            pn = this.Parent;
            link_uc = pn.getLinkCost(this.VirtualLinks.PhysicalLink);
            node_uc = pn.getNodeCost(this.getDCPI);
            %% Accurate Model
            % A slice cannot decide how to devide the static cost between slices.
            % one method is devision by usage, using the physic network's load data
            % of each slice. So this function only calculate the dynamic part of the
            % cost, and the slices should further calculate the static cost outside
            % this method.
            rc = dot(link_uc, link_load) + dot(node_uc, node_load);
        end
    end
    
    methods (Access = {?Slice, ?CloudNetwork})
        %% Linear constraint without consdiering the bound constraint.
        % |As| takes the following form
        %
        % $$\left[ \begin{array}{ccccc}
        %   I_1     & H_s &      &          &       \\
        %   I_2     &     & H_s  &          &       \\
        %   \vdots  &     &      &  \ddots  &       \\
        %   I_F     &     &      &          & H_s
        % \end{array} \right]
        % \left[ \begin{array}{c}x\\z_1\\z_2\\\vdots\\z_f\\\vdots\\z_F\end{array}\right] $$
        %
        % where
        %
        % $$I_f = {\left[ \begin{array}{cccc}
        % \alpha_f &          &        &         \\
        %          & \alpha_f &        &         \\
        %          &          & \ddots &         \\
        %          &          &        & \alpha_f
        % \end{array} \right]}_{P\times P},$$
        % $$H_s = {\left[ \begin{array}{cccccccccc}
        % -h_{1,1} & \cdots & -h_{N,1} &          &        &          &        &         &        & \\
        %          &        &          & -h_{1,2} & \cdots & -h_{N,2} &        &         &        & \\
        %          &        &          &          &        &          & \ddots &         &        & \\
        %          &        &          &          &        &          &        &-h_{1,P} & \cdots & -h_{N,P}
        % \end{array} \right]}_{P\times NP},$$
        % $$ x = \left[\begin{array}{c}x_1\\x_2\\ \vdots\\x_P\end{array}\right]$$
        % $$ z_f = \left[\begin{array}{c}z_{1,1,f}\\ \vdots\\ z_{N,1,f}\\z_{1,2,f}\\
        %   \vdots\\ z_{N,2,f}\\ \vdots\\z_{N,P,F}\end{array}\right]$$
        %
        % According to the martix formulation, the number of non-zero elements in |As|
        % is equal to |F*(P+nnz(Hs))|, where |nnz(Hs)| is equal to the number of
        % nonzero elements in |I_node_path|.
        %
        % NOTE: it is not necessary to evaluate As_res each time when visiting it. so, we
        % define a normal function to update the property |As_res|.
        function As = getAs_res(this)
            NC = this.NumberDataCenters;
            NP = this.NumberPaths;
            NV = this.NumberVNFs;
            nnz_As = NV*(NP+nnz(this.I_node_path));
            num_lcon = NP*NV;
            % old_As = this.As_res
            As = spalloc(num_lcon, this.num_vars, nnz_As);
            row_index = 1:NP;
            col_index = NP+(1:NC);
            for f = 1:NV
                if nargin >= 2 && ~isempty(this.alpha_f) 
                    % used when treat all VNFs as one function.
                    for p = 1:NP
                        As(p,p) = this.alpha_f(this.path_owner(p)); %#ok<SPRIX>
                    end
                else
                    af = this.Parent.VNFTable{this.VNFList(f),{'ProcessEfficiency'}};
                    As(row_index,1:NP) = af * eyesparse(NP); %#ok<SPRIX>
                end
                for p = 1:NP
                    As(row_index(p), col_index) = -this.I_node_path(:,p); %#ok<SPRIX>
                    col_index = col_index + NC;
                end
                row_index = row_index + NP;
            end
            this.As_res = As;
        end

    end
  
    methods (Access = protected)
        function [x, fval, exitflag] = ...
                obj_fun(this, x0, As, bs, Aeq, beq, lbs, ub, lb, fmincon_opt, method)
            if strfind(method, 'price') % 'price', 'slice-price'
                [x, fval, exitflag] = fmincon(@(x)Slice.fcnProfit(x,this), ...
                    x0, As, bs, Aeq, beq, lbs, ub, lb, fmincon_opt);
            else  % 'normal', 'single-function', ...
                [x, fval, exitflag] = ...
                    fmincon(@(x)Slice.fcnSocialWelfare(x,this), ...
                    x0, As, bs, Aeq, beq, lbs, ub, lb, fmincon_opt);
            end
        end
    end
end

%% Methods     
% * *priceOptimalFlowRate* : find the optimal flow rate that maximizing the net profit of
% the network slice.
%
%      [rate, net_profit, node_load, link_load] = priceOptimalFlowRate(this, x0, options)
%
% * *priceOptimalFlowRateCompact* : find the optimal flow rate that maximizing the net profit
% of the network slice.
%
%      x = priceOptimalFlowRateCompact(this, x0)
%
% The optimizition procedure in this method remove the unnecessary components from the
% independent variable |x|, so that the problem scale is cut down.
%
% * *fcnProfit* |static| : Evalute the objective function and gradient.
%
%      [profit, grad]= fcnProfit(var_x, S)
%
% |grad|: the gradient value of the objective function.
% The upper bound number of non-zero elements in the gradient vector: the gradient on path
% variable is nonzeros, so there is |P| components; whether the gradient on node variable
% is zeros is decided by the node-path incidence matrix, i.e. |nnz(I_node_path)*F|.
%
% * *fcnHessian* |static| : Hessian matrix of the Largrangian.
%
%      hess = fcnHessian(var_x, ~, S)
%
% Since the problem only contains linear constraint, the hessian matrix of the
% Largrangian is equal to the second derivatives of the objective function, and the
% Largrangian multipliers $\lambda$ takes no effect.
% The Hessian matrix contains only $P^2$ nonzeros elements on the diagonal,
% which is the second derviatives on path variables.