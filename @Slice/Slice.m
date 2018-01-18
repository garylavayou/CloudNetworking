%%  Network Slice
% Support network resource allocation scheme.
%%
classdef Slice < VirtualNetwork & EventReceiver
    % Specify the properties that can only be modified by Physcial Network directly
    properties (SetAccess = {?Slice,?CloudNetwork,?SliceFlowEventDispatcher})
        temp_vars = struct; % [x,z]: temporary variables that will not be rounded.
        Variables;          % [x,z]: final results from |temp_vars|, with rounding.
        weight;             % weight for the slice user's utility.
        As_res;             % coefficient matrix of the processing constraints
        prices;
    end
    
    properties (Access = protected)
        % coefficient used to compute node resource consumption.
        % is a function of node-path incident matrix (I_node_path).
        Hrep;
        % coefficient used to compute VNF intance capacity
        % is a function of node-path incident matrix (I_node_path).
        Hdiag;              
        x0;                 % start point
        link_load;          % temporary results of link load.
        node_load;          % temporary results of node load.
        flow_rate;          % temporary results of flow rate.
        I_active_variables;  % indicator of active variables
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
            if isfield(slice_data, 'LinkPrice')
                this.VirtualLinks.Price = slice_data.LinkPrice;
            else
                this.VirtualLinks.Price  = zeros(this.NumberVirtualLinks,1);
            end
            %%%
            % Data center node price
            if isfield(slice_data, 'NodePrice')
                this.VirtualDataCenters.Price = slice_data.NodePrice;
            else
                this.VirtualDataCenters.Price = zeros(this.NumberDataCenters,1);
            end
            
            if isfield(slice_data,'Weight')
                this.weight = slice_data.Weight;
            end
            
            this.options = structmerge(this.options, ...
                getstructfields(slice_data, 'SlicingMethod', 'default', 'dynamic-slicing'),...
                getstructfields(slice_data, 'PricingPolicy', 'ignore'));

            this.getAs_res;
            this.getHrep;
            this.getHdiag;
        end
        
        function finalize(this, node_price, link_price)
            if this.NumberFlows == 0
                this.Variables.x = [];
                this.Variables.z = [];
                this.VirtualDataCenters{:,'Load'} = 0;
                this.VirtualLinks{:,'Load'} = 0;
                this.VirtualDataCenters{:,'Capacity'} = 0;
                this.VirtualLinks{:,'Capacity'} = 0;
            else
                this.Variables.x = this.temp_vars.x;
                this.Variables.z = this.temp_vars.z;
                % subclass inctance may have different implementation of _postProcessing_, it
                % is recommended, that the subclass maintains the default behavior of this
                % function.
                this.postProcessing();
                this.VirtualDataCenters.Load = this.getNodeLoad;
                this.VirtualLinks.Load = this.getLinkLoad;
                this.VirtualDataCenters.Capacity = this.VirtualDataCenters.Load;
                this.VirtualLinks.Capacity = this.VirtualLinks.Load;
                this.FlowTable.Rate = this.getFlowRate;
                this.setPathBandwidth;
            end
            if nargin >= 2 && ~isempty(node_price)
                this.VirtualDataCenters.Price = node_price(this.getDCPI);
            end
            if nargin >= 3 && ~isempty(link_price)
                this.VirtualLinks.Price = link_price(this.VirtualLinks.PhysicalLink);
            end
        end
    end
    
    methods
        % equation:
        %    $\sum_{p,f}{b_{np}z_{npf}} \le V_{n}$
        %
        % For the _node capacity constraint_, the coefficient matrix is filled row-by-row.
        % On row i, the non-zero elements located at (i-1)+(1:NC:((NP-1)*NC+1)) for the
        % first |NC*NP| columns, and then the first |NC*NP| columns are duplicated for
        % |NV| times, resulting in |NC*NP*NV| columns.
        %
        % The base matrix is the same as in <Hdiag>.
        function H = getHrep(this)
            NP = this.NumberPaths;
            NC = this.NumberDataCenters;
            NV = this.NumberVNFs;
            
            H = spalloc(NC, this.num_varz, nnz(this.I_node_path)*NV);
            col_index = (1:NC:((NP-1)*NC+1))';      % the vnf variable index related to the first node
            col_index = repmat(col_index, 1, NV);   % derive other node's vnf variable index
            for v = 2:NV
                col_index(:,v) = col_index(:,v-1) + NC*NP;
            end
            col_index = col_index(:);
            for n = 1:NC    % all path that use node n is counted
                H(n, col_index) = repmat(this.I_node_path(n,:), 1, NV);  %#ok<SPRIX>
                col_index = col_index + 1;
            end
            this.Hrep = H;
        end

        % equation:
        %    $\sum_{p}{b_{np}z_{npf}} \le v_{nf}$
        %
        % For the node capacity constraint for VNF |f|, the coeffiect matrix Hs is
        % filled block-by-block. In the |k|th block (NC*NC), we put the |k|th column
        % of H_np to the diagnal of the block. Then |H_np| is duplicated into a larger
        % diagonal corresponding to all VNFs.
        %
        % The capacity of VNF instance |VNFCapacity|, has been calculate after
        % allocating resource to the slice. Some component of |VNFCapacity| might be
        % zero, since VNF f may not be located at node n.
        %
        % See also <Hrep>.
        function Hs = getHdiag(this)
            global DEBUG; %#ok<NUSED>
            NC = this.NumberDataCenters;
            NP = this.NumberPaths;
            Hs_np = spalloc(NC, NC*NP, nnz(this.I_node_path));
            col_index = 1:NC;
            for p = 1:NP
                Hs_np(1:NC,col_index) = diag(this.I_node_path(:,p));  %#ok<SPRIX>
                col_index = col_index + NC;
            end
            Hs = block_diag(Hs_np, this.NumberVNFs);
            this.Hdiag = Hs;
        end

        function n = get.num_vars(this)
            n = (this.NumberVNFs*this.NumberDataCenters+1)*this.NumberPaths;
        end
        
        function n = get.num_varz(this)
            n = this.NumberVNFs*this.NumberDataCenters*this.NumberPaths;
        end
        
        function n = get.num_lcon_res(this)
            n = size(this.As_res,1);
        end
        
        function setPathBandwidth(this, x)
            if nargin == 1
                x = this.Variables.x;
            end
            p = 1;
            for i = 1:this.NumberFlows
                pathlist = this.FlowTable{i,'Paths'}.paths;
                for l = 1:length(pathlist)
                    pathlist{l}.bandwidth = x(p);
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
           
        %% Capacity and Load
        % Public interface for network to inquire the resource occupation of the slices.
        % In class <Slice>, <getLinkCapacity> and <getNodeCapacity> are equal to the
        % protected methods <getLinkLoad> and <getNodeLoad> respectively. But in
        % subclasses of <Slice>, the slice load might be less than its capacity, so that
        % the two group of methods return different results.
        function c = getLinkCapacity(this, path_vars)
            if nargin == 1
                c = this.getLinkLoad;
            else
                c = this.getLinkLoad(path_vars);
            end
        end
        
        function c = getNodeCapacity(this, node_vars)
            if nargin == 1
                c = this.getNodeLoad;
            else
                c = this.getNodeLoad(node_vars);
            end
        end
        
        function vc = getVNFCapacity(this, z)
            %       znpf = reshape(full(this.Variables.z), this.NumberDataCenters, ...
            %       this.NumberPaths, this.NumberVNFs);
            %       znpf = znpf.* full(this.I_node_path);  % compatible arithmetic operation
            %       this.VNFCapacity = reshape(sum(znpf,2), this.NumberDataCenters*this.NumberVNFs,1);
            if nargin <= 1
                z = this.Variables.z;
            end
            vc = this.Hdiag * z;
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
        
        function b = checkFeasible(this, vars, opt_opts)
            if nargin <= 1 || isempty(vars)
                vars = [this.Variables.x; this.Variables.z];
            else
                vars = vars(1:this.num_vars);
            end
            if nargin >=3 && isfield(opt_opts, 'ConstraintTolerance')
                b = isempty(find(this.As_res*vars>opt_opts.ConstraintTolerance,1));
            else
                b = isempty(find(this.As_res*vars>1e-10,1));
            end
        end
                
        [payment, grad, pseudo_hess] = fcnLinkPricing(this, link_price, link_load);
        [payment, grad, pseudo_hess] = fcnNodePricing(this, node_price, node_load);    
        
        function eventhandler(~, ~, ~)
        end
    end

    methods(Static)
        %% TODO: move to Network and split it according to type
        slice_template = loadSliceTemplate(index);
    end
    
    methods (Static,Access=protected)
        % Objective function and gradient
        [profit, grad]= fcnProfit(vars, slice, options);
        %% Evaluate the objective function and gradient
        % only active independent variables are passed into the objective function.
        % Considering the constraint's coefficient matrix, if the corresponding column of the
        % coefficient matrix for a variable is all zero, then this variable is inactive and can be
        % directly set as 0. So we can remove it from the optimization problem.
        %
        % NOTE: we can also remove the all-zero rows of the coefficient matrix, which do not
        % influence the number of variables. See also <optimalFlowRate>.
        function [profit, grad]= fcnProfitCompact(act_vars, slice, options)
            vars = zeros(options.num_orig_vars,1);
            vars(slice.I_active_variables) = act_vars;
            
            % we extend the active variables by adding zeros to the inactive ones.
            if nargout <= 1
                profit = Slice.fcnProfit(vars, slice, options);
            else
                [profit, grad] = Slice.fcnProfit(vars, slice, options);
                % eliminate the inactive variable's derivatives.
                grad = grad(slice.I_active_variables);
            end
        end
        
        [profit, grad] = fcnSocialWelfare(x_vars, S);
        function [profit, grad] = fcnSocialWelfareCompact(act_vars, slice)
            vars = zeros(options.num_orig_vars,1);
            vars(slice.I_active_variables) = act_vars;
            
            if nargout <= 1
                profit = Slice.fcnSocialWelfare(vars, slice);
            else
                [profit, grad] = Slice.fcnSocialWelfare(vars, slice);
                grad = grad(slice.I_active_variables);
            end
        end

        % Hessian matrix
        hess = fcnHessian(var_x, ~, slice, options);
        %% Compact form of Hessian matrix of the Lagrangian
        function hess = fcnHessianCompact(act_vars, lambda, slice, options) %#ok<INUSL>
            vars = zeros(options.num_orig_vars,1);
            vars(slice.I_active_variables) = act_vars;
            if nargin == 3
                hess = Slice.fcnHessian(vars, [], slice);
            else
                hess = Slice.fcnHessian(vars, [], slice, options);
            end
            hess = hess(slice.I_active_variables, slice.I_active_variables);
        end
        
        function interpretExitflag(exitflag, message)
            global DEBUG INFO;
            if nargin <= 1
                message = '';
            else
                message = strtok(message, newline);
            end
            switch exitflag
                case 0
                    if ~isempty(DEBUG) && DEBUG
                        warning(message);    % max-iteration number exceed.
                    elseif ~isempty(INFO) && INFO
                        cprintf('SystemCommands', '%s\n', message);
                    end
                case 1
                    if ~isempty(INFO) && INFO
                        fprintf('(%d) %s\n', exitflag, message);
                    end
                case 2
                    if ~isempty(DEBUG) && DEBUG
                        warning('%s(%d)', message, exitflag);
                    elseif ~isempty(INFO) && INFO
                        cprintf('Comment', '%s(%d)\n', message, exitflag);
                    end
                case -3
                    error('error: Objective function unbounded below (%d). %s', exitflag, message);
                otherwise
                    error('error: Abnormal exit (%d). %s', exitflag, message);
            end
        end
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
                node_load = this.VirtualDataCenters.Capacity;
            end
            if nargin <= 2 || isempty(link_load)
                link_load = this.VirtualLinks.Capacity;
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
            if rc == 0
                warning('zero slice cost.');
            end
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
        function As = getAs_res(this, flow_owner, alpha_f)
            NC = this.NumberDataCenters;
            NP = this.NumberPaths;
            NV = this.NumberVNFs;       % For 'single-function', NV=1.
            nnz_As = NV*(NP+nnz(this.I_node_path));
            num_lcon = NP*NV;
            As = spalloc(num_lcon, this.num_vars, nnz_As);
            row_index = 1:NP;
            for f = 1:NV
                if nargin >= 3   % used when treat all VNFs as one function.
                    for p = 1:NP
                        As(p,p) = alpha_f(flow_owner(this.path_owner(p))); %#ok<SPRIX>
                    end
                else
                    af = this.Parent.VNFTable{this.VNFList(f),{'ProcessEfficiency'}};
                    As(row_index,1:NP) = af * eyesparse(NP); %#ok<SPRIX>
                end
                row_index = row_index + NP;
            end
            col_index = 1:NC;
            Hst = zeros(NP, NC*NP);     % |Hst| is a staircase-like diagnoal.
            for p = 1:NP
                Hst(p, col_index) = -this.I_node_path(:,p);
                col_index = col_index + NC;
            end
            As(:, (NP+1):end) = block_diag(Hst, NV);
            this.As_res = As;
        end

    end
  
    methods (Access = protected)
        % Called by _optimalFlowRate_, if 'Form=compact', options should be speicifed with
        % 'num_orig_vars'.
        function [x, fval] = optimize(this, params, options)
            global DEBUG;
            if contains(options.SlicingMethod, 'price') % 'price', 'slice-price'
                if strcmpi(options.Form, 'compact')
                    [x_compact, fval, exitflag, output] = ...
                        fmincon(@(x)Slice.fcnProfitCompact(x, this, options), ...
                        params.x0, params.As, params.bs, params.Aeq, params.beq, ...
                        params.lb, params.ub, [], options.fmincon_opt);
                else
                    [x, fval, exitflag, output] = ...
                        fmincon(@(x)Slice.fcnProfit(x, this, options), ...
                        params.x0, params.As, params.bs, params.Aeq, params.beq, ...
                        params.lb, params.ub, [], options.fmincon_opt);
                end
            else  % 'normal', 'single-function', ...
                if strcmpi(options.Form, 'compact')
                    [x_compact, fval, exitflag, output] = ...
                        fmincon(@(x)Slice.fcnSocialWelfareCompact(x,this), ...
                        params.x0, params.As, params.bs, params.Aeq, params.beq, ...
                        params.lb, params.ub, [], options.fmincon_opt);
                else
                    [x, fval, exitflag, output] = ...
                        fmincon(@(x)Slice.fcnSocialWelfare(x,this), ...
                        params.x0, params.As, params.bs, params.Aeq, params.beq, ...
                        params.lb, params.ub, [], options.fmincon_opt);
                end
            end
            this.interpretExitflag(exitflag, output.message);
            if ~isempty(DEBUG) && DEBUG
                assert(this.checkFeasible(x, ...
                    struct('ConstraintTolerance', options.fmincon_opt.ConstraintTolerance)), ...
                    'error: infeasible solution.');
            end
            
            if strcmpi(options.Form, 'compact')
                x = spalloc(this.num_vars, 1, nnz(x_compact));
                x(this.I_active_variables) = x_compact;
            end
        end
        
        %%%
        % Flow processing constraint might be violated, due to rounding of x,z. Three
        % candidate method can be applied for post processing:
        %       (1) round the x components, so that flow processing demand is meet;
        %       (2) drop the items in x that violates the constraint;
        %       (3) recover those z components corresponding to violated constraints. This
        %           would result in more small components.
        % NOTE: (2) and (3) cannot retain the feasibility of the solution, since it ignore
        % small violations.
        % To support post-processing, we return the index of violated-constriants.
        %
        % NOTE: due to rounding to zero, the capacity constraints will not be violated.
        function [tf, vars] = postProcessing(this)
            global DEBUG INFO;
            var_x = this.Variables.x;
            var_z = this.Variables.z;
            tol_zero = this.Parent.options.NonzeroTolerance;
            var_x(var_x<tol_zero*max(var_x)) = 0;
            var_z(var_z<tol_zero*max(var_z)) = 0;
            A1 = this.As_res(:,1:this.NumberPaths);
            A2 = this.As_res(:,(this.NumberPaths+1):end);
            v1 = A1*var_x;
            v2 = A2*var_z;
            index_violate = find(v1+v2>0);      % re = v1+v2
            if ~isempty(index_violate)
                message = sprintf('%s: Maximal violation is %.4f.', calledby, full(max(v1+v2)));
                if ~isempty(DEBUG) && DEBUG
                    warning(message); %#ok<SPWRN>
                elseif ~isempty(INFO) && INFO
                    cprintf('SystemCommands', '%s\n', message);
                end
                
                NP = this.NumberPaths;
                NV = this.NumberVNFs;
                post_process = this.Parent.options.PostProcessing;
                switch post_process
                    case 'round'
                        b_violate = false(size(v1));
                        b_violate(index_violate) = true;
                        pid = 1:NP:((NV-1)*NP+1);
                        for i = 1:NP
                            p_violate = find(b_violate(pid));
                            if ~isempty(p_violate)
                                t = v1(pid(p_violate))./abs(v2(pid(p_violate)));
                                var_x(i) = var_x(i) ./ max(t);
                            end
                            pid = pid + 1;
                        end
                        % The rounding error may still leads to postive error, therefore
                        % we set a relatively small tolerance, i.e. 1e-10. 
                        assert(this.checkFeasible([var_x; var_z], ...
                            struct('ConstraintTolerance', 1e-10)), ...
                            '[Post Processing]: failed infeasible solution.');
                    case {'drop','recover'}
                        %% Ignore the tiny components
                        % Use the differentiation of |v1-v2|/|v1+v2| is more flexible than use only
                        % |v1-v2|.
                        %   # if |v1-v2| is remarkable, both form can eaily identify the difference;
                        %   # If both v1 and v2 is large while |v1-v2| is small, the former can
                        %     identify the small value regardless of the tolerance, while the later
                        %     relies on the tolerance;
                        %   # Otherwise, if both v1 and v2 is small (most small value should have been
                        %     rounding in the former process), then both form rely on tolerance.
                        %
                        %% ISSUE
                        % |tol_zero| should not be too large compared to the absolute
                        % value of variables. It leads to inaccurate solution (the
                        % violation still exists), which cannot pass the feasible assertion.
                        % Therfore, *'drop' is not recommend*.  So |tol_zero| should be
                        % smaller, so that all violated ones can be
                        % identified.
                        %   # If we use the 'drop' option, some large components of the
                        %     solution might be discarded (even if the absolute value of
                        %     residual is large),  which degenerate the solution to the
                        %     sub-optmal.
                        %   # If we use the 'recover' option, then some trivial components
                        %     will be recovered. However the quality of solution declines.
                        %
                        % Since the |tol_zero| might lead to miss some violated components, we
                        % do not perform feasible assertion after processing. We Since the
                        % constraint tolerance in optimization is 10^-6 by default, the
                        % residual error's magnitude may be higher than 10^-6.
                        % Therefore, we have to set |tol_zero| to 10^-3~10^-4 (tuned according
                        % to the real data).
                        re = v1 + v2;
                        re(re<tol_zero) = 0;
                        mid_v = (abs(v1)+abs(v2))/2;
                        nz_index = find(mid_v~=0);
                        tol_con = this.Parent.options.ConstraintTolerance;
                        b_violate = (re(nz_index)./mean(mid_v(nz_index))) > tol_con; 
                        % => mean(mid_v(nz_index)) | max(mid_v(nz_index))
                        vi_idx = nz_index(b_violate);
                        pidx = mod(vi_idx-1, NP)+1;
                        if strcmpi(post_process, 'recover')
                            fidx = ceil(vi_idx/NV);
                            for i = 1:length(vi_idx)
                                zidx = (1:this.NumberDataCenters)+ ...
                                    ((fidx(i)-1)*this.NumberDataCenters*NP + ...
                                    (pidx(i)-1)*this.NumberDataCenters);
                                var_z(zidx) = this.temp_vars.z(zidx);
                            end
                        else
                            pidx = unique(pidx);
                            var_x(pidx) = 0;
                        end
                    otherwise
                        error('%s: invalid processing option (%s).', calledby, post_process);
                end
            end
            if nargout >= 2
                vars = [var_x; var_z];
            end
            this.Variables.x = var_x;
            this.Variables.z = var_z;
            tf = true;
        end
        
        function temp_vars = get_temp_variables(this)
            temp_vars = [this.temp_vars.x; this.temp_vars.z];
        end

        %%
        % <getLinkLoad> and <getNodeLoad> are only used inside the <Slice> class. The
        % substrate network cares how much resources are ocuppied by the slice, that is
        % the capacity of resources. See also <getLinkCapacity> and <getNodeCpacity>.
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
            v_n = this.Hrep*node_vars;
            % node_load = sum(f, node_vars(:,:,f).*I_node_path).
            %             NC = this.NumberDataCenters;
            %             NP = this.NumberPaths;
            %             v_n = zeros(NC,1);
            %             np = NC * NP;
            %             z_index = 1:np;
            %             for i = 1:this.NumberVNFs
            %                 node_vars_fi = reshape(node_vars(z_index), NC, NP);
            %                 v_n = v_n + sum(this.I_node_path.*node_vars_fi,2);
            %                 z_index = z_index + np;
            %             end
            %             assert(isempty(find(abs(this.Hrep*node_vars-v_n)>10^-6,1)), 'error: unequal node load');
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

    end
    
    methods (Access = {?Slice, ?CloudNetwork, ?SliceFlowEventDispatcher})
        [utility, node_load, link_load] = priceOptimalFlowRate(this, x0, new_opts);
    end
    methods (Access = {?Slice, ?CloudNetwork})
        [profit,cost] = optimalFlowRate( this, new_opts );
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
%      [profit, grad]= fcnProfit(vars, slice, options)
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