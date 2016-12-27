%%  Network Slice

%%
classdef Slice < handle
    
    properties (SetAccess = {?PhysicalNetwork})
        Identifier;     %
        Type;           % Type from slice template
        Topology;       % topology information of the slice
        VirtualLinks;   % information of virtual links in the slice
        VirtualNodes;   % information of virtual nodes in the slice
        FlowTable;      % flow information in the slice
        VNFList;        % List of virtual network functions in the slice
        Variables;      % [x,z]
        
        I_edge_path;    % Edge-Path Incidence Matrix
        I_node_path;    % Node-Path Incidence Matrix
        I_flow_path;    % Flow-Path Incidence Matrix
        
        NumberVirtualNodes;    % Number of nodes in the slice
        NumberVirtualLinks;    % Number of edges in the slice
        NumberFlows;    % Number of flows in the slice
        NumberPaths;    % Number of paths in the slice
        NumberVNFs;     % Number of virtual network functions in the slice
        
        I_path_function;
        alpha_f;
    end
    
    methods
        function this = Slice(slice_data)
            if nargin == 0
                return;
            elseif isa(slice_data, 'Slice')
                copy(this, slice_data);
            end
            this.Topology = DirectedGraph(slice_data.adjacent);
            if isfield(slice_data, 'Type')
                this.Type = slice_data.Type;
            end
            this.VirtualLinks = array2table(slice_data.link_map_S2P, ...
                'VariableNames', {'PhysicalLink'});
            if isfield(slice_data, 'link_price')
                this.VirtualLinks.Price = slice_data.link_price;
            end
            if isfield(slice_data, 'link_capacity')
                this.VirtualLinks.Capacity = slice_data.link_capacity;
            end
            this.VirtualLinks.Load = zeros(this.NumberVirtualLinks, 1);
            this.PhysicalLinkMap = array2table(slice_data.link_map_P2S,...
                'VariableNames', {'VirtualLink'});
            this.VirtualNodes = array2table(slice_data.node_map_S2P,...
                'VariableNames', {'PhysicalNode'});
            if isfield(slice_data, 'node_price')
                this.VirtualNodes.Price = slice_data.node_price;
            end
            if isfield(slice_data, 'node_capacity')
                this.VirtualNodes.Capacity = slice_data.node_capacity;
            end
            this.VirtualNodes.Load = zeros(this.NumberVirtualNodes, 1);
            this.PhyscialNodeMap = array2table(slice_data.node_map_P2S,...
                'VariableNames', {'VirtualNode'});
            this.FlowTable = slice_data.flow_table;
            this.VNFList = slice_data.VNFList;
            if isfield(slice_data,'Weight')
                this.weight = slice_data.Weight;
            end
            if isfield(slice_data,'Identifier')
                this.Identifier = slice_data.Identifier;
            end
            
            this.Parent = slice_data.parent;
            if isfield(slice_data, 'alpha_f')
                this.initializeState(slice_data.alpha_f);
                this.alpha_f = slice_data.alpha_f;
            else
                this.initializeState;
            end
        end
        
        function n = get.NumberVirtualNodes(this)
            n = this.Topology.NumberNodes;
        end
        
        function m = get.NumberVirtualLinks(this)
            m = this.Topology.NumberEdges;
        end
        
        function f = get.NumberFlows(this)
            f = height(this.FlowTable);
        end
        
        function p = get.NumberPaths(this)
            p = 0;
            for i=1:this.NumberFlows
                p = p + this.FlowTable.Paths(i).Width;
            end
        end
        
        function n = get.NumberVNFs(this)
            n = length(this.VNFList);
        end
        %         function p = get.Parent(this)
        %             p = this.Parent;
        %         end
        function n = get.num_vars(this)
            n = (this.NumberVNFs*this.NumberVirtualNodes+1)*this.NumberPaths;
        end
        
        function n = get.num_lcon_res(this)
            n = size(this.As_res,1);
        end
        
        function pid = getLocalPathId(this, path)
            pid = path.id - this.FlowTable.Paths(1).paths{1}.id + 1;
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
        
        function vn = getNodeLoad(this, node_vars)
            % retrive the node load of the slice, given the node variables.
            if nargin == 1
                node_vars = this.Variables.z;
            end
            %%
            % |node_vars| is index by |(node,path,function)|.
            % node_load = sum(f, node_vars(:,:,f).*I_node_path).
            NN = this.NumberVirtualNodes;
            NP = this.NumberPaths;
            NV = this.NumberVNFs;
            vn = zeros(NN,1);
            np = NN * NP;
            z_index = 1:np;
            for i = 1:NV
                node_vars_fi = reshape(node_vars(z_index), NN, NP);
                vn = vn + sum(this.I_node_path.*node_vars_fi,2);
                z_index = z_index + np;
            end
            
            %% Alternative way to compute the node load.
            %             col_index = (1:NN:((NP-1)*NN+1))';
            %             col_index = repmat(col_index, 1, NV);
            %             for c = 2:NV
            %                 col_index(:,c) = col_index(:,c-1) + NN*NP;
            %             end
            %             col_index = col_index(:);
            %             As = zeros(NN, NP*NV);
            %             for row_index = 1:NN
            %                 As(row_index, col_index) = repmat(this.I_node_path(row_index,:),1, NV);
            %                 col_index = col_index + 1;
            %             end
            %             vn = As*node_vars;
        end
   
        function pc = getPathCost(this, lambda_e, lambda_n)
            link_uc = this.Parent.getLinkField('UnitCost', this.VirtualLinks.PhysicalLink); % the virtual links's unit cost
            node_uc = this.Parent.getNodeField('UnitCost', this.VirtualNodes.PhysicalNode); % the virtual nodes's unit cost
            delta = this.Parent.delta;
            epsilon = this.Parent.unitStaticNodeCost;
            N = this.Parent.NumberNodes;
            phis_n = epsilon*delta*(N-1)/this.Parent.totalNodeCapacity;
            phis_l = epsilon*(1-delta)*(N-1)/this.Parent.totalLinkCapacity;
            pc = cell(this.NumberFlows, 1);
            p = 1;
            for i = 1:this.NumberFlows
                num_paths = this.FlowTable{i,'Paths'}.Width;
                pc{i} = zeros(num_paths, 1);
                for j = 1:num_paths
                    eid = this.I_edge_path(:,p)~=0;
                    pc{i}(j) = pc{i}(j) + sum(link_uc(eid)+lambda_e(eid)) + phis_l*nnz(eid);
                    nid = this.I_node_path(:,p)~=0;
                    pc{i}(j) = pc{i}(j) + min(node_uc(nid)+lambda_n(nid)) + phis_n;
                    p = p + 1;
                end
            end
        end
        
        %% Slice Cost
        % When compute the static cost, the Capacity of all physical nodes and links is
        % included.
        % |epsilon/pn.NumberSlices| is a constant. To keep consistence with other methods,
        % this part should not be ignored. See also _getNetworkCost_ and _getStaticCost_.
        % The calculation is not absolutely precise, since it cannot be decide that the
        % static cost should be arributed to which slices.
        %
        % When calculate network cost as a single slice, |pn.NumberSlices| is not
        % applicable. So this method can not be used, use _getNetworkCost_ instead. 
        function sc = getSliceCost(this, node_load, link_load, model)
            if nargin <= 1 || isempty(node_load)
                node_load = this.VirtualNodes.Load;
            end
            if nargin <= 2 || isempty(link_load)
                link_load = this.VirtualLinks.Load;
            end
            if nargin <=3
                model = 'Approximate';
            end
            
            pn = this.Parent;
            link_uc = pn.getLinkField('UnitCost', this.VirtualLinks.PhysicalLink); % the virtual links's unit cost
            node_uc = pn.getNodeField('UnitCost', this.VirtualNodes.PhysicalNode); % the virtual nodes's unit cost
            epsilon = pn.unitStaticNodeCost;
            N = pn.NumberNodes;
            phis_n = epsilon*pn.delta*(N-1)/pn.totalNodeCapacity;
            phis_l = epsilon*(1-pn.delta)*(N-1)/pn.totalLinkCapacity;
            
            if strcmp(model, 'Approximate')
                %% TODO, when the slice has no invalid parent (Physical Network)
                % pn.NumberSlices is invalid.
                sc = dot(link_uc, link_load) + dot(node_uc, node_load) ...
                    + (phis_n*sum(node_load)+phis_l*sum(link_load)+epsilon/pn.NumberSlices);
            elseif strcmp(model, 'Accurate')
                %% Accurate Model
                % A slice cannot decide how to devide the static cost between slices.
                % one method is devision by useage, using the physic network's load data
                % of each slice. So this function only calculate the dynamic part of the
                % cost, and the slices should further calculate the static cost outside
                % this method. 
                sc = dot(link_uc, link_load) + dot(node_uc, node_load);
            end
        end
        
        %% Get the node cost on a path.
        function [nc, nid] = getPathNodeCost(this, pid, lambda_n)
           node_uc = this.Parent.getNodeField('UnitCost', this.VirtualNodes.PhysicalNode);
           epsilon = this.Parent.unitStaticNodeCost;
           N = this.Parent.NumberNodes;
           phis_n = epsilon*this.Parent.delta*(N-1)/this.Parent.totalNodeCapacity;
           nid = find(this.I_node_path(:,pid)~=0);
           nc = node_uc(nid) + lambda_n(nid) + phis_n;           
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
                    warning('checkFeasible: Maximal violation is %.4f.\n', max(ub));
                end
            end
        end
        
        
        [utility, node_load, link_load] = priceOptimalFlowRate(this, x0, options);
        x = priceOptimalFlowRateCompact(this, x0);
        [fval, node_load, link_load] = subproblemNetSocialWelfare( this, lambda );
        [fval, node_load, link_load] = subproblemNetSocialWelfare2( this, lambda );
    end
    
    % Specify the properties that can only be modified by Physcial Network directly
    properties (SetAccess = {?PhysicalNetwork})
        Parent;
    end
    
    properties (GetAccess = {?PhysicalNetwork})
        PhysicalLinkMap;
        PhyscialNodeMap;    %
        path_owner;         % the associate flow of path;
        num_vars;           % number of optimization variables in the problem
        num_lcon_res;       % number of linear constraints
        As_res;             % coefficient matrix of the processing constraints
        x_path;             % allocated bandwidth of paths
        z_npf;              % allocated resource on nodes
        x0;                 % start point
        I_active_variable;  % indicator of active variables
        weight;             % weight for the slice user's utility.
        link_load;          % temporary results of link load.
        node_load;          % temporary results of node load.
        flow_rate;          % temporary results of flow rate.
    end
    
    
    methods (Access = public, Static)
        
        slice_template = loadSliceTemplate(index);

        % Objective function and gradient
        [profit, grad]= fcnProfit(var_x, S);
        [profit, grad]= fcnProfitCompact(act_var_x, S);
        [fval,  grad] = subproblemObjective(var_x, lambda, S);
        [profit, grad] = fcnNetProfit(x_vars, S, options);

        % Hessian matrix
        hess = fcnHessian(var_x, ~, S);
        hess = fcnHessianCompact(act_var_x, ~, S);
        
    end
    
    methods (Access = private)
        function initializeState(this, alpha_f, options)
            if nargin <= 2
                options = 1;
            end
            NN = this.NumberVirtualNodes;
            NP = this.NumberPaths;
            NL = this.NumberVirtualLinks;
            NF = this.NumberFlows;
            NV = this.NumberVNFs;
            if options == 1
            this.I_node_path = sparse(NN, NP);
            this.I_edge_path = sparse(NL, NP);
            this.I_flow_path = sparse(NF, NP);
            this.path_owner = zeros(NP,1);
            pid = 0;
            for fid=1:NF
                path_list = this.FlowTable{fid,{'Paths'}};
                for j = 1:path_list.Width
                    pid = pid + 1;
                    this.I_flow_path(fid,pid) = 1;
                    this.path_owner(pid) = fid;
                    path = path_list.paths{j};
                    for k = 1:(path.Length-1)
                        e = path.Link(k);
                        eid = this.Topology.IndexEdge(e(1),e(2));
                        this.I_edge_path(eid, pid) = 1;
                        this.I_node_path(e(1), pid) = 1;
                    end
                    this.I_node_path(e(2), pid) = 1;
                end
            end
            end
            
            %% Linear constraint without consdiering the bound constraint.
            % |As| takes the following form
            %
            % $$\left[ \begin{array}{ccccc}
            %   I_1     & H_s &      &          &       \\
            %   I_2     &     & H_s  &          &       \\
            %   \vdots  &     &      &  \ddots  &       \\
            %   I_F     &     &      &          & H_s
            % \end{array} \right]
            % \left[ \begin{array}{c}x\\z_1\\z_2\\z_i\\z_F\end{array}\right] $$
            %
            % where
            %
            % $$I_f = \left[ \begin{array}{cccc}
            % \alpha_f &          &        &         \\
            %          & \alpha_f &        &         \\
            %          &          & \ddots &         \\
            %          &          &        & \alpha_f
            % \end{array} \right],$$
            % $$H_s = \left[ \begin{array}{cccccccccc}
            % -h_{1,1} & \cdots & -h_{N,1} &          &        &          &        &         &        & \\
            %          &        &          & -h_{1,2} & \cdots & -h_{N,2} &        &         &        & \\
            %          &        &          &          &        &          & \ddots &         &        & \\
            %          &        &          &          &        &          &        &-h_{1,P} & \cdots & -h_{N,P}
            % \end{array} \right],$$
            % $$ x = \left[\begin{array}{c}x_1\\x_2\\ \vdots\\x_P\end{array}\right]$$
            % $$ z_f = \left[\begin{array}{c}z_{1,1,1f}\\ \vdots\\ z_{N,1,f}\\z_{1,2,f}\\
            %   \vdots\\ z_{N,2,f}\\ \vdots\\z_{N,P,F}\end{array}\right]$$
            %
            % According to the martix formulation, the number of non-zero elements in |As|
            % is equal to |F*(P+nnz(Hs))|, where |nnz(Hs)| is equal to the number of
            % nonzero elements in |I_node_path|.
            nnz_As = NV*(NP+nnz(this.I_node_path));
            num_lcon = this.NumberPaths*NV;
            this.As_res = spalloc(num_lcon, this.num_vars, nnz_As);
            row_index = 1:NP;
            col_index = NP+(1:NN);
            for f = 1:NV
                if nargin >= 2 && ~isempty(alpha_f) % used when treat all VNFs as one function.
                    for p = 1:NP
                        this.As_res(p,p) = alpha_f(this.path_owner(p));
                    end
                else
                    af = this.Parent.VNFTable{this.VNFList(f),{'ProcessEfficiency'}};
                    this.As_res(row_index,1:NP) = af * eyesparse(NP);
                end
                for p = 1:NP
                    this.As_res(row_index(p), col_index) = -this.I_node_path(:,p);
                    col_index = col_index + NN;
                end
                row_index = row_index + NP;
            end
        end
        
        function copy(this, slice)
            this.Identifier = slice.Identifier;     
            this.Topology = DirectedGraph(slice.Topology);
            this.VirtualLinks = slice.VirtualLinks;   
            this.VirtualNodes = slice.VirtualNodes;   
            this.VNFList = slice.VNFList;        
            this.Variables = slice.Variables;      
            this.I_edge_path = slice.I_edge_path;
            this.I_node_path = slice.I_node_path;
            this.I_flow_path = slice.I_flow_path;    
            this.I_path_function = slice.I_path_function;
            this.alpha_f = slice.alpha_f;
            
            this.Parent = slice.Parent;
            
            this.PhysicalLinkMap = slice.PhysicalLinkMap;
            this.PhyscialNodeMap = slice.PhyscialNodeMap;
            this.path_owner = slice.path_owner;         
            this.As_res = slice.As_res;             
            this.x_path = slice.x_path;             
            this.z_npf = slice.z_npf;              
            this.x0 = slice.x0;                 
            this.I_active_variable = slice.I_active_variable;  
            this.weight = slice.weight;             
            this.link_load = slice.link_load;          
            this.node_load = slice.node_load;          
            this.flow_rate = slice.flow_rate;   
            
            this.FlowTable = slice.FlowTable;
            for f = 1:height(slice.FlowTable)
                % path_list is handle object, is should be copyed to the new table.
                path_list = PathList(slice.FlowTable{f,'Paths'});
                this.FlowTable{f,'Paths'} = path_list;
            end            
        end
    end
   
end

%% Properties
% * *Topology*: Normally, the network slice will not run shrtest path algorithm, so the
% absolute value of the adjacent matrix of Topology does not matter. On the other hand,
% the link and node capacity of the slice is also not determined until the substrate
% network allocate the resource to the slice.
% * *VirtualLinks* : fields include _PhysicalLink_, _Price_, _Load_.
% * *VitrualNodeMap* : fields include _PhysicalNode_, _Price_, _Load_.
%
% * *NumberVirtualNodes* |get|
%
% * *NumberVirtualLinks* |get|
%
% * *NumberFlows* |get|
%
% * *NumberPaths* |get|

%% Methods
% * *getLocalPathId* : find the path's local identifier.
%
%      pid = getLocalPathId(slice, path)
%
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
% * *fcnProfitCompact* |static| : Evalute the objective function and gradient.
%
%      [profit, grad]= fcnProfitCompact(act_var_x, S)
%
% Only active independent variables are passed into the objective function.
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
%
% * *fcnHessianCompact* |static| : Compact form of Hessian matrix of the Largrangian.
%
%      hess = fcnHessianCompact(act_var_x, ~, S)