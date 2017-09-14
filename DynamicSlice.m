classdef DynamicSlice < Slice & EventSender & EventReceiver
    %DynamicSlice Event-driven to dynamic configure slice.
    %   Detailed explanation goes here

    properties
        FlowArrivalRate;
        FlowServiceInterval;
    end
    properties(Access={?DynamicSlice,?DynamicNetwork})
        b_ondepart = false;
        Vnf;            % Capaciy of VNF instances on each node, this is configured by inter-slicing
        Hs;
        % As_res;       % has been defined in <Slice>.
        old_variables;  % last one configuration
        changed_index;
        %%% define reconfiguraton cost for link and node variables.
        % **flow reassginement cost*, i.e., the flow-path reconfiguration cost, which is
        % denpendent on length of paths.
        % **node reassginement cost*,
        %
        % If a resource price is p, we assume that the reconfigruation cost
        % coefficient of this resource is ¦Èp. If we adopt a varying pricing, like
        % p=av+b, then the coeffcient is ¦È(av+b).
        edge_reconfig_cost;
        node_reconfig_cost;
    end
    
    events 
        AddFlowSucceed;
        AddFlowFailed;
        RemoveFlowSucceed;
        RemoveFlowFailed;          % NOT used.
    end
    methods
        function eventhandler(this, source, eventData) %#ok<INUSL>
            global DEBUG; %#ok<NUSED>
            %             target = eventData.targets;
            % where target should be the owner of arrving flows
            sl = eventData.slice;    
            if sl ~= this  % filter message
                return;
            end
            switch eventData.EventName
                case 'FlowArrive'
                    fidx = this.OnAddingFlow(eventData.flow);
                    % notify slice
                    if isempty(fidx)
                        notify(this, 'AddFlowFailed');
                    else
                        ev = eventData.event;
                        eventData = FlowEventData(ev, sl, fidx);
                        notify(this, 'AddFlowSucceed', eventData);
                    end
                case 'FlowDepart'
                    % notify slice
                    flow_id = eventData.flow;
                    eventData = FlowEventData(eventData.event, sl, flow_id);
                    this.OnRemovingFlow(flow_id);
                    notify(this, 'RemoveFlowSucceed', eventData);
                otherwise
                    error('error: cannot hand event %s.', eventData.EventName);
            end
        end
        %% TODO: only update the colums arriving/removing 
        % since the resource changes after each slice configuration, the matrix needs be
        % computed each time. Therefore, we define the access functio to evalue it.
        % This decoupt the computation of incident matrix and As_res, in
        % <initializeState>. 
        %         function As = getAs_res(this)
        %         end
        
        function Hs = get.Hs(this)
            % For the node capacity constraint for VNF |f|, the coeffiect matrix Hs is
            % filled block-by-block. In the |k|th block (NC*NC), we put the |k|th column
            % of H_np to the diagnal of the block. Then |H_np| is duplicated into a larger
            % diagnal corresponding to all VNFs.
            %
            % The capacity of VNF instance Vnf, has been calculate after allocating
            % resource to the slice. Some component of Vnf might be zero, since VNF f may
            % not be located at node n.
            global DEBUG; %#ok<NUSED>
            NC = this.NumberDataCenters;
            NP = this.NumberPaths;
            Hs_np = spalloc(NC, NC*NP, nnz(this.I_node_path));
            col_index = 1:NC;
            for p = 1:NP
                Hs_np(1:NC,col_index) = diag(this.I_node_path(:,p));
                col_index = col_index + NC;
            end
            Hs = block_diag(Hs_np, this.NumberVNFs);
        end
        
        %% fast slice reconfiguration when flow arriving and depaturing
        % * *flow*: flow table entries.
        %
        % *Update state*: when single flow arrive or depart, update the state record,
        % including: |I_node_path, I_edge_path, I_flow_path, path_owner| and |As|.
        % Since only add/remove entry for only one flow, this operation is more efficient
        % than <initializeState>.
        function fidx = OnAddingFlow(this, flow)  
            % temporarily allocated flow id.
            global DEBUG; %#ok<NUSED>
            num_exist_flows = height(this.FlowTable);
            num_new_flows = height(flow);
            num_exist_paths = this.NumberPaths;
            fidx = (1:num_new_flows) + num_exist_flows;
            flow.Identifier = (1:num_new_flows) + this.FlowTable{end,'Identifier'}; % temporary identifier
            if num_exist_flows == 0
                pid = 0;
            else
                pid = this.FlowTable{end, 'Paths'}.paths{end}.local_id;
            end
            
            old_state = this.get_state; % if fail to adding flow, recover flow table
            
            this.FlowTable = [this.FlowTable; flow];
            
            changed_path_index((num_exist_paths+1):this.NumberPaths) = true;
            for fid = fidx
                for p = 1:length(this.FlowTable{fid, 'Paths'})
                    pid = pid + 1;
                    this.I_flow_path(fid, pid) = 1;
                    this.path_owner(pid) = fid;
                    path = this.FlowTable{fid, 'Paths'}.paths{p};
                    path.local_id = pid;
                    for k = 1:(path.Length-1)
                        e = path.Link(k);
                        eid = this.Topology.IndexEdge(e(1),e(2));
                        this.I_edge_path(eid, pid) = 1;
                        dc_index = this.VirtualNodes{e(1),'DataCenter'};
                        if dc_index~=0
                            this.I_node_path(dc_index, pid) = 1;
                        end
                    end
                    dc_index = this.VirtualNodes{e(2),'DataCenter'}; % last node
                    if dc_index~=0
                        this.I_node_path(dc_index, pid) = 1;
                    end
                end
            end
            this.getAs_res;        % update As_res

            % used for all methods to calculate the reconfiguration cost. and used for
            % fast reconfiguration during optimization.
            this.changed_index.x = sparse(changed_path_index);
            this.identify_change(this.changed_index.x);
            this.old_variables.x = zeros(this.NumberPaths,1);
            this.old_variables.z = zeros(this.num_varz,1);
            this.old_variables.x(~this.changed_index.x) = old_state.variables.x;
            this.old_variables.z(~this.changed_index.z) = old_state.variables.z;

            %% TODO, designed the output structure for multiple slices.
            global g_results event_num;
            options = this.Parent.options;
            switch options.Method
                case 'reconfig'
                    % provide 'method' and 'model' to customize the <optimalFlowRate>
                    options.Method = 'slice';
                    options.Model = 'FixedCost'; 
                    profit = this.optimalFlowRate(options);
                    g_results.profit(event_num,1) = ...
                        profit - this.getSliceCost('quadratic-price');
                    g_results.solution(event_num,1) = this.Variables;
                    [   g_results.cost(event_num,1), ...
                        g_results.num_reconfig(event_num,1),...
                        g_results.rat_reconfig(event_num,1),...
                        g_results.num_vars(event_num,1)] ...
                        = this.get_reconfig_cost;
                case 'fastconfig'
                    profit = this.fastReconfigure('add');
                    g_results.profit(event_num,1) = ...
                        profit - this.getSliceCost('quadratic-price');
                    g_results.solution(event_num,1) = this.Variables;
                    [   g_results.cost(event_num,1),...
                        g_results.num_reconfig(event_num,1),...
                        g_results.rat_reconfig(event_num,1),...
                        g_results.num_vars(event_num,1)]...
                        = this.get_reconfig_cost;
                case 'fastconfig2'
                otherwise
                    error('NetworkSlicing:UnsupportedMethod', ...
                        'error: unsupported method (%s) for network slicing.', ...
                        options.Method) ;
            end
            g_results.num_flows(event_num,1) = this.NumberFlows;
            % if failed
            % fidx = [];
        end
        
        function tf = OnRemovingFlow(this, fid)
            global DEBUG; %#ok<NUSED>
            fidx = this.FlowTable.Identifier == fid;
            if isempty(find(fidx==true,1))
                warning('flow (identifier %d) not in slice (identifer %d).',...
                    fid, this.Identifier);
            end
            old_state = this.get_state; % if fail to adding flow, recover flow table
            changed_path_index = false(this.NumberPaths,1);
            flow_index = find(fidx);
            for fid = flow_index
                for p = 1:length(this.FlowTable{fid, 'Paths'})
                    pid = this.FlowTable{fid, 'Paths'}.paths{p}.local_id;
                    changed_path_index(pid) = true;
                end
            end
            % no need to compare the changed variables, its difference is constant
            % (0-x0)
            % used for all methods to calculate the reconfiguration cost.
            this.changed_index.x = sparse(changed_path_index);
            this.identify_change(this.changed_index.x); % must be called before flows are removed from the table
            this.old_variables.x = old_state.variables.x(~this.changed_index.x);
            this.old_variables.z = old_state.variables.z(~this.changed_index.z);
            % re-allocate local path ids;
            %             num_rm_paths = nnz(changed_path_index);
            %             for fid = fidx(end):this.NumberFlows
            %                 for p = 1:length(this.FlowTable{fid, 'Paths'})
            %                     path = this.FlowTable{fid, 'Paths'}.paths{p};
            %                     path.local_id = path.local_id - num_rm_paths;
            %                 end
            %             end
            % update flow table and local path id.
            pid = this.FlowTable{flow_index(1), 'Paths'}.paths{1}.local_id-1;
            fidt = flow_index(1):(this.NumberFlows-length(flow_index));
            this.FlowTable(fidx,:) = []; 
            this.path_owner(changed_path_index) = [];
            for fid = fidt
                path_list = this.FlowTable{fid,{'Paths'}};
                for j = 1:path_list.Width
                    pid = pid + 1;
                    path_list.paths{j}.local_id = pid;
                    this.path_owner(pid) = fid;
                end
            end
            this.I_edge_path(:, changed_path_index) = [];
            this.I_node_path(:, changed_path_index) = [];
            this.I_flow_path(:, changed_path_index) = [];
            this.I_flow_path(fidx, :) = [];
            this.getAs_res;

            %% TODO
            % perform resource allocation
            global g_results event_num;
            options = this.Parent.options;
            switch options.Method
                case 'reconfig'
                    options.Method = 'slice';
                    options.Model = 'FixedCost';
                    profit = this.optimalFlowRate(options);
                    g_results.profit(event_num,1) = ...
                        profit - this.getSliceCost('quadratic-price');
                    g_results.solution(event_num,1) = this.Variables;
                    [   g_results.cost(event_num,1),...
                        g_results.num_reconfig(event_num,1),...
                        g_results.rat_reconfig(event_num,1),...
                        g_results.num_vars(event_num,1)]...
                        = this.get_reconfig_cost;
                case 'fastconfig'
                    profit = this.fastReconfigure('remove');
                    g_results.profit(event_num,1) = profit - ...
                        this.getSliceCost('quadratic-price');
                    g_results.solution(event_num,1) = this.Variables;
                    [   g_results.cost(event_num,1),...
                        g_results.num_reconfig(event_num,1),...
                        g_results.rat_reconfig(event_num,1),...
                        g_results.num_vars(event_num,1)]...
                        = this.get_reconfig_cost;
                case 'fastconfig2'
                otherwise
                   error('NetworkSlicing:UnsupportedMethod', ...
                       'error: unsupported method (%s) for network slicing.', ...
                       options.Method);
            end
            g_results.num_flows(event_num,1) = this.NumberFlows;
            %             flowtable = this.FlowTable(~rmidx,:);
            tf = true;
        end
        
        %%%
        % override <Slice.getSliceCost>.
        %   cost = getSliceCost(this, pricing_policy, reconfig_cost_model)
        % |reconfig_cost_model|: reconfigure cost model, including: |'constant'|,
        % |'linear'|, and |'none'|. If |'none'| is adopted, then the reconfiguration cost
        % is not counted.
        function cost = getSliceCost(this, pricing_policy, reconfig_cost_model)
            if nargin<=1 || isempty(pricing_policy)
                pricing_policy = 'quadratic-price';
            end
            if nargin <= 2 || isempty(reconfig_cost_model)
                reconfig_cost_model = 'linear';
            end
            
            link_price = this.VirtualLinks.Price;
            link_load = this.VirtualLinks.Capacity;
            node_price = this.VirtualDataCenters.Price;
            node_load = this.VirtualDataCenters.Capacity;
            if strcmpi(pricing_policy, 'quadratic-price')
                link_payment = this.fcnLinkPricing(link_price, link_load);
                node_payment = this.fcnNodePricing(node_price, node_load);
                cost = link_payment + node_payment;
            else
                cost = dot(link_price, link_load) + dot(node_price, node_load);
            end
            if ~strcmpi(reconfig_cost_model, 'none')
                cost = cost + this.get_reconfig_cost(reconfig_cost_model);
            end
        end
        
        %% Find the index of variables for the newly added/removed flow
        % assume the newly added/removed flow index is u, the the index of
        % corrsponding flow and VNF allocation variables is
        % [x_u, z_npf], for all p in p_u.
        %
        % NOTE: This function should be called after adding new flow, or before removing
        % departing flow.
        function identify_change(this, changed_path_index)  
            global DEBUG; %#ok<NUSED>
            base_z_index = false(this.NumberDataCenters, this.NumberPaths);
            base_z_index(:,changed_path_index) = true;
            base_zidx = find(base_z_index);
            this.changed_index.z = spalloc(this.num_varz, 1, ...
                nnz(base_z_index)*this.NumberVNFs);
            id_offset = 0;
            for f = 1:this.NumberVNFs
                this.changed_index.z(id_offset+base_zidx) = true;
                id_offset = id_offset + numel(base_z_index);
            end
        end
        
        %% TODO: when consifer the reconfiguration cost, the path length should be consdiered with path cost 
        function profit = fastReconfigure(this, action)
            NL = this.NumberVirtualLinks;
            NN = this.NumberDataCenters;
            NV = this.NumberVNFs;
            NP = this.NumberPaths;

            %%% Formulate input for convex optimization (fmincon).
            % The problem has multiple inequalities, and the lowerbounds for variables.
            global DEBUG; %#ok<NUSED>
            H = this.Hs;
            As_res = this.As_res;        % update As_res
            nnz_As = nnz(As_res) + nnz(H) + nnz(this.I_edge_path) + ... 
                 + 4*NP + 4*this.num_varz; % 2 I_x, 2 I_tx, 2 I_z, 2 I_tz
            num_lcon = this.num_lcon_res + NN*NV + NL + 2*NP + 2*this.num_varz;
            As = spalloc(num_lcon, 2*this.num_vars, nnz_As);
            As(1:this.num_lcon_res,1:this.num_vars) = As_res;
            row_offset = this.num_lcon_res;
            As(row_offset+(1:NN*NV), (1:this.num_varz)+NP) = H;
            row_offset = row_offset + NN*NV;
            As(row_offset+(1:NL), 1:NP) = this.I_edge_path;
            row_offset = row_offset + NL;
            As(row_offset+(1:NP), 1:NP) = eye(NP);
            As(row_offset+(1:NP), this.num_vars+(1:NP)) = -eye(NP);
            row_offset = row_offset + NP;
            As(row_offset+(1:NP), 1:NP) = -eye(NP);
            As(row_offset+(1:NP), this.num_vars+(1:NP)) = -eye(NP);
            row_offset = row_offset + NP;
            As(row_offset+(1:this.num_varz), (NP+1):this.num_vars) = eye(this.num_varz);
            As(row_offset+(1:this.num_varz), (this.num_vars+NP+1):end) = -eye(this.num_varz);
            row_offset = row_offset + this.num_varz;
            As(row_offset+(1:this.num_varz), (NP+1):this.num_vars) = -eye(this.num_varz);
            As(row_offset+(1:this.num_varz), (this.num_vars+NP+1):end) = -eye(this.num_varz);

            bs = [sparse(this.num_lcon_res,1);
                this.Vnf(:);
                this.VirtualLinks.Capacity;
                this.old_variables.x;       % which have been pre-processed, so it can be
                -this.old_variables.x;      % compared with the current states.
                this.old_variables.z;
                -this.old_variables.z];
            lbs = sparse(2*this.num_vars,1);
                
            var0 = [this.old_variables.x;
                this.old_variables.z;
                this.old_variables.x;
                this.old_variables.z];
            %% Reconfiguration Cost in Problem Formulation
            % comparing the old variables with new variables to decide the reconfiguration
            % cost.
            % Since the two vector has different number of elements, we should comparing
            % it accordingly, and set aside the variables of new arriving/departing flow.
            
            %% Perform optimization
            fmincon_opt = optimoptions(@fmincon);
            fmincon_opt.Algorithm = 'interior-point';
            fmincon_opt.HessianFcn = @(x,lambda)Slice.fcnHessian(x,lambda,this);
            fmincon_opt.SpecifyObjectiveGradient = true;
            [x, fval, exitflag] = fmincon(@(x)DynamicSlice.fcnFastConfigProfit(x,this), ...
                var0, As, bs, [], [], lbs, [], [], fmincon_opt);
            % x is a local solution to the problem when exitflag is positive.
            if exitflag == 0
                warning('reaching maximum number of iterations.');
            elseif exitflag < 0
                error('abnormal exit with flag %d.',exitflag);
            elseif exitflag ~= 1
                warning('(exitflag = %d) local optimal solution found.', exitflag);
            end
            options = this.Parent.options;
            options.Action = action;
            if ~this.checkFeasible(x, options)
                error('error: infeasible solution.');
            end
            %%%
            % The simplify process in <optimalFlowRate> is not needed, since the objective function will
            % force those variables to be zero.
            this.x_path = x(1:NP);
            this.z_npf = x((this.NumberPaths+1):this.num_vars);
            this.flow_rate = this.getFlowRate(this.x_path);
            if nargout == 1    % final results
                this.Variables.x = this.x_path;
                this.Variables.z = this.z_npf;
                %                 this.Variables.z(this.Variables.z<10^-3) = 0;
                %                 this.Variables.x(this.Variables.x<10^-3) = 0;
                tol_zero = 10^-4;
                this.Variables.x(this.Variables.x<tol_zero*max(this.Variables.x)) = 0;
                this.Variables.z(this.Variables.z<tol_zero*max(this.Variables.z)) = 0;
                options.Display = 'final';
                if ~this.checkFeasible([], options)
                    warning('optimalFlowRate: the rounding of variables with small quantity will make the solution infeasible.');
                end
                this.setPathBandwidth;
                this.FlowTable.Rate = this.getFlowRate;
                this.VirtualLinks.Load = this.getLinkLoad;
                this.VirtualDataCenters.Load = this.getNodeLoad;
                profit = -fval;
            end
        end  
        
        % optimize only on the subgraph, which carries the flows that intersection with
        % the arriving/departuring flow. this can significantly reduce the problem scale
        % when there is lots of flow, and the numbe of hops of the flow is small.
        %         function profit = fastReconfigureCompact(this, action)
        %         end
        
        function b = checkFeasible(this, vars, options)
            if nargin <= 1
                vars = [];
            end
            if nargin <= 2
                options = [];
            end
            if isempty(vars)
                b = checkFeasible@Slice(this, vars, options);
            else
                b = checkFeasible@Slice(this, vars(1:this.num_vars), options);
            end
            if b
                % TODO add other check coditions.
            end
        end
    end

    methods
        function this = DynamicSlice(slice_data)
            this@Slice(slice_data);
            if isfield(slice_data, 'Flow')
                if isfield(slice_data.Flow, 'ArrivalRate') && ...
                        slice_data.Flow.ArrivalRate>0 && slice_data.Flow.ArrivalRate<inf
                    this.FlowArrivalRate = slice_data.Flow.ArrivalRate;
                end
                if isfield(slice_data.Flow, 'ServiceInterval') && ...
                    slice_data.Flow.ServiceInterval>0 && slice_data.Flow.ServiceInterval<inf
                    this.FlowServiceInterval = slice_data.Flow.ServiceInterval;
                end
            end
        end
        
        function tf = isDynamicFlow(this)
            if isempty(this.FlowArrivalRate) || isempty(this.FlowServiceInterval)
                tf = false;
            else
                tf = true;
            end
        end
        %%%
        % *Create new flows*
        % Creating new flows in the slice could guarantee no extra node or link would be
        % needed. If we enable new flows from new locations, we should create the flow
        % in the network.
%         function ft = createflow(this, slice)
%             graph = slice.Topology;
%             switch slice.Options.FlowPattern
%                 case {FlowPattern.RandomSingleFlow, FlowPattern.RandomMultiFlow}
%                     ft = this.generateFlowTable(graph,
%                 case FlowPattern.RandomInterDataCenter
%                 case FlowPattern.RandomInterBaseStation
%                 case FlowPattern.RandomDataCenter2BaseStation
%                 otherwise
%                     error('error: cannot handle the flow pattern <%s>.', ...
%                         slice.Options.FlowPattern.char);
%             end
%         end
                
    end
    methods (Access = {?DynamicSlice, ?DynamicNetwork})
        function setVnfCapacity(this)
            znpf = reshape(this.Variables.z, this.NumberDataCenters, ...
                this.NumberPaths, this.NumberVNFs);
            znpf = znpf.* full(this.I_node_path);  % compatible arithmetic operation
            this.Vnf = reshape(sum(znpf,2), this.NumberDataCenters, this.NumberVNFs);
        end
        
        % Return: reconfiguration cost, number of reconfiguration, ration of
        % reconfigurations, and total number of variables.
        function [c,n,r,t] = get_reconfig_cost(this, model)
            if isempty(this.old_variables)
                c = 0;
                if nargout >=2
                    n = 0;
                end
                if nargout >= 3
                    r = 0;
                end
                if nargout >= 4
                    t = 0;
                end
                return;
            end
            if nargin <= 1
                model = 'linear'; % model = {linear|constant}
            end 
            % old_variables is valid after adding/removing flows.
            path_reconfig_cost = (this.I_edge_path)' * this.edge_reconfig_cost;
            diff_x = abs(this.Variables.x-this.old_variables.x);
            diff_z = abs(this.Variables.z-this.old_variables.z);
            mid_x = 1/2*(abs(this.Variables.x)+abs(this.old_variables.x));
            mid_z = 1/2*(abs(this.Variables.z)+abs(this.old_variables.z));
            nz_index_x = mid_x~=0;
            nz_index_z = mid_z~=0;
            norm_diff_z = zeros(size(diff_z)); 
            norm_diff_z(nz_index_z) = abs(diff_z(nz_index_z)./mid_z(nz_index_z));
            norm_diff_x = zeros(size(diff_x)); 
            norm_diff_x(nz_index_x) = abs(diff_x(nz_index_x)./mid_x(nz_index_x));  
            tol_vec = 10^-4;
            switch model
                case 'linear'
                    sum_tz = sum(reshape(diff_z, ...
                        this.NumberDataCenters, this.NumberPaths*this.NumberVNFs),2);
                    c = dot(diff_x, path_reconfig_cost)...
                        + dot(sum_tz, this.node_reconfig_cost);
                case 'constant'
                    sum_tz = sum(reshape(norm_diff_z>tol_vec, ...
                        this.NumberDataCenters, this.NumberPaths*this.NumberVNFs),2);
                    c = dot(norm_diff_x>tol_vec, path_reconfig_cost)...
                        + dot(sum_tz, this.node_reconfig_cost);
                otherwise
                    error('error: invalid model.');
            end
            % the difference between current solution to previous solution
            if nargout >= 2
                paths_length = sum(this.I_edge_path,1);
                n = nnz(norm_diff_z>tol_vec) + dot(paths_length,norm_diff_x>tol_vec);
            end
            if nargout >= 3
                t = (nnz(this.I_edge_path)+this.num_varz);
                r = n/t;
            end
        end
    end
    
    methods (Access = private)
        
        %% TODO 
        % get and set the network state, including flow_table, incident matrix, etc.
        function s = get_state(this)
            s.flow_table = this.FlowTable;
            s.I_node_path = this.I_node_path;
            s.I_edge_path = this.I_edge_path;
            s.I_flow_path = this.I_flow_path;
            s.path_owner = this.path_owner;
            s.As_res = this.As_res;
            s.variables = this.Variables;
        end
        function set_state(this, s)
            this.FlowTable = s.flow_table;
            this.I_node_path = s.I_node_path;
            this.I_edge_path = s.I_edge_path;
            this.I_flow_path = s.I_flow_path;
            this.path_owner = s.path_owner;
            this.As_res = s.As_res;
            this.variables = s.Variables;
        end
    end
    
    methods(Static)
        function th = THETA(t)
            persistent var_theta;
            if nargin >= 1
                var_theta = t;
            end
            th = var_theta;
        end
        %% 
        % The resource consumption cost is constant, since the slice will not
        % release/request resouce to/from the substrate network, and the resource prices
        % during the reconfiguration does not change.
        %
        % NOTE: hession matrix of the objective is only has non-zero values for x
        % (super-linear). See also <fcnHessian>.
        function [profit, grad] = fcnFastConfigProfit(vars, S, options)
            num_vars = S.num_vars;
            NP = S.NumberPaths;
            NV = S.NumberVNFs;
            if isempty(vars)
                var_path = S.Variables.x;
            else
                var_path = vars(1:NP);
                % var_node = vars((S.NumberPaths+1):S.num_vars);
                % |var_tx| and |var_tz| are auxilliary variables to transform L1 norm to
                % linear constraints and objective part.
                var_tx = vars(num_vars+(1:NP));
                var_tz = vars((num_vars+NP+1):end);
            end
            % 
            %             link_load = S.getLinkLoad;
            %             node_load = S.getNodeLoad;
            flow_rate = S.getFlowRate(var_path);
            
            
            %% objective value
            % |edge_reconfig_cost| is a vector, and we known edge-path incident matrix, so
            % we can calculate the |path_reconfig_cost|.
            if nargin>=3 && isfield(options, 'objective') && ...
                    strcmp(options.objective, 'ProfitWithReconfigurationCost')
                reconfig_cost = S.get_reconfig_cost;
                profit = -S.weight*sum(fcnUtility(flow_rate)) + reconfig_cost;
                warning('.......');
            else
                profit = -S.weight*sum(fcnUtility(flow_rate));
                %%%
                % calculate reconfigure cost by |tx| and |tz|.
                path_reconfig_cost = (S.I_edge_path)' * S.edge_reconfig_cost;
                sum_tz = sum(reshape(var_tz, S.NumberDataCenters, NP*NV),2);
                profit = profit + ...
                    dot(var_tx, path_reconfig_cost) + dot(sum_tz, S.node_reconfig_cost);
            end
            %             profit = profit - S.constant_profit; => move to SliceEx
            
            % If there is only one output argument, return the real profit (positive)
            if nargout <= 1
                profit = -profit;
            else
                %% Gradient
                % The partial derivatives are computed by dividing the variables into four
                % parts, i.e., $x,z,t_x,t_z$. Since |z| does not appear in objective
                % function, the corrsponding derivatives is zero. 

                % The partial derivatives of x
                grad = spalloc(length(vars), 1, NP+num_vars);
                for p = 1:NP
                    i = S.path_owner(p);
                    grad(p) = -S.weight/(1+S.I_flow_path(i,:)*var_path) + grad(p) ; %#ok<SPRIX>
                end
                %%%
                % The partial derivatives of t_x is the vector |path_reconfig_cost|;
                % The partial derivatives of t_z is duplicaton of |node_reconfig_cost|,
                % since the node reconfiguration cost is only depend on node.
                t_index = num_vars+(1:NP);
                grad(t_index) = path_reconfig_cost;
                grad((num_vars+1+NP):end) = repmat(S.node_reconfig_cost, NP*NV, 1);                
            end
        end
    end
end

