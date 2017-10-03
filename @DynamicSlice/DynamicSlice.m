classdef DynamicSlice < Slice & EventSender & EventReceiver
    %DynamicSlice Event-driven to dynamic configure slice.
    %   Detailed explanation goes here


    properties
        FlowArrivalRate;
        FlowServiceInterval;
        interval_events = 20;
        interval_time;
        portion_adhoc_flows_threshold = 0.1;
        b_adhoc_flows = fasle(1000,1);  % only record a window of flows.
        num_total_arrive = 0;
    end
    properties(Access={?DynamicSlice,?DynamicNetwork})
        b_ondepart = false;

        VNFCapacity;            % Capaciy of VNF instances on each node, this is configured by inter-slicing
        % Add to table VirtualLinks, this do not change between two dimensioning
        % procedure. 
        %         edge_reconfig_cost;     
        % Add to table VirtualDataCenters, this do not change between two dimensioning
        % procedure. 
        %         node_reconfig_cost;     
        vnf_reconfig_cost;
        
        Hdiag;
        % As_res;       % has been defined in <Slice>.
        %%% define reconfiguraton cost for link and node variables.
        % **flow reassginement cost*, i.e., the flow-path reconfiguration cost, which is
        % denpendent on length of paths.
        % **node reassginement cost*,
        %
        % If a resource price is p, we assume that the reconfigruation cost
        % coefficient of this resource is ¦Èp. If we adopt a varying pricing, like
        % p=av+b, then the coeffcient is ¦È(av+b).
    end
    properties(Access=private)
        z_reconfig_cost;    % for z_npf, used in optimization, the value updates during each fast reconfiguration
        x_reconfig_cost;    % for x_p, used in optimization, the value updates during each fast reconfiguration
        prev_vnf_capacity;  % last time's VNF instance capcity;
        old_variables;      % last one configuration
        changed_index;
    end
    properties(Dependent, Access=protected)
        num_varv;
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
            this.num_received_events = this.num_received_events + 1;
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
    
    methods
        function n = get.num_varv(this)
            n = this.NumberDataCenters*this.NumberVNFs;
        end
        function Hs = get.Hdiag(this)
            % For the node capacity constraint for VNF |f|, the coeffiect matrix Hs is
            % filled block-by-block. In the |k|th block (NC*NC), we put the |k|th column
            % of H_np to the diagnal of the block. Then |H_np| is duplicated into a larger
            % diagonal corresponding to all VNFs.
            %
            % The capacity of VNF instance |VNFCapacity|, has been calculate after
            % allocating resource to the slice. Some component of |VNFCapacity| might be
            % zero, since VNF f may not be located at node n.
            %
            % See also <Slice.Hrep>.
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
        
        %%%
        % * *getSliceCost*
        % Get the cost of creating a slice, overriding <Slice.getSliceCost>. The cost
        % including resource consumption cost and reconfiguration cost.
        %   cost = getSliceCost(this, pricing_policy, reconfig_cost_model)
        % |reconfig_cost_model|: reconfigure cost model, including: |'const'|,
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
        
        function b = checkFeasible(this, vars, options) %#ok<INUSD>
            if nargin <= 1
                vars = [];
            end
            if nargin <= 2
                options = struct; %#ok<NASGU>
            end
            if isempty(vars)
                b = checkFeasible@Slice(this, vars);
            else
                b = checkFeasible@Slice(this, vars(1:this.num_vars));
            end
            %             if b
            %                 % TODO add other check conditions.
            %                 switch options.Action
            %                     case 'add'
            %                     case 'remove'
            %                 end
            %             end
        end
    end

    methods (Access = {?DynamicSlice, ?DynamicNetwork})
        function setVnfCapacity(this)
            znpf = reshape(full(this.Variables.z), this.NumberDataCenters, ...
                this.NumberPaths, this.NumberVNFs);
            znpf = znpf.* full(this.I_node_path);  % compatible arithmetic operation
            this.VNFCapacity = reshape(sum(znpf,2), this.NumberDataCenters*this.NumberVNFs,1);
            if ~strcmpi(this.Parent.options.Method, 'fastconfig')
                this.prev_vnf_capacity = this.getVNFInstanceCapacity;
                if ~isequal(this.VNFCapacity(:), this.prev_vnf_capacity)
                    error('arithmetic error.')
                end
            end
            %%%
            % Equal to Hs * z_npf;
            %   Hs = this.getHs;
            %   this.VNFCapacity = Hs * this.Variables.z;
        end
        %%%
        % 
        function vc = getVNFInstanceCapacity(this, z)
            if nargin <= 1
                z = this.Variables.z;
            end
            vc = this.Hdiag * z;
        end
        
        % model = {'linear'|'const'|'none'}
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
            % old_variables is valid after adding/removing flows.
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
            %%
            % reconfiguration cost of VNF capacity
            options = getstructfields(this.Parent.options, {'Method', 'DiffNonzeroTolerance'});            
            if ~strcmpi(options.Method, 'fastconfig')
                vnf_capacity = this.VNFCapacity(:);
                diff_v = abs(vnf_capacity - this.prev_vnf_capacity);
                mid_v = 1/2*(abs(vnf_capacity) + abs(this.prev_vnf_capacity));
                nz_index_v = mid_v~=0;
                norm_diff_v = zeros(size(diff_v));
                norm_diff_v(nz_index_v) = abs(diff_v(nz_index_v)./mid_v(nz_index_v));
            end
            %%
            % reconfigration cost of flow reassignment.
            tol_vec = options.DiffNonzeroTolerance;
            switch model
                case 'linear' 
                    c = dot(diff_x, this.x_reconfig_cost) + ...
                        dot(diff_z, this.z_reconfig_cost);
                    if ~strcmpi(options.Method, 'fastconfig')
                        c = c + dot(diff_v, this.vnf_reconfig_cost);
                    end
                case 'const'
                    % logical array cannot be used as the first argument of dot.
                    c = dot(this.x_reconfig_cost, norm_diff_x>tol_vec)...
                        + dot(this.z_reconfig_cost, norm_diff_z>tol_vec);
                    if ~strcmpi(options.Method, 'fastconfig')
                        c = c + dot(this.vnf_reconfig_cost, norm_diff_v>tol_vec);
                    end
                otherwise
                    error('error: invalid model.');
            end
            % the difference between current solution to previous solution
            if nargout >= 2
                paths_length = sum(this.I_edge_path,1);
                n = nnz(norm_diff_z>tol_vec) + dot(paths_length,norm_diff_x>tol_vec);
                if ~strcmpi(options.Method, 'fastconfig')
                    % include the reconfiguration of VNF instances.
                    n = n + nnz(norm_diff_v>tol_vec);
                end
            end
            if nargout >= 3
                %   Number of edge variables, VNF assignment variables
                t = nnz(this.I_edge_path)+nnz(this.I_node_path)*this.NumberVNFs+...
                    this.num_varv;
                r = n/t;
            end
        end
    end
    
    methods (Access = private)
        function s = get_state(this)
            s.flow_table = this.FlowTable;
            s.I_node_path = this.I_node_path;
            s.I_edge_path = this.I_edge_path;
            s.I_flow_path = this.I_flow_path;
            s.path_owner = this.path_owner;
            s.variables = this.Variables;
            s.VNFCapacity = this.VNFCapacity;
        end
        function set_state(this, s)
            this.FlowTable = s.flow_table;
            this.I_node_path = s.I_node_path;
            this.I_edge_path = s.I_edge_path;
            this.I_flow_path = s.I_flow_path;
            this.path_owner = s.path_owner;
            this.variables = s.Variables;
            this.VNFCapacity = s.VNFCapacity;
        end
    end
    
    methods (Access = protected)
        profit = fastReconfigure2(this, action, options);
        profit = fastReconfigure(this, action, options);
        tf = executeMethod(this, action);
        
        %% fast slice reconfiguration when flow arriving and depaturing
        % * *flow*: flow table entries.
        %
        % *Update state*: when single flow arrive or depart, update the state record,
        % including: |I_node_path, I_edge_path, I_flow_path, path_owner| and |As|.
        % Since only add/remove entry for only one flow, this operation is more efficient
        % than <initializeState>.
        function fidx = OnAddingFlow(this, flow)  
            % temporarily allocated flow id.
            global DEBUG total_iter_num; %#ok<NUSED>
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

            % used for all methods to calculate the reconfiguration cost. and used for
            % fast reconfiguration during optimization.
            this.changed_index.x = sparse(changed_path_index);
            this.identify_change(this.changed_index.x);
            this.old_variables.x = zeros(this.NumberPaths,1);
            this.old_variables.z = zeros(this.num_varz,1);
            this.old_variables.x(~this.changed_index.x) = old_state.variables.x;
            this.old_variables.z(~this.changed_index.z) = old_state.variables.z;

            %% TODO, designed the output structure for multiple slices.
            tf = this.executeMethod('add');
            if tf == false
                fidx = [];
                this.set_state(old_state);
            end
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

            tf = this.executeMethod('remove');
            if tf == false
                this.set_state(old_state);
            end
        end
        
           
        % |parameters|: include fields x0, As, bs, Aeq, beq, lbs, ub, lb;
        % |options|: include fields fmincon_opt, costmodel;
        function [x, fval, exitflag] = optimize(this, params, options)
            switch options.CostModel
                case 'fixcost'
                    if isfield(options, 'Form') && strcmpi(options.Form, 'compact')
                        [x_compact, fval, exitflag] = ...
                            fmincon(@(x)DynamicSlice.fcnSocialWelfareCompact(x, this, ...
                            getstructfields(options, 'CostModel')), ...
                            params.x0, params.As, params.bs, params.Aeq, params.beq, ...
                            params.lb, params.ub, [], options.fmincon_opt);
                        x = zeros(this.num_vars, 1);
                        x(this.I_active_variables) = x_compact;
                    else
                        [x, fval, exitflag] = ...
                            fmincon(@(x)DynamicSlice.fcnSocialWelfare(x, this, ...
                            getstructfields(options, 'CostModel')), ...
                            params.x0, params.As, params.bs, params.Aeq, params.beq, ...
                            params.lb, params.ub, [], options.fmincon_opt);
                    end
                otherwise
                    [x, fval, exitflag] = optimize@Slice(this, params, ...
                        rmfield(options, 'CostModel'));
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
        
        %% TODO
        % optimize only on the subgraph, which carries the flows that intersection with
        % the arriving/departuring flow. this can significantly reduce the problem scale
        % when there is lots of flow, and the numbe of hops of the flow is small.
        %         function profit = fastReconfigureCompact(this, action)
        %         end
       
 
    end
    
    methods(Static)
        %%%
        % Emulation of a static property.
        function th = THETA(t)
            persistent var_theta;
            if nargin >= 1
                var_theta = t;
            end
            th = var_theta;
        end
        %%%
        % Objective value and gradient of the objective function for fast reconfiguration.
        % This is the compact form, which does not consider those in-active variables
        % (corresponds to $h_np=0$). 
        % See also <Slice.fcnHessian> and <SliceEx.fcnHessianCompact>.
        function [profit, grad] = fcnFastConfigProfitCompact(vars, S, options)
            var_x = vars(1:options.num_varx);
            num_orig_vars = options.num_varx + options.num_varz; % |num_vars| counts for |x,z,[v]|
            num_vars = length(vars);
            if isfield(options, 'num_varv')
                num_orig_vars = num_orig_vars + options.num_varv;
                var_tv = vars((num_vars-options.num_varv+1):end);
            end
            var_tx = vars(num_orig_vars+(1:options.num_varx));
            var_tz = vars(num_orig_vars+options.num_varx+(1:options.num_varz));
            flow_rate = S.getFlowRate(var_x);
            profit = -S.weight*sum(fcnUtility(flow_rate));
            profit = profit + ...
                dot(var_tx, S.x_reconfig_cost) + dot(var_tz, S.z_reconfig_cost);
            if isfield(options, 'num_varv')
                profit = profit + dot(var_tv, S.vnf_reconfig_cost);
            end
            grad = spalloc(length(vars), 1, options.num_varx+num_orig_vars);
            for p = 1:options.num_varx
                i = S.path_owner(p);
                grad(p) = -S.weight/(1+S.I_flow_path(i,:)*var_x); %#ok<SPRIX>
            end
            %%%
            % The partial derivatives of t_x is the vector |x_reconfig_cost|;
            % The partial derivatives of t_z is duplicaton of |z_reconfig_cost|,
            % since the node reconfiguration cost is only depend on node.
            grad(num_orig_vars+(1:options.num_varx)) = S.x_reconfig_cost;
            grad(num_orig_vars+options.num_varx+(1:options.num_varz)) = S.z_reconfig_cost;
            if isfield(options, 'num_varv')
                grad((num_vars-options.num_varv+1):end) = S.vnf_reconfig_cost;                
            end
        end
        %%%
        % Hessian matrix (compact form) of objective function. See also <Slice.fcnHessian>
        % and <Slice.fcnHessianCompact>.
        function hess = fcnHessianCompact(vars, ~, S, options)
            var_x = vars(1:options.num_varx);
            num_vars = length(vars);
            hess = spalloc(num_vars, num_vars, options.num_varx^2);   % non-zero elements less than | options.num_varx^2|
            for p = 1:options.num_varx
                i = S.path_owner(p);
                hess(p,1:options.num_varx) = S.weight *...
                    S.I_flow_path(i,:)/(1+(S.I_flow_path(i,:)*var_x))^2; %#ok<SPRIX>
            end
        end
        %%
        % The resource consumption cost is constant, since the slice will not
        % release/request resouce to/from the substrate network, and the resource prices
        % during the reconfiguration does not change. Therefore, the resource consumption
        % cost it not counted in the objective function. The objective function contains
        % the user utility and the reconfiguration cost.
        %
        % NOTE: hession matrix of the objective is only has non-zero values for x
        % (super-linear). See also <fcnHessian>.
        function [profit, grad] = fcnFastConfigProfit(vars, S, options)
            num_vars = S.num_vars;
            NP = S.NumberPaths;
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
            flow_rate = S.getFlowRate(var_path);
            
            %% objective value
            if nargin>=3 && isfield(options, 'Objective') && ...
                    strcmp(options.Objective, 'ProfitWithReconfigurationCost')
                reconfig_cost = S.get_reconfig_cost;
                profit = -S.weight*sum(fcnUtility(flow_rate)) + reconfig_cost;
                warning('.......');
            else
                profit = -S.weight*sum(fcnUtility(flow_rate));
                %%%
                % calculate reconfiguration cost by |tx| and |tz| (L1 approximation),
                % which is similar to calculating by |x-x0| and |z-z0|.
                profit = profit + ...
                    dot(var_tx, S.x_reconfig_cost) + dot(var_tz, S.z_reconfig_cost);
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
                    grad(p) = -S.weight/(1+S.I_flow_path(i,:)*var_path); %#ok<SPRIX>
                end
                %%%
                % The partial derivatives of t_x is the vector |x_reconfig_cost|;
                % The partial derivatives of t_z is duplicaton of |z_reconfig_cost|,
                % since the node reconfiguration cost is only depend on node.
                t_index = num_vars+(1:NP);
                grad(t_index) = S.x_reconfig_cost;
                grad((num_vars+1+NP):end) = S.z_reconfig_cost;
            end
        end
        [profit, grad] = fcnSocialWelfare(x_vars, S, options);
        %         function [profit, grad] = fcnSocialWelfareCompact(act_vars, S, options)
        %             vars = zeros(S.num_vars,1);
        %             vars(S.I_active_variables) = act_vars;
        %
        %             if nargin == 2
        %                 options = struct;
        %             end
        %             if nargout <= 1
        %                 profit = DynamicSlice.fcnSocialWelfare(vars,s,options);
        %             else
        %                 [profit, grad] = DynamicSlice.fcnSocialWelfare(vars,s,options);
        %             end
        %         end
    end
end

