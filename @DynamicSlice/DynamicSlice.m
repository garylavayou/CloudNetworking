classdef DynamicSlice < Slice & EventSender 
    %DynamicSlice Event-driven to dynamic configure slice.
    properties
        FlowArrivalRate;
        FlowServiceInterval;
        time;           % .{Current, LastDimensioning, Interval}
        portion_adhoc_flows_threshold = 0.1;
        b_adhoc_flows = zeros(1000,1);  % only record a window of flows.
        num_total_arrive = 0;
    end
    properties(Access={?DynamicSlice,?DynamicNetwork})
        b_ondepart = false;
        % the variable |b_derive_vnf| decide if update VNF instance capacity.
        %    |b_derive_vnf=true|: derive VNF instance capacity from the optimization
        %    results of |Variables.z|;
        %    |b_derive_vnf=false|: apply |Variables.v| as VNF instance capacity.
        % For first time slice dimensioning,
        %   VNF capacity is the sum of VNF instance assigment (sum of z_npf);
        % For later reallocation (ReF2), the VNF capacity should be set as the optimized
        % variables (|this.Variables.v|, which might be larger than sum of |z_npf|, due to
        % reconfiguration cost).
        b_derive_vnf = true;
        b_dim = false;
        old_net_state;
        net_changes = struct();
        vnf_reconfig_cost;
        
        % As_res;       % has been defined in <Slice>.
        %% define reconfiguraton cost for link and node variables.
        % **flow reassginement cost*, i.e., the flow-path reconfiguration cost, which is
        % denpendent on length of paths.
        % **node reassginement cost*,
        %
        % If a resource price is p, we assume that the reconfigruation cost
        % coefficient of this resource is ¦Èp. If we adopt a varying pricing, like
        % p=av+b, then the coeffcient is ¦È(av+b).
        
        a = 0.8;        % a for the history, should have a larger weight (EMA)
    end
    properties(Access=private)
        z_reconfig_cost;    % for z_npf, used in optimization, the value updates during each fast reconfiguration
        x_reconfig_cost;    % for x_p, used in optimization, the value updates during each fast reconfiguration
        old_variables;      % last one configuration, last time's VNF instance capcity is the |v| field;
        old_state;
        changed_index;
        topts;              % used in optimization, avoid passing extra arguments.
        reject_index;
        lower_bounds = struct([]);
    end
    properties(Dependent, Access=protected)
        num_varv;
        
        % Final results of VNF instance capacity on each node, this is configured by
        % inter/intra-slicing, the values are stored in |Variables.v|; This property is a
        % wrapper.
        VNFCapacity;
    end
    events
        AddFlowSucceed;
        AddFlowFailed;
        RemoveFlowSucceed;
        RemoveFlowFailed;          % NOT used.
        RequestDimensioning;    % request slice dimensioning at once
        DeferDimensioning;      % defer slice dimensioning until the network makes dicision.
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
            this.time.Current = eventData.Time;
            this.event.Count = this.event.Count + 1;
            this.event.RecentCount = this.event.RecentCount + 1;
            switch eventData.EventName
                case 'FlowArrive'
                    fidx = this.OnAddingFlow(eventData.flow);
                    % notify slice
                    if isempty(fidx)
                        notify(this, 'AddFlowFailed');
                    else
                        %%
                        % update the portion of coming adhoc flows
                        if this.getOption('Adhoc')
                            for i = 1:length(fidx)
                                if this.FlowTable{fidx(i), 'Type'} == FlowType.Adhoc
                                    record_flow_index = ...
                                        mod(this.num_total_arrive, length(this.b_adhoc_flows)) + 1;
                                    this.b_adhoc_flows(record_flow_index) = 1;
                                end
                                this.num_total_arrive = this.num_total_arrive + 1;
                            end
                        end
                        ev = eventData.event;
                        eventData = FlowEventData(ev, sl, fidx);
                        notify(this, 'AddFlowSucceed', eventData);
                    end
                case 'FlowDepart'
                    % notify slice
                    flow_id = eventData.flow;
                    this.net_changes = struct();
                    eventData = FlowEventData(eventData.event, sl, flow_id);
                    this.OnRemovingFlow(flow_id);
                    notify(this, 'RemoveFlowSucceed', eventData);
                    %%
                    % *After removing flow*: To avoid frequent resource release, a
                    % resource utilization threshold should be set. Only when the resource
                    % utilization ration is lower than the threshold, the resource will be
                    % released. 
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
            this.time = struct('Current', 0, 'LastDimensioning', 0, ...
                'DimensionInterval', 0, 'ConfigureInterval', 0);
            % Now 'ConfigureInterval' is constant, yet it also can be updated via EMA.
            this.time.ConfigureInterval = 1/(2*this.FlowArrivalRate);   
            this.event.RecentCount = 0;
            % Interval for performing dimensioning should be configurable.
            this.options = structmerge(this.options, ...
                getstructfields(slice_data, ...
                {'TimeInterval', 'EventInterval', 'Trigger', 'Adhoc'}, 'ignore'));
            if isfield(this.options, 'EventInterval')
                this.time.DimensionInterval = this.options.EventInterval/(2*this.FlowArrivalRate); % both arrival and departure events
            elseif isfield(this.options, 'TimeInterval')
                this.time.DimensionInterval = this.options.TimeInterval;
            elseif isfield(this.options, 'Trigger')
                this.time.DimensionInterval = 10/this.FlowArrivalRate;
            end
            if ~isfield(slice_data, 'ReconfigScaler')
                this.options.ReconfigScaler = 1;
            else
                this.options.ReconfigScaler = slice_data.ReconfigScaler;
            end
        end
        
        function finalize(this, node_price, link_price)
            finalize@Slice(this, node_price, link_price);
            if ~isempty(this.lower_bounds)
                this.VirtualLinks.Capacity = this.temp_vars.c;
            end
            if ~this.b_derive_vnf ...
                    && isfield(this.temp_vars,'v') && ~isempty(this.temp_vars.v)
                if this.NumberFlows > 0
                    this.Variables.v = this.temp_vars.v;
                    this.postProcessing(struct('VNFCapacity', true, 'OnlyVNFCapacity', true))
                    % Override the capacity setting in the super class.
                    this.VirtualDataCenters.Capacity = sum(reshape(this.VNFCapacity, ...
                        this.NumberDataCenters,this.NumberVNFs),2);
                end
            end
            %                 if strcmp(this.options.PricingPolicy, 'quadratic-price')
            %% Reconfiguration Cost
            % intra-slice reconfiguration: the flow reassignment cost and the VNF
            % instance reconfiguration cost is denpendtent on the resource
            % consummption in the slice;
            % inter-slice reoncfiguration: the VNF node resource allocation cost is
            % dependent on the total load of the substrat network.
            if this.NumberFlows > 0
                [~, this.VirtualLinks.ReconfigCost] = this.fcnLinkPricing(...
                    this.VirtualLinks.Price, this.VirtualLinks.Capacity);
                this.VirtualLinks.ReconfigCost = ...
                    (DynamicSlice.ETA/this.time.ConfigureInterval) * this.VirtualLinks.ReconfigCost;
                % here the |node_price| is the price of all data centers.
                [~, this.VirtualDataCenters.ReconfigCost] = this.fcnNodePricing(...
                    this.VirtualDataCenters.Price, this.VirtualDataCenters.Capacity);
                this.VirtualDataCenters.ReconfigCost = ...
                    (DynamicSlice.ETA/this.time.ConfigureInterval) * this.VirtualDataCenters.ReconfigCost;
            end
            %                 else
            %                 end
            if this.b_derive_vnf
                % We initialize the reconfiguration cost coefficient when the slice is
                % added. The coefficients is required to compute the cost after the
                % optimizaton at <DynamicCloudNetwork.onAddingSlice>.
                % After redimensioning, the former cost coefficients are still needed to
                % calculate the reconfiguration cost. So in that case, we update it before
                % performing reconfigure/redimensioning.
                this.b_derive_vnf = false;
                this.x_reconfig_cost = (this.I_edge_path)' * this.VirtualLinks.ReconfigCost;
                this.z_reconfig_cost = repmat(this.VirtualDataCenters.ReconfigCost, ...
                    this.NumberPaths*this.NumberVNFs, 1);
                this.vnf_reconfig_cost = this.Parent.options.VNFReconfigCoefficient * ...
                    repmat(this.VirtualDataCenters.ReconfigCost, this.NumberVNFs, 1);
                this.VNFCapacity = this.getVNFCapacity;
                this.old_variables = this.Variables;
                this.get_state;
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
        %             switch slice.options.FlowPattern
        %                 case {FlowPattern.RandomSingleFlow, FlowPattern.RandomMultiFlow}
        %                     ft = this.generateFlowTable(graph,
        %                 case FlowPattern.RandomInterDataCenter
        %                 case FlowPattern.RandomInterBaseStation
        %                 case FlowPattern.RandomDataCenter2BaseStation
        %                 otherwise
        %                     error('error: cannot handle the flow pattern <%s>.', ...
        %                         slice.options.FlowPattern.char);
        %             end
        %         end
        
    end
    
    methods
        function c = getLinkCapacity(this, path_vars)
            if nargin == 1
                c = this.VirtualLinks.Capacity;
            elseif ~isempty(this.lower_bounds)
                c = this.temp_vars.c;
            else
                c = getLinkCapacity@Slice(this, path_vars);
            end
        end
        function c = getNodeCapacity(this, node_vars)
            if nargin == 1
                c = this.VirtualDataCenters.Capacity;
            elseif ~isempty(this.lower_bounds)
                c = sum(reshape(this.temp_vars.v, ...
                    this.NumberDataCenters,this.NumberVNFs),2);
            else
                c = getNodeCapacity@Slice(this, node_vars);
            end
        end
        
        function c = get.VNFCapacity(this)
            c = this.Variables.v;
        end
        function set.VNFCapacity(this, value)
            if ~isfield(this.Variables, 'v') || isempty(this.Variables.v)
                this.Variables.v = value(:);    % transform to column vector
            else
                this.Variables.v(:) = value;    % automate reshape
            end
        end
        
        function n = get.num_varv(this)
            n = this.NumberDataCenters*this.NumberVNFs;
        end
    end
    
    methods
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
                reconfig_cost_model = 'none';
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
        
        %%%
        % At present this method inherit the superclass. See also <Slice.checkFeasible>.
        %         function b = checkFeasible(this, vars, options)
        %             if nargin <= 1
        %                 vars = [];
        %             else
        %                 vars = vars(1:this.num_vars);
        %             end
        %             if nargin <= 2
        %                 options = struct;
        %             end
        %             b = checkFeasible@Slice(this, vars, options);
        %             %             if b
        %             %                 % TODO add other check conditions.
        %             %                 switch options.Action
        %             %                     case 'add'
        %             %                     case 'remove'
        %             %                 end
        %             %             end
        %         end
        
        function b = isAdhocFlow(this)
            num_total = min(this.num_total_arrive,length(this.b_adhoc_flows));
            if num_total == 0 || nnz(this.b_adhoc_flows)/num_total < this.portion_adhoc_flows_threshold
                b = true;
            else
                b = false;
            end
        end
    end
    
    methods (Access = {?DynamicSlice, ?DynamicNetwork})
        % model = {'linear'|'const'|'none'}
        % Return: reconfiguration cost, number of reconfiguration, ration of
        % reconfigurations, and total number of variables.
        function c = get_reconfig_cost(this, model)
            if nargin == 1
                model = {'const'};
            elseif ischar(model)
                model = {model};
            end
            stat_names = cell(length(model),1);
            for i = 1:length(model)
                switch model{i}
                    case 'const'
                        stat_names{i} = 'Cost';
                    case 'linear'
                        stat_names{i} = 'LinearCost';
                    otherwise
                        error('[DynamicSlice]: invalid cost model <%s>.', model{i});
                end
            end
            stat = this.get_reconfig_stat(stat_names);
            c = zeros(length(stat_names),1);
            for i = 1:length(stat_names)
                c(i) = stat.(stat_names{i});
            end
        end
        
        function stat = get_reconfig_stat(this, stat_names)
            options = getstructfields(this.Parent.options, ...
                {'Method', 'DiffNonzeroTolerance', 'NonzeroTolerance'});    
            if nargin == 1
                stat_names = {'All'};
            elseif ischar(stat_names)
                stat_names = {stat_names};
            end
            % |old_variables| is still valid after adding/removing flows.
            % |old_variables| have more elements than |this.Variables.x| when removing
            % flows.
            if length(this.old_variables.x) <= length(this.Variables.x)
                new_x = this.Variables.x;
                tI_node_path = this.I_node_path;
                % adding flow, |path| will increase, and |edge| might(might not) increase;
                tI_edge_path = this.I_edge_path;
            else
                new_x = zeros(size(this.old_variables.x));
                new_x(~this.changed_index.x) = this.Variables.x;
                % removing flow: |path| will decrease, and |edge| might(might not) decrease;
                tI_node_path = this.old_state.I_node_path;
                tI_edge_path = this.old_state.I_edge_path;
            end
            diff_x = abs(new_x-this.old_variables.x);
            mid_x = 1/2*(abs(new_x)+abs(this.old_variables.x));
            nz_index_x = mid_x~=0;
            diff_x_norm = zeros(size(diff_x));
            % Here, we assume that all tiny variables (smaller than NonzeroTolerance) have been
            % eleminated. So that the following operation can identify changes.
            diff_x_norm(nz_index_x) = abs(diff_x(nz_index_x)./mid_x(nz_index_x));
            if length(this.old_variables.z) <= length(this.Variables.z)
                new_z = this.Variables.z;
            else
                new_z = zeros(size(this.old_variables.z));
                new_z(~this.changed_index.z) = this.Variables.z;
            end
            diff_z = abs(new_z-this.old_variables.z);
            mid_z = 1/2*(abs(new_z)+abs(this.old_variables.z));
            nz_index_z = mid_z~=0;
            diff_z_norm = zeros(size(diff_z));
            diff_z_norm(nz_index_z) = abs(diff_z(nz_index_z)./mid_z(nz_index_z));
            %%%
            % Reconfiguration cost of VNF capacity.
            % No reconfiguration of VNF instance for 'fastconfig'.
            if length(this.old_variables.v) <= length(this.Variables.v)
                new_vnf_capacity = this.VNFCapacity(:);
            else
                new_vnf_capacity = zeros(length(this.old_variables.v),1);
                new_vnf_capacity(~this.changed_index.v) = this.VNFCapacity(:);
            end
            diff_v = abs(new_vnf_capacity-this.old_variables.v(:));
            mid_v = 1/2*(abs(new_vnf_capacity)+abs(this.old_variables.v(:)));
            nz_index_v = mid_v~=0;
            diff_v_norm = zeros(size(diff_v));
            diff_v_norm(nz_index_v) = abs(diff_v(nz_index_v)./mid_v(nz_index_v));
            tol_vec = options.DiffNonzeroTolerance;
            stat = table;
            for i = 1:length(stat_names)
                %% Reconfiguration Cost
                if contains(stat_names{i},{'All', 'Cost'},'IgnoreCase',true)
                    % logical array cannot be used as the first argument of dot.
                    stat.Cost = dot(this.x_reconfig_cost, diff_x_norm>tol_vec)...
                        + dot(this.z_reconfig_cost, diff_z_norm>tol_vec)...
                        + dot(this.vnf_reconfig_cost, diff_v_norm>tol_vec);
                end
                if contains(stat_names{i},{'All', 'LinearCost'},'IgnoreCase',true)
                    stat.LinearCost = dot(diff_x, this.x_reconfig_cost) + ...
                        dot(diff_z, this.z_reconfig_cost) + ...
                        dot(diff_v, this.vnf_reconfig_cost);
                    stat.LinearCost = stat.LinearCost * this.options.ReconfigScaler;
                end
                %% Number of Variables
                if contains(stat_names{i},{'All', 'NumberVariables'},'IgnoreCase',true)
                    stat.NumberVariables = ...
                        nnz(tI_edge_path)+nnz(tI_node_path)*this.NumberVNFs;
                    stat.NumberVariables = stat.NumberVariables + nnz(mid_v);
                end            
                %% Number of reconfigurations between current and previous solution
                % NOTE: resource allocation and release will not take place at the same
                % time. 
                %   Number of edge variables, VNF assignment variables
                %   the reconfiguration of VNF instances.
                if contains(stat_names{i},{'All', 'NumberReconfigVariables'},'IgnoreCase',true)
                    paths_length = sum(tI_edge_path,1);
                    stat.NumberReconfigVariables = nnz(diff_z_norm>tol_vec) ...
                        + dot(paths_length,diff_x_norm>tol_vec)...
                        + nnz(diff_v_norm>tol_vec);
                end
                %% Number of Flows
                if contains(stat_names{i},{'All', 'NumberFlows'},'IgnoreCase',true)
                    % For convenience of comparison, we store the number of flows including the
                    % removed one.
                    stat.NumberFlows = max(this.NumberFlows, height(this.old_state.flow_table));
                end
                %% Time
                if contains(stat_names{i},{'All', 'Time'},'IgnoreCase',true)
                    stat.Time = this.time.Current;
                end
                %% Number of Reconfigured Flows
                if contains(stat_names{i},{'All', 'NumberReconfigFlows'},'IgnoreCase',true)
                    np = numel(diff_x_norm);
                    nv = this.NumberVNFs;
                    nn = numel(diff_z_norm)/(nv*np);
                    diff_path = (diff_x_norm>tol_vec)' + ...
                        sum(sum(reshape(diff_z_norm>tol_vec,nn,np,nv),3),1);
                    if length(this.path_owner) <= length(this.old_state.path_owner)
                        stat.NumberReconfigFlows = ...
                            numel(unique(this.old_state.path_owner(diff_path~=0)));
                    else
                        stat.NumberReconfigFlows = ...
                            numel(unique(this.path_owner(diff_path~=0)));
                    end
                end
            end
        end
    end
    
    methods (Access = private)
        function s = get_state(this)
            this.old_state.flow_table = this.FlowTable;
            this.old_state.I_node_path = this.I_node_path;
            this.old_state.I_edge_path = this.I_edge_path;
            this.old_state.I_flow_path = this.I_flow_path;
            this.old_state.path_owner = this.path_owner;
            this.old_state.variables = this.Variables;
            this.old_state.x_reconfig_cost = this.x_reconfig_cost;
            this.old_state.z_reconfig_cost = this.z_reconfig_cost;
            this.old_state.vnf_reconfig_cost = this.vnf_reconfig_cost;
            this.old_state.As_res = this.As_res;
            if nargout >= 1
                s = this.old_state;
            end
            this.old_state.link_load = this.VirtualLinks.Load;
            this.old_state.link_capacity = this.VirtualLinks.Capacity;
            this.old_state.node_load = this.VirtualDataCenters.Load;
            this.old_state.node_capacity = this.VirtualDataCenters.Capacity;
            this.old_state.vnf_load = this.getVNFCapacity;
            this.old_state.vnf_capacity = this.VNFCapacity;
        end
        function set_state(this)
            this.FlowTable = this.old_state.flow_table;
            % path is shallow copyed, its local identifier might have been changed.
            pid = 1;
            for fid = fidx
                for p = 1:this.FlowTable{fid, 'Paths'}.Width
                    this.FlowTable{fid, 'Paths'}.paths{p}.local_id = pid;
                    pid = pid + 1;
                end
            end
            this.I_node_path = this.old_state.I_node_path;
            this.I_edge_path = this.old_state.I_edge_path;
            this.I_flow_path = this.old_state.I_flow_path;
            this.path_owner = this.old_state.path_owner;
            this.Variables = this.old_state.variables;
            this.x_reconfig_cost = this.old_state.x_reconfig_cost;
            this.z_reconfig_cost = this.old_state.z_reconfig_cost;
            this.vnf_reconfig_cost = this.old_state.vnf_reconfig_cost;
            this.As_res = this.old_state.As_res;
        end
        function convert(this, x, ~)
            NP = this.NumberPaths;
            this.temp_vars.x = x(1:NP);
            this.temp_vars.z = x((NP+1):this.num_vars);
            offset = this.num_vars;
            this.temp_vars.v = x(offset+(1:this.num_varv));
            offset = offset + this.num_varv;
            this.temp_vars.tx = x(offset+(1:NP));
            this.temp_vars.tz = x(offset+((NP+1):this.num_vars));
            offset = offset + this.num_vars;
            this.temp_vars.tv = x(offset+(1:this.num_varv));
            if ~isempty(this.lower_bounds)
                offset = offset + this.num_varv;
                this.temp_vars.c = x(offset+(1:this.NumberVirtualLinks));
            end
            if nargin >= 3
                this.Variables = getstructfields(this.temp_vars, {'x','z','v','c'}, 'ignore');
            end
        end
        function [profit, cost] = handle_zero_flow(this, new_opts)
            cost = this.getSliceCost(new_opts.PricingPolicy);
            profit = -cost;
            this.Variables.x = [];
            this.Variables.z = [];
            this.VirtualLinks{:,'Load'} = 0;
            this.VirtualDataCenters{:,'Load'} = 0;
            % When the number of flows reduced to 0, we do not change the VNF instance
            % capcity, so that the reconfiguration cost is reduced. However, at the
            % next time when a flow arrive, the reconfiguration cost might be larger.
            % Another method is to set the VNF instance capacity to 0.
            %                 if nargout >= 1
            %                     this.VNFCapacity = 0;
            %                 end
        end
    end
 
    
    methods (Access = {?Slice, ?CloudNetwork, ?SliceFlowEventDispatcher})
        [utility, node_load, link_load] = priceOptimalFlowRate(this, x0, new_opts);
    end
    methods (Access = {?Slice, ?CloudNetwork})
        %%
        % |new_opts|:
        % * *Model*: |fixcost|.
        function [profit, cost] = optimalFlowRate( this, new_opts )
            if nargin <= 1
                new_opts = struct;
            end
            if this.NumberFlows == 0
                [profit, cost] = this.handle_zero_flow(new_opts);
            else
                [profit,cost] = optimalFlowRate@Slice( this, new_opts );
                if isfield(new_opts, 'CostModel') && strcmpi(new_opts.CostModel, 'fixcost')
                    profit = profit - cost;
                end
                if nargout >= 1
                    % When output argument specified, we finalize the VNF capacity.
                    % After reconfiguration VNF instance capcity has changed.
                    this.VNFCapacity = this.getVNFCapacity;
                end
            end
        end
    end
    methods (Access = protected)
        [profit,cost] = fastReconfigure2(this, action, options);
        [profit,cost] = fastReconfigure(this, action, options);
        [exitflag,fidx] = executeMethod(this, action);
                
        %% Deep Copy
        function this = copyElement(ds)
            this = copyElement@Slice(ds);
            %%
            % The copyed version may not have the same targets as the copy source. We can
            % mannually update the target/listener list using AddListener/RemoveListener.
            %{
              temp = copyElement@EventSender(ds);
              this.targets = temp.targets;
              this.listeners = temp.listeners;
            %}
            this.listeners = ListArray('event.listener');
            this.targets = ListArray('EventReceiver');
        end
        
        %% fast slice reconfiguration when flow arriving and depaturing
        % * *flow*: flow table entries.
        %
        % *Update state*: when single flow arrive or depart, update the state record,
        % including: |I_node_path, I_edge_path, I_flow_path, path_owner| and |As|.
        % Since only add/remove entry for only one flow, this operation is more efficient
        % than <initializeState>.
        function fidx = OnAddingFlow(this, flows)
            % temporarily allocated flow id.
            global DEBUG total_iter_num; %#ok<NUSED>
            num_exist_flows = height(this.FlowTable);
            num_new_flows = height(flows);
            num_exist_paths = this.NumberPaths;
            fidx = (1:num_new_flows) + num_exist_flows;
            if num_exist_flows == 0
                flows{:,'Identifier'} = 1:num_new_flows;
            else
                flows{:,'Identifier'} = (1:num_new_flows) + this.FlowTable{end,'Identifier'}; % temporary identifier
            end
            if num_exist_flows == 0
                pid = 0;
            else
                pid = this.FlowTable{end, 'Paths'}.paths{end}.local_id;
            end
            
            this.get_state; % if fail to adding flow, recover flow table
            
            this.FlowTable = [this.FlowTable; flows];
            changed_path_index((num_exist_paths+1):this.NumberPaths) = true;
            for fid = fidx
                for p = 1:this.FlowTable{fid, 'Paths'}.Width
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
            
            %% Previous State
            % In _onAddingFlow_ and _onRemovingFlow_, the slice topology does not change.
            % The data centers and NVF type do not changes, so there is no new/deleted VNF
            % instance variables.
            %
            % Record the state changes after flow is added to the flow table to maintain
            % the full set of variables, different from that in _<onRemovingFlow>_.
            this.identify_change(changed_path_index);
            this.old_variables.x = zeros(this.NumberPaths,1);
            this.old_variables.z = zeros(this.num_varz,1);
            this.old_variables.x(~this.changed_index.x) = this.old_state.variables.x;
            this.old_variables.z(~this.changed_index.z) = this.old_state.variables.z;
            this.x_reconfig_cost = (this.I_edge_path)' * this.VirtualLinks.ReconfigCost;
            this.z_reconfig_cost = repmat(this.VirtualDataCenters.ReconfigCost, ...
                this.NumberPaths*this.NumberVNFs, 1);
            this.topts.x_reconfig_cost = this.options.ReconfigScaler*this.x_reconfig_cost;
            this.topts.z_reconfig_cost = this.options.ReconfigScaler*this.z_reconfig_cost;
            this.topts.old_variables_x = this.old_variables.x;
            this.topts.old_variables_z = this.old_variables.z;
            [ef, ~] = this.executeMethod('add');
            %% Failure handling
            % Zero-rate: now we leave it unhandled, the zero-rate flow will stay in
            % the slice. It may be allocated with resource at later stage. On the
            % other hand, we can remove the zero-rate flows from the slice (reject).
            if ef < 0
                fidx = [];
                this.set_state();
            end            
        end
        
        function ef = OnRemovingFlow(this, fid)
            global DEBUG; %#ok<NUSED>
            fidx = this.FlowTable.Identifier == fid;
            if isempty(find(fidx==true,1))
                warning('flow (identifier %d) not in slice (identifer %d).',...
                    fid, this.Identifier);
            end
            this.get_state; % if fail to adding flow, recover slice data
            changed_path_index = false(this.NumberPaths,1);
            flow_index = find(fidx);
            for fid = flow_index
                for p = 1:this.FlowTable{fid, 'Paths'}.Width
                    pid = this.FlowTable{fid, 'Paths'}.paths{p}.local_id;
                    changed_path_index(pid) = true;
                end
            end
            %% Previous State
            % |old_variables|(to be renamed) is used for all methods to calculate the
            % reconfiguration cost. 
            %
            % When removing flows, the deleted variables' value is treated as zeros, its
            % difference is constant (0-x0). Reconfiguration cost corresponding to the
            % previous state, should be updated before performing optimization.
            % *Link reconfigure cost* is a vector, and we know edge-path incident matrix,
            % so we can calculate the |x_reconfig_cost| for all paths.
            % 
            % Record changes before flows are removed from flow table, _identify_change_
            % is called. 
            this.identify_change(changed_path_index); 
            this.old_variables.x = this.old_state.variables.x;
            this.old_variables.z = this.old_state.variables.z;
            this.x_reconfig_cost = (this.I_edge_path)' * this.VirtualLinks.ReconfigCost;
            this.z_reconfig_cost = repmat(this.VirtualDataCenters.ReconfigCost, ...
                this.NumberPaths*this.NumberVNFs, 1);
            this.topts.x_reconfig_cost = ...
                this.options.ReconfigScaler*this.x_reconfig_cost(~this.changed_index.x);
            this.topts.z_reconfig_cost = ...
                this.options.ReconfigScaler*this.z_reconfig_cost(~this.changed_index.z);
            this.topts.old_variables_x = this.old_variables.x(~this.changed_index.x);
            this.topts.old_variables_z = this.old_variables.z(~this.changed_index.z);
            % Before performing optimization, there is not situation that VNF instances
            % will be removed. So there is no need to copy the |old_variable.v| like x/z.
            
            %% Update flow table and local path id.
            % After removing flows, re-allocate local identifier to subsequent flows (flow
            % ID) and paths (path ID). 
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
            %% Update incident matrix
            % remove the path/flow-related rows/columns.
            this.I_edge_path(:, changed_path_index) = [];
            this.I_node_path(:, changed_path_index) = [];
            this.I_flow_path(:, changed_path_index) = [];
            this.I_flow_path(fidx, :) = [];
            
            [ef,~] = this.executeMethod('remove');
            if ef < 0
                this.set_state();
            end
        end
        
        
        % |parameters|: include fields x0, As, bs, Aeq, beq, lbs, ub, lb;
        % |options|: include fields fmincon_opt, CostModel, PricingPolicy (if
        %       CostModel='fixcost'); 
        function [x, fval] = optimize(this, params, options)
            if isfield(options, 'CostModel') && strcmpi(options.CostModel, 'fixcost')
                if isfield(options, 'Form') && strcmpi(options.Form, 'compact')
                    [x_compact, fval, exitflag, output] = ...
                        fmincon(@(x)DynamicSlice.fcnSocialWelfareCompact(x, this, ...
                        getstructfields(options, {'CostModel', 'num_orig_vars'})), ...
                        params.x0, params.As, params.bs, params.Aeq, params.beq, ...
                        params.lb, params.ub, [], options.fmincon_opt);
                    x = zeros(this.num_vars, 1);
                    x(this.I_active_variables) = x_compact;
                else
                    [x, fval, exitflag, output] = ...
                        fmincon(@(x)DynamicSlice.fcnSocialWelfare(x, this, ...
                        getstructfields(options, 'CostModel')), ...
                        params.x0, params.As, params.bs, params.Aeq, params.beq, ...
                        params.lb, params.ub, [], options.fmincon_opt);
                end
                this.interpretExitflag(exitflag, output.message);
                assert(this.checkFeasible(x, ...
                    struct('ConstraintTolerance', options.fmincon_opt.ConstraintTolerance)), ...
                    'error: infeasible solution.');
            else
                [x, fval] = optimize@Slice(this, params, rmstructfields(options, 'CostModel'));     
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
            this.changed_index.x = logical(sparse(changed_path_index));
            base_z_index = false(this.NumberDataCenters, this.NumberPaths);
            base_z_index(:,changed_path_index) = true;
            if ~isempty(fieldnames(this.net_changes))
                base_z_index(this.net_changes.DCIndex,:) = true;
            end
            base_z_index = sparse(base_z_index);
            this.changed_index.z = repmat(base_z_index(:), this.NumberVNFs, 1);
            if ~isempty(fieldnames(this.net_changes))
                this.changed_index.v = sparse(this.NumberDataCenters, this.NumberVNFs);
                this.changed_index.v(this.net_changes.DCIndex,:) = true;
                this.changed_index.v = logical(this.changed_index.v(:));
            end
        end
        
        %% TODO
        % optimize only on the subgraph, which carries the flows that intersection with
        % the arriving/departuring flow. this can significantly reduce the problem scale
        % when there is lots of flow, and the numbe of hops of the flow is small.
        %         function profit = fastReconfigureCompact(this, action)
        %         end
        function temp_vars = get_temp_variables(this, bfull)
            if nargin >= 2 && bfull
                temp_vars = [this.temp_vars.x; this.temp_vars.z; this.temp_vars.v; ...
                    this.temp_vars.tx; this.temp_vars.tz; this.temp_vars.tv];
            else
                temp_vars = [this.temp_vars.x; this.temp_vars.z; this.temp_vars.v];
            end
        end
        
        function release_resource_description(this)
            %% Release Unused Resource Description
            % release unused resource description: datacenters, nodes, links;
            %   a datacenter is unused when its capcacity is zero;
            %   a node is unused only when all adjecent links are unused;
            %   a link is unused when its capcity is zero and no paths use it;
            %
            this.get_state;
            b_removed_dcs = this.VirtualDataCenters.Capacity <= eps;
            b_removed_links = this.VirtualLinks.Capacity <= eps;
            for i = 1:this.NumberFlows
                pathlist = this.FlowTable{i, 'Paths'};
                for j = 1:pathlist.Width
                    path = pathlist.paths{j};
                    h = path.node_list(1:(end-1));
                    t = path.node_list(2:end);
                    eid = this.Topology.IndexEdge(h,t);
                    b_removed_links(eid) = false;
                end
            end
            b_removed_nodes = this.Topology.Remove([], b_removed_links);
            
            %% Update variables
            % see also <onRemoveFlow>.
            this.net_changes.NodeIndex = b_removed_nodes;
            this.net_changes.DCIndex = b_removed_dcs;
            this.net_changes.LinkIndex = b_removed_links;
            this.identify_change(false(this.NumberPaths,1));
            this.I_node_path(b_removed_dcs, :) = [];
            this.I_edge_path(b_removed_links, :) = [];
            this.Variables.z(this.changed_index.z) = [];
            this.Variables.v(this.changed_index.v) = [];
            %             this.vnf_reconfig_cost(this.changed_index.v) = [];
            
            %% Update resources
            this.PhysicalNodeMap{:,'VirtualNode'} = 0;
            this.VirtualNodes(b_removed_nodes,:) = [];
            this.PhysicalNodeMap{this.VirtualNodes.PhysicalNode,'VirtualNode'} = ...
                (1:this.NumberVirtualNodes)';
            this.PhysicalLinkMap{:,'VirtualLink'} = 0;
            this.VirtualLinks(b_removed_links,:) = [];
            this.PhysicalLinkMap{this.VirtualLinks.PhysicalLink,'VirtualLink'} = ...
                (1:this.NumberVirtualLinks)';
            dc_node_index = this.VirtualDataCenters.VirtualNode;
            this.VirtualDataCenters(b_removed_dcs, :) = [];
            this.VirtualNodes.DataCenter(dc_node_index(b_removed_dcs)) = 0;
            dc_node_index(b_removed_dcs) = [];
            this.VirtualNodes.DataCenter(dc_node_index) = 1:this.NumberDataCenters;
        end
        
        function [b, vars] = postProcessing(this, options)
            if nargin >= 2 && isfield(options, 'VNFCapacity')
                tol_zero = this.Parent.options.ConstraintTolerance;
                this.Variables.v(this.Variables.v<tol_zero*max(this.Variables.v)) = 0;
                % [TODO] May need post processing for VNF capacity constraint;
                if isfield(options, 'OnlyVNFCapacity') && options.OnlyVNFCapacity
                    return;
                end
            end
            if nargout <= 1
                b = postProcessing@Slice(this);
            else
                [b, vars] = postProcessing@Slice(this);
            end
        end        
    end
    
    methods(Static)
        %%%
        % Emulation of a static property.
        function th = ETA(t)
            %% Reconfiguration Cost Coefficient
            persistent var_theta;
            if nargin >= 1
                var_theta = t;
            end
            th = var_theta;
        end
    end
    methods(Static, Access = protected)    
        [profit, grad]= fcnProfit(vars, slice, options);
        hess = fcnHessian(var_x, ~, slice, options);
        % Inherit <Slice.fcnProfit>, dynamically calling <fcnProfit> of <Slice> or the
        % subclasses. 
        % Note: options must provide the 'num_orig_vars' field.        
        function [profit, grad] = fcnProfitCompact(act_vars, slice, options)
            vars = zeros(options.num_orig_vars,1);
            vars(slice.I_active_variables) = act_vars;
            
            % we extend the active variables by adding zeros to the inactive ones.
            if nargout <= 1
                profit = DynamicSlice.fcnProfit(vars, slice, options);
            else
                [profit, grad] = DynamicSlice.fcnProfit(vars, slice, options);
                % eliminate the inactive variable's derivatives.
                grad = grad(slice.I_active_variables);
            end
        end
        
        %%%
        % Objective value and gradient of the objective function for fast reconfiguration.
        % This is the compact form, which does not consider those in-active variables
        % (corresponds to $h_np=0$).
        % Used by <DynamicSlice.fastReconfigure> and <DynamicSlice.fastReconfigure2>.
        % See also <Slice.fcnHessian> and <SliceEx.fcnHessianCompact>.
        function [profit, grad] = fcnFastConfigProfitCompact(vars, slice, options)
            var_x = vars(1:options.num_varx);
            num_basic_vars = options.num_varx + options.num_varz; % |num_vars| counts for |x,z,[v]|
            num_vars = length(vars);
            if isfield(options, 'num_varv')
                num_basic_vars = num_basic_vars + options.num_varv;
                var_tv = vars((num_vars-options.num_varv+1):end);
            end
            var_tx = vars(num_basic_vars+(1:options.num_varx));
            var_tz = vars(num_basic_vars+options.num_varx+(1:options.num_varz));
            flow_rate = slice.getFlowRate(var_x);
            profit = -slice.weight*sum(fcnUtility(flow_rate));
            profit = profit + dot(var_tx, slice.topts.x_reconfig_cost) + ...
                dot(var_tz, slice.topts.z_reconfig_cost);
            if isfield(options, 'num_varv')
                profit = profit + dot(var_tv, slice.topts.vnf_reconfig_cost);
            end
            
            grad = spalloc(num_vars, 1, options.num_varx+num_vars/2);
            for p = 1:options.num_varx
                i = slice.path_owner(p);
                grad(p) = -slice.weight/(1+slice.I_flow_path(i,:)*var_x); %#ok<SPRIX>
            end
            %%%
            % The partial derivatives of t_x is the vector |x_reconfig_cost|;
            % The partial derivatives of t_z is duplicaton of |z_reconfig_cost|,
            % since the node reconfiguration cost is only depend on node.
            grad(num_basic_vars+(1:options.num_varx)) = slice.topts.x_reconfig_cost;
            grad(num_basic_vars+options.num_varx+(1:options.num_varz)) = ...
                slice.topts.z_reconfig_cost;
            if isfield(options, 'num_varv')
                grad((num_vars-options.num_varv+1):end) = slice.topts.vnf_reconfig_cost;
            end
        end
        %%%
        % Hessian matrix (compact form) of objective function. See also <Slice.fcnHessian>
        % and <Slice.fcnHessianCompact>.
        % Called by <priceOptimalFlowRate>, specify method.
        function hess = fcnHessianCompact(act_vars, lambda, S, options)
            if isfield(options, 'Method') && contains(options.Method, 'dimconfig')
                vars = zeros(options.num_orig_vars,1);
                vars(S.I_active_variables) = act_vars;
                hess = DynamicSlice.fcnHessian(vars, lambda, S, options);
                hess = hess(S.I_active_variables, S.I_active_variables);
            else
                var_x = act_vars(1:options.num_varx);
                num_vars = length(act_vars);
                hess = spalloc(num_vars, num_vars, options.num_varx^2);   % non-zero elements less than | options.num_varx^2|
                for p = 1:options.num_varx
                    i = S.path_owner(p);
                    hess(p,1:options.num_varx) = S.weight *...
                        S.I_flow_path(i,:)/(1+(S.I_flow_path(i,:)*var_x))^2; %#ok<SPRIX>
                end
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
        function [profit, grad] = fcnFastConfigProfit(vars, slice, options)
            NP = slice.NumberPaths;
            var_path = vars(1:NP);
            num_basic_vars = slice.num_vars;
            num_vars = length(vars);
            % var_node = vars((S.NumberPaths+1):S.num_vars);
            % |var_tx| and |var_tz| are auxilliary variables to transform L1 norm to
            % linear constraints and objective part.
            if nargin >= 3  % for _fastReconfigure2_, including |v|, and |tv|
                num_basic_vars = num_basic_vars + options.num_varv;
                var_tv = vars((num_vars-options.num_varv+1):end);
            end
            var_tx = vars(num_basic_vars+(1:NP));
            var_tz = vars(num_basic_vars+NP+(1:slice.num_varz));
            flow_rate = slice.getFlowRate(var_path);
            %% objective value
            profit = -slice.weight*sum(fcnUtility(flow_rate));
            %%%
            % calculate reconfiguration cost by |tx| and |tz| (L1 approximation),
            % which is similar to calculating by |x-x0| and |z-z0|.
            profit = profit + dot(var_tx, slice.topts.x_reconfig_cost) + ...
                dot(var_tz, slice.topts.z_reconfig_cost);
            if nargin >= 3 % for _fastConfigure2_
                profit = profit + dot(var_tv, slice.topts.vnf_reconfig_cost);
            end
            % If there is only one output argument, return the real profit (positive)
            if nargout <= 1
                profit = -profit;
            else
                %% Gradient
                % The partial derivatives are computed by dividing the variables into four
                % parts, i.e., $x,z,t_x,t_z$. Since |z| does not appear in objective
                % function, the corrsponding derivatives is zero.                
                % The partial derivatives of x
                grad = spalloc(num_vars, 1, NP+num_vars/2);
                for p = 1:NP
                    i = slice.path_owner(p);
                    grad(p) = -slice.weight/(1+slice.I_flow_path(i,:)*var_path); %#ok<SPRIX>
                end
                %%%
                % The partial derivatives of t_x is the vector |x_reconfig_cost|;
                % The partial derivatives of t_z is duplicaton of |z_reconfig_cost|,
                % since the node reconfiguration cost is only depend on node.
                grad(num_basic_vars+(1:NP)) = slice.topts.x_reconfig_cost;
                grad(num_basic_vars+NP+(1:slice.num_varz)) = slice.topts.z_reconfig_cost;
                if nargin >= 3
                    grad((num_vars-slice.num_varv+1):end) = slice.topts.vnf_reconfig_cost;
                end
            end
        end
        
        
        [profit, grad] = fcnSocialWelfare(x_vars, S, options);
        function [profit, grad] = fcnSocialWelfareCompact(act_vars, slice, options)
            vars = zeros(options.num_orig_vars,1);
            vars(slice.I_active_variables) = act_vars;
            
            if nargin == 2
                options = struct;
            end
            %% 
            % Here, we do not wish the subclasses dynamically override this static method,
            % so we MUST use class name to access the method.
            if nargout <= 1
                profit = DynamicSlice.fcnSocialWelfare(vars,slice,options);
            else
                [profit, grad] = DynamicSlice.fcnSocialWelfare(vars,slice,options);
                grad = grad(slice.I_active_variables);
            end
        end
    end
end
