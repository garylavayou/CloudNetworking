classdef DynamicSlice < Slice & EventSender
    %DynamicSlice Event-driven to dynamic configure slice.
    properties (Constant)
        GLOBAL_OPTIONS = StaticProperties;
    end
    
    properties(Constant, Access = private)
        NUM_MEAN_BETA = 15;
        ENABLE_DYNAMIC_NORMALIZER = true;
        GET_BETA_METHOD = 'Average';    % 'Average', 'ExponetialMovingAverage';
    end
    
    properties
        FlowArrivalRate;
        FlowServiceInterval;
        time;           % .{Current, LastDimensioning, Interval}
        portion_adhoc_flows_threshold = 0.1;
        b_adhoc_flows = zeros(1000,1);  % only record a window of flows.
        num_total_arrive = 0;
        raw_beta;
        raw_cost;
        raw_costv;
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
        b_dim = true;      % reset each time before reconfigurtion.
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
        % coefficient of this resource is ��p. If we adopt a varying pricing, like
        % p=av+b, then the coeffcient is ��(av+b).
        
        a = 0.8;        % a for the history, should have a larger weight (EMA)
        
        %% options for slice dimensioning schedule algorithm
        % 'omega_upper':
        % 'omega_lower':
        % 'alpha':
        % 'series_length'
        % 'trend_length'
        sh_options = struct('omega_upper', 0.97, 'omega_lower', 0.85, ...
            'alpha', 0.3, 'series_length', 40, 'trend_length', 20);
        %% data for slice dimensioning schdule algorithm
        % 'omegas':
        % 'profits'
        % 'omega_trend':
        % 'profit_trend':
        sh_data = struct('omegas', [] , 'profits', [], ...
            'omega_trend', struct('ascend', [], 'descend', []),...
            'profit_trend', struct('ascend', [], 'descend', []),...
            'index', 0, 'reserve_dev', 0);
        invoke_method = 0;
    end
    properties(Access=private)
        z_reconfig_cost;    % for z_npf, used in optimization, the value updates during each fast reconfiguration
        x_reconfig_cost;    % for x_p, used in optimization, the value updates during each fast reconfiguration
        old_variables;      % last one configuration, last time's VNF instance capcity is the |v| field;
        old_state;
        changed_index;
        diff_state;         % reset each time before reconfiguration.
        topts;              % used in optimization, avoid passing extra arguments.
        reject_index;
        lower_bounds = struct([]);
        upper_bounds = struct([]);
        max_flow_rate;
    end
    properties(Dependent)
        UnitReconfigureCost;
    end
    properties(Dependent, SetAccess = protected)
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
                'DimensionInterval', 0, 'ConfigureInterval', 0, 'DimensionIntervalModified', 0);
            % Now 'ConfigureInterval' is constant, yet it also can be updated via EMA.
            this.time.ConfigureInterval = 1/(2*this.FlowArrivalRate);
            this.event.RecentCount = 0;
            % Interval for performing dimensioning should be configurable.
            this.options = structmerge(this.options, ...
                getstructfields(slice_data, ...
                {'TimeInterval', 'EventInterval', 'Trigger', 'Adhoc'}, ...
                'ignore'));
            this.options = structmerge(this.options, ...
                getstructfields(slice_data, 'ReconfigMethod' , 'error'));
            switch this.options.ReconfigMethod 
                case {ReconfigMethod.DimconfigReserve, ReconfigMethod.DimconfigReserve0,...
                        ReconfigMethod.FastconfigReserve}
                    this.options = structmerge(this.options, ...
                        getstructfields(slice_data, 'bReserve', 'default', 2));
                case ReconfigMethod.DimBaseline
                    this.options = structmerge(this.options, ...
                        getstructfields(slice_data, 'bReserve', 'default', 0));
                    this.sh_options.omega_upper = 1;
                    this.sh_options.omega_lower = 0.9;
                case ReconfigMethod.Dimconfig
                    this.options = structmerge(this.options, ...
                        getstructfields(slice_data, 'bReserve', 'default', 0));
                    this.sh_options.omega_upper = 1;
                otherwise
                    this.options = structmerge(this.options, ...
                        getstructfields(slice_data, 'bReserve', 'default', 1));
            end
            switch this.options.ReconfigMethod
                case {ReconfigMethod.DimBaseline, ReconfigMethod.Dimconfig}
                    u = 0.2;
                case ReconfigMethod.DimconfigReserve0
                    u = 0.1;
                otherwise
                    u = 0;
            end
            if this.options.ReconfigMethod == ReconfigMethod.DimconfigReserve
                this.sh_options = structmerge(this.sh_options, ...
                    getstructfields(slice_data, 'UtilizationVariance', 'default', 0.05));
            else
                this.sh_options.UtilizationVariance = u;
            end
            if this.options.bReserve
                if this.options.ReconfigMethod == ReconfigMethod.DimconfigReserve0
                    this.options.Reserve = 0;
                else
                    this.options = structmerge(this.options, ...
                        getstructfields(slice_data, 'Reserve','default', 0.9));
                end
            end
            if isfield(this.options, 'EventInterval')
                this.time.DimensionInterval = this.options.EventInterval/(2*this.FlowArrivalRate); % both arrival and departure events
                this.time.MinDimensionInterval = 1/5*this.time.DimensionInterval;
            elseif isfield(this.options, 'TimeInterval')
                this.time.DimensionInterval = this.options.TimeInterval;
            elseif isfield(this.options, 'Trigger')
                this.time.DimensionInterval = 10/this.FlowArrivalRate;
            end
            this.time.DimensionIntervalModified = this.time.DimensionInterval;
            if ~DynamicSlice.ENABLE_DYNAMIC_NORMALIZER
                if ~isfield(slice_data, 'ReconfigScaler')
                    this.options.ReconfigScaler = 1;
                else
                    this.options.ReconfigScaler = slice_data.ReconfigScaler;
                end
            end
            if isfield(slice_data, 'penalty')
                this.options.penalty = slice_data.penalty;
            end
        end
        
        function finalize(this, node_price, link_price)
            global DEBUG g_results;
            finalize@Slice(this, node_price, link_price);
            if isfield(this.temp_vars,'v') && ~isempty(this.temp_vars.v)
                if this.NumberFlows > 0
                    % Override the capacity setting in the super class.
                    this.VirtualDataCenters.Capacity = full(sum(reshape(this.VNFCapacity, ...
                        this.NumberDataCenters,this.NumberVNFs),2));
                else
                  warning('zero flow, Variables.v not initialized.');  
                end
            elseif this.NumberFlows > 0
                this.Variables.v = this.getVNFCapacity;
            else
                warning('zero flow, Variables.v not initialized.');
            end
            if isfield(this.temp_vars,'c') && ~isempty(this.temp_vars.c)
                this.VirtualLinks.Capacity = full(this.temp_vars.c);
            end

            %                 if strcmp(this.options.PricingPolicy, 'quadratic-price')
            %% Reconfiguration Cost
            % intra-slice reconfiguration: the flow reassignment cost and the VNF
            % instance reconfiguration cost is denpendent on the resource
            % consummption in the slice;
            % inter-slice reoncfiguration: the VNF node resource allocation cost is
            % dependent on the total load of the substrat network.
            if this.NumberFlows > 0
                [~, this.VirtualLinks.ReconfigCost] = this.fcnLinkPricing(...
                    this.VirtualLinks.Price, this.VirtualLinks.Capacity);
                eta = DynamicSlice.GLOBAL_OPTIONS.get('eta');
                this.VirtualLinks.ReconfigCost = ...
                    (eta/this.time.ConfigureInterval) * this.VirtualLinks.ReconfigCost;
                % here the |node_price| is the price of all data centers.
                [~, this.VirtualDataCenters.ReconfigCost] = this.fcnNodePricing(...
                    this.VirtualDataCenters.Price, this.VirtualDataCenters.Capacity);
                this.VirtualDataCenters.ReconfigCost = ...
                    (eta/this.time.ConfigureInterval) * this.VirtualDataCenters.ReconfigCost;
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
                this.old_variables = this.Variables;
                this.get_state;
                this.max_flow_rate = this.FlowTable.Rate;
                if this.ENABLE_DYNAMIC_NORMALIZER
                    %% Initialize L1-approximation normalizer
                    % change of magnitude:
                    %   x -> 0; 0 -> x:    |x|      : small portion
                    %   x -> x1:           |x-x1|   : large portion
                    % we choose beta to be (1/2) of the magnitude.
                    field_names = {'x', 'z', 'v'};
                    cost_entries = {'x_reconfig_cost', 'z_reconfig_cost', 'vnf_reconfig_cost'} ;
                    if strcmpi(this.GET_BETA_METHOD, 'LeastSquare')
                        for i = 1:2
                            this.raw_cost.const.(field_names{i}) = ...
                                dot(this.(cost_entries{i}), this.Variables.(field_names{i})~=0);
                            this.raw_cost.linear.(field_names{i}) = ...
                                dot(this.(cost_entries{i}), this.Variables.(field_names{i}))/2;
                        end
                        this.raw_costv.const = ...
                            dot(this.vnf_reconfig_cost, this.Variables.v~=0);
                        this.raw_costv.linear = ...
                            dot(this.vnf_reconfig_cost, this.Variables.v)/4;
                        this.raw_costv.bInit = true;
                    else
                        for i =1:3
                            cost = dot(this.(cost_entries{i}), this.Variables.(field_names{i}))/2;
                            cost0 = dot(this.(cost_entries{i}), this.Variables.(field_names{i})~=0);
                            this.raw_beta.(field_names{i}) = cost0/cost;
                        end
                    end
                end
                
                stat = this.get_reconfig_stat();
                options = getstructfields(this.options, 'PricingPolicy', 'default', 'quadratic');
                options.bFinal = true;
                stat.Profit = this.getProfit(options);
                stat.ReconfigType = ReconfigType.Dimensioning;
                stat.ResourceCost = this.getSliceCost(options.PricingPolicy);   % _optimizeResourcePriceNew_ use 'quadratic-price'
                stat.FairIndex = (sum(this.FlowTable.Rate))^2/(this.NumberFlows*sum(this.FlowTable.Rate.^2));
                if this.options.ReconfigMethod >= ReconfigMethod.Dimconfig
                    this.sh_data.profits = stat.Profit;
                    this.sh_data.omegas = stat.Utilization;
                end
                if ~isempty(DEBUG) && DEBUG
                    disp(stat);
                end
                g_results = stat;       % The first event.
                this.b_dim = 0;
            else
                if this.options.ReconfigMethod >= ReconfigMethod.Dimconfig && ...
                        this.options.ReconfigMethod < ReconfigMethod.DimBaseline% && false
                    this.temp_vars = rmstructfields(this.temp_vars, {'x','z','tx','tz'});
                    this.update_reconfig_cost(this.sh_options.action);
                    new_opt = struct;
                    if this.options.ReconfigMethod == ReconfigMethod.DimconfigReserve
                        new_opt.bEnforceReserve = true;
                    end
                    this.fastReconfigure(this.sh_options.action, new_opt);
                end
            end
            [omega, sigma_o] = this.utilizationRatio();
            if this.options.ReconfigMethod >= ReconfigMethod.Dimconfig
                this.sh_options.sigma = sigma_o + var_sigma(omega);
            end
            %% Nest function: fiiting the variation of utilization
            % fitting y = ax^2+bx+c
            function s = var_sigma(omega, b) %#ok<INUSD>
                %                 A = [1 1 1;
                %                     0.9^2 0.9 1;
                %                     0.8^2 0.8 1;
                %                     0.5^2 0.5 1;
                %                     0 0 1];
                %                 if nargin == 1
                %                     if this.NumberFlows <= 60
                %                         b = [0.3; 0.2; 0.15; 0.1; 0.3];
                %                     elseif this.NumberFlows <= 120
                %                         b = [0.2; 0.14; 0.1; 0.05; 0.2];
                %                     elseif this.NumberFlows <= 180
                %                         b = [0.18; 0.12; 0.09; 0.04; 0.18];
                %                     elseif this.NumberFlows <= 240
                %                         b = [0.15; 0.1; 0.08; 0.03; 0.15];
                %                     else
                %                         b = [0.12; 0.09; 0.075; 0.025; 0.12];
                %                     end
                %                 end
                %                 c = A\b;
                %                 f = @(x) c(1)*x^2+c(2)*x+c(3);
                %                 s = f(omega);
                %%
                s = 0.2*sqrt(100/this.NumberFlows);
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
        %% Get link and node capacity
        % <getLinkCapacity>
        % <getNodeCapacity>
        % Due to resource reservation or reconfiguration cost constraint, the link/node
        % capcity might be larger than the demand.
        function c = getLinkCapacity(this, isfinal)
            if nargin == 1 || isfinal
                c = this.VirtualLinks.Capacity;
            else
                if this.invoke_method == 0
                    c = getLinkCapacity@Slice(this, isfinal);
                else
                    c = this.temp_vars.c;
                end
            end
        end
        
        function c = getNodeCapacity(this, isfinal) % isfinal=true by default
            if nargin == 1 || isfinal
                c = this.VirtualDataCenters.Capacity;
            else
                if this.invoke_method == 0
                    c = getNodeCapacity@Slice(this, isfinal);
                else
                    % since we do not reconfigure VNF capacity through fast slice reconfiguration,
                    % the sum of VNF capacity equals to the node capacity.
                    c = sum(reshape(this.temp_vars.v, ...
                        this.NumberDataCenters,this.NumberVNFs),2);
                end
            end
        end
        
        function c = get.VNFCapacity(this)
            if ~isfield(this.Variables, 'v') || isempty(this.Variables.v)
                warning('VNF capacity not set, set to VNF load.');
                c = this.getVNFCapacity;
                this.Variables.v = c;
            else
                c = this.Variables.v;
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
            switch pricing_policy
                case {'quadratic-price', 'quadratic'}
                    link_payment = this.fcnLinkPricing(link_price, link_load);
                    node_payment = this.fcnNodePricing(node_price, node_load);
                    cost = link_payment + node_payment;
                case 'linear'
                    cost = dot(link_price, link_load) + dot(node_price, node_load);
                otherwise
                    error('%s: invalid pricing policy', calledby);
            end
            if ~strcmpi(reconfig_cost_model, 'none')
                cost = cost + this.get_reconfig_cost(reconfig_cost_model, true);
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
        
        function b = getbeta(this)
            field_names = {'x','z','v'};
            if strcmpi(this.GET_BETA_METHOD, 'LeastSquare')
                for i=1:2
                    b.(field_names{i}) = sum([this.raw_cost.const.(field_names{i})].*[this.raw_cost.linear.(field_names{i})])./ ...
                        sum([this.raw_cost.linear.(field_names{i})].^2);
                end
                b.v = sum(this.raw_costv.const.*this.raw_costv.linear)./ ...
                    sum(this.raw_costv.linear.^2);
            else
                if strcmpi(this.GET_BETA_METHOD, 'ExponentiaLMovingAverage')
                    a0 = 0.25;
                    for i=1:2
                        %                 b.(field_names{i}) = mean(this.raw_beta.(field_names{i}));
                        n = length(this.raw_beta.(field_names{i}));
                        alpha = zeros(n,1);
                        alpha(2:n) = a0 * (1-a0).^(n-2:-1:0);
                        alpha(1) = (1-a0)^(n-1);
                        b.(field_names{i}) = dot(alpha, this.raw_beta.(field_names{i}));
                    end
                else
                    for i=1:2
                        b.(field_names{i}) = mean(this.raw_beta.(field_names{i}));
                    end
                end
                % v_nf >= sum_{p}{z_npf}
                n_path = nnz(this.I_node_path)/this.NumberDataCenters;
                b.v = b.z/n_path;
            end
        end
    end
    
    methods (Access = {?DynamicSlice, ?DynamicNetwork})
        %% calculate the difference of slice state
        % Arguments:
        %   # isfinal: If the slice is not in the final stage, the difference should be
        %     update with each call. Otherwise, the difference can be saved for later use.
        %     see also <Slice.isFinal>; 
        function ds = diffstate(this, isfinal)
            if isfinal && ~isempty(this.diff_state)
                if nargout == 1
                    ds = this.diff_state;
                end
                return;
            end            
            
            % At the beginging, |old_variables| equals to |Variables|.
            % |old_variables| is still valid after adding/removing flows.
            % |old_variables| have more elements than |this.Variables.x| when removing
            % flows.
            if length(this.old_variables.x) <= length(this.temp_vars.x)  % only ==
                if isfinal
                    new_x = this.Variables.x;
                else
                    new_x = this.temp_vars.x;
                end
                %% 
                % In the optimization, we ignore the reconfiguration cost of the new flow,
                % while in the results, we should consider the new flow's reconfiguration
                % cost. See also <update_reconfig_cost>. 
                ds.tI_node_path = this.I_node_path;
                % adding flow, |path| will increase, and |edge| might(might not) increase;
                ds.tI_edge_path = this.I_edge_path;
            else
                new_x = zeros(size(this.old_variables.x));
                if isfinal
                    new_x(~this.changed_index.x) = this.Variables.x;
                else
                    new_x(~this.changed_index.x) = this.temp_vars.x;
                end
                %% The removed flow has been considered with the reconfguration cost
                % The removed flow's reconfiguration cost is constant and does not appear
                % in the optimization, see also <update_reconfig_cost>.
                %
                % removing flow: |path| will decrease, and |edge| might(might not) decrease;
                ds.tI_node_path = this.old_state.I_node_path;
                ds.tI_edge_path = this.old_state.I_edge_path;
            end
            ds.diff_x = sparse(new_x-this.old_variables.x);
            mid_x = 1/2*(abs(new_x)+abs(this.old_variables.x));
            nz_index_x = mid_x~=0;
            ds.diff_x_norm = sparse(length(ds.diff_x),1);
            % Here, we assume that all tiny variables (smaller than NonzeroTolerance) have been
            % eleminated. So that the following operation can identify changes.
            ds.diff_x_norm(nz_index_x) = abs(ds.diff_x(nz_index_x)./mid_x(nz_index_x));

            if length(this.old_variables.z) <= length(this.temp_vars.z)
                if isfinal
                    new_z = this.Variables.z;
                else
                    new_z = this.temp_vars.z;
                end
            else
                new_z = zeros(size(this.old_variables.z));
                if isfinal                
                    new_z(~this.changed_index.z) = this.Variables.z;
                else
                    new_z(~this.changed_index.z) = this.temp_vars.z;
                end
            end
            ds.diff_z = sparse(new_z-this.old_variables.z);
            mid_z = 1/2*(abs(new_z)+abs(this.old_variables.z));
            nz_index_z = mid_z~=0;
            ds.diff_z_norm = sparse(length(ds.diff_z),1);
            ds.diff_z_norm(nz_index_z) = abs(ds.diff_z(nz_index_z)./mid_z(nz_index_z));
            %%%
            % Reconfiguration cost of VNF capacity.
            % No reconfiguration of VNF instance for 'fastconfig'.
            if isfield(this.temp_vars, 'v')
                if length(this.old_variables.v) <= length(this.temp_vars.v)
                    if isfinal
                        new_vnf_capacity = this.VNFCapacity(:);
                    else
                        new_vnf_capacity = this.temp_vars.v;
                    end
                else
                    new_vnf_capacity = zeros(length(this.old_variables.v),1);
                    if isfinal
                        new_vnf_capacity(~this.changed_index.v) = this.VNFCapacity(:);
                    else
                        new_vnf_capacity(~this.changed_index.v) = this.temp_vars.v;
                    end
                end
                ds.diff_v = new_vnf_capacity-this.old_variables.v(:);
                ds.mid_v = 1/2*(abs(new_vnf_capacity)+abs(this.old_variables.v(:)));
                nz_index_v = ds.mid_v~=0;
                ds.diff_v_norm = zeros(length(ds.diff_v),1);
                ds.diff_v_norm(nz_index_v) = abs(ds.diff_v(nz_index_v)./ds.mid_v(nz_index_v));
            end
            if isfinal
                this.diff_state = ds;
            end
        end
        
        % model = {'linear'|'const'|'none'}
        % Return: reconfiguration cost, number of reconfiguration, ration of
        % reconfigurations, and total number of variables.
        function [total_cost, reconfig_cost] = get_reconfig_cost(this, model, isfinal)
            if nargin <= 1
                model = 'const';
            end
            if strcmpi(model, 'const')
                isfinal = true;
            elseif strcmpi(model, 'linear') && nargin <= 2
                isfinal = false;
            end
            options = getstructfields(this.Parent.options, ...
                {'DiffNonzeroTolerance', 'NonzeroTolerance'});
            tol_vec = options.DiffNonzeroTolerance;
            s = this.diffstate(isfinal);
            
            if strcmpi(model, 'const')
                % logical array cannot be used as the first argument of dot.
                reconfig_cost.x = dot(this.x_reconfig_cost, s.diff_x_norm>tol_vec);
                reconfig_cost.z = dot(this.z_reconfig_cost, s.diff_z_norm>tol_vec);
                if isfield(s, 'diff_v_norm')
                    reconfig_cost.v = dot(this.vnf_reconfig_cost, s.diff_v_norm>tol_vec);
                end
            elseif strcmpi(model, 'linear')
                reconfig_cost.x = dot(abs(s.diff_x), this.x_reconfig_cost);
                reconfig_cost.z = dot(abs(s.diff_z), this.z_reconfig_cost);
                if isfield(s, 'diff_v')
                    reconfig_cost.v = dot(abs(s.diff_v), this.vnf_reconfig_cost);
                end
                % The linear cost coefficient (*_reconfig_cost) is the original
                % one, which is not scaled. In orde to compute the scaled cost as the
                % problem formulation, we additionally mutiply the scaler.
                if this.ENABLE_DYNAMIC_NORMALIZER
                    beta = this.getbeta();
                    reconfig_cost.x = beta.x * reconfig_cost.x;
                    reconfig_cost.z = beta.z * reconfig_cost.z;
                    if isfield(s, 'diff_v')
                        reconfig_cost.v = beta.v * reconfig_cost.v;
                    end
                else
                    reconfig_cost.x = this.options.ReconfigScaler * reconfig_cost.x;
                    reconfig_cost.z = this.options.ReconfigScaler * reconfig_cost.z;
                    if isfield(s, 'diff_v')
                        reconfig_cost.v = this.options.ReconfigScaler * reconfig_cost.v;
                    end
                end
            else
                error('%s: invalid cost model <%s>.', calledby, model);
            end
            total_cost = reconfig_cost.x + reconfig_cost.z;
            if isfield(reconfig_cost, 'v')
                total_cost = total_cost + reconfig_cost.v;
            end
        end
        
        function stat = get_reconfig_stat(this, stat_names)
            if nargin == 1
                stat_names = {'All'};
            elseif ischar(stat_names)
                stat_names = {stat_names};
            end
            options = getstructfields(this.Parent.options, ...
                {'DiffNonzeroTolerance', 'NonzeroTolerance'});
            s = this.diffstate(true);
            
            stat = table;
            tol_vec = options.DiffNonzeroTolerance;
            for i = 1:length(stat_names)
                %% Time
                if contains(stat_names{i},{'All', 'Time'},'IgnoreCase',true)
                    stat.Time = this.time.Current;
                end
                if contains(stat_names{i},{'All', 'Utilization'},'IgnoreCase',true)
                    stat.Utilization = this.utilizationRatio();
                end
                %% Reconfiguration Cost
                if contains(stat_names{i},{'All', 'Cost'},'IgnoreCase',true)
                    stat.Cost = this.get_reconfig_cost('const');
                end
                if contains(stat_names{i},{'All', 'LinearCost'},'IgnoreCase',true)
                    stat.LinearCost = this.get_reconfig_cost('linear', true);
                end
                %% Number of Variables
                if contains(stat_names{i},{'All', 'NumberVariables'},'IgnoreCase',true)
                    stat.NumberVariables = ...
                        nnz(s.tI_edge_path)+nnz(s.tI_node_path)*this.NumberVNFs;
                    if isfield(s, 'mid_v')
                        stat.NumberVariables = stat.NumberVariables + nnz(s.mid_v);
                    end
                end
                %% Number of reconfigurations between current and previous solution
                % NOTE: resource allocation and release will not take place at the same
                % time.
                %   Number of edge variables, VNF assignment variables
                %   the reconfiguration of VNF instances.
                if contains(stat_names{i},{'All', 'NumberReconfigVariables'},'IgnoreCase',true)
                    paths_length = sum(s.tI_edge_path,1);
                    stat.NumberReconfigVariables = nnz(s.diff_z_norm>tol_vec) ...
                        + dot(paths_length,s.diff_x_norm>tol_vec);
                    if isfield(s, 'diff_v_norm')
                        stat.NumberReconfigVariables = stat.NumberReconfigVariables + ...
                            nnz(s.diff_v_norm>tol_vec);
                    end
                end
                %% Number of Flows
                if contains(stat_names{i},{'All', 'NumberFlows'},'IgnoreCase',true)
                    % For convenience of comparison, we store the number of flows including the
                    % removed one.
                    stat.NumberFlows = max(this.NumberFlows, height(this.old_state.flow_table));
                end
                %% Number of Reconfigured Flows
                if contains(stat_names{i},{'All', 'NumberReconfigFlows'},'IgnoreCase',true)
                    np = numel(s.diff_x_norm);
                    nv = this.NumberVNFs;
                    nn = numel(s.diff_z_norm)/(nv*np);
                    diff_path = (s.diff_x_norm>tol_vec)' + ...
                        sum(sum(reshape(full(s.diff_z_norm>tol_vec),nn,np,nv),3),1);
                    if length(this.path_owner) <= length(this.old_state.path_owner)
                        stat.NumberReconfigFlows = ...
                            numel(unique(this.old_state.path_owner(diff_path~=0)));
                    else
                        stat.NumberReconfigFlows = ...
                            numel(unique(this.path_owner(diff_path~=0)));
                    end
                end
                if contains(stat_names{i},{'All', 'ResourceCost'},'IgnoreCase',true)
                    options = structmerge(options, ...
                        getstructfields(this.options, 'PricingPolicy', 'default', 'quadratic'));
                    stat.ResourceCost = this.getSliceCost(options.PricingPolicy);
                end
                if contains(stat_names{i},{'All', 'FairIndex'},'IgnoreCase',true)
                    stat.FairIndex = (sum(this.FlowTable.Rate))^2/...
                        (this.NumberFlows*sum(this.FlowTable.Rate.^2));
                end
                if contains(stat_names{i},{'All', 'Interval'},'IgnoreCase',true)
                    if this.b_dim
                        stat.Interval = this.time.DimensionIntervalModified;
                    else
                        stat.Interval = this.time.ConfigureInterval;
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
        
        %% Convert the optimizatio results to temporary variables
        % temp variables might include: x,z,v,tx,tz,tv,c (see also
        % <Dynamic.priceOptimalFlowRate>);
        function convert(this, x, ~)
            NP = this.NumberPaths;
            this.temp_vars.x = x(1:NP);
            this.temp_vars.z = x((NP+1):this.num_vars);
            offset = this.num_vars;
            this.temp_vars.v = x(offset+(1:this.num_varv));
            offset = offset + this.num_varv;
            if this.invoke_method >= 2
                this.temp_vars.tx = x(offset+(1:NP));
                this.temp_vars.tz = x(offset+((NP+1):this.num_vars));
                offset = offset + this.num_vars;
                this.temp_vars.tv = x(offset+(1:this.num_varv));
                offset = offset + this.num_varv;
            end
            this.temp_vars.c = x(offset+(1:this.NumberVirtualLinks));
            %             if nargin >= 3
            %                 this.Variables = getstructfields(this.temp_vars, {'x','z','v','c'}, 'ignore');
            %             end
        end
        function [profit, cost] = handle_zero_flow(this, new_opts)
            cost = this.getSliceCost(new_opts.PricingPolicy);
            profit = -cost;
            this.temp_vars.x = [];
            this.temp_vars.z = [];
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
        
        %%
        % dynamic adjust the normalizer for L1-approximation of the reconfiguration cost.
        % This method is used when 'ENABLE_DYNAMIC_NORMALIZER' is set to true.
        function postl1normalizer(this)
            [~, cost] = this.get_reconfig_cost('const');
            [~, linear_cost] = this.get_reconfig_cost('linear', true);
            field_names = {'const', 'linear'};
            b = this.getbeta();
            if strcmpi(this.GET_BETA_METHOD, 'LeastSquare')
                arr_cost = [cost, linear_cost];
                fn = fieldnames(arr_cost(2));
                for j = 1:length(fn)
                    arr_cost(2).(fn{j}) = arr_cost(2).(fn{j})/b.(fn{j});
                end
                for i = 1:2
                    if length(this.raw_cost.(field_names{i})) < this.NUM_MEAN_BETA
                        this.raw_cost.(field_names{i})(end+1) = ...
                            getstructfields(arr_cost(i), {'x','z'});
                    else
                        this.raw_cost.(field_names{i}) = ...
                            [this.raw_cost.(field_names{i})(2:end), ...
                            getstructfields(arr_cost(i), {'x','z'})];
                    end
                end
                if isfield(arr_cost, 'v')   
                    if arr_cost(1).v < 10^-3 || arr_cost(2).v < 10^-3
                        return;
                    end
                    if this.raw_costv.bInit
                        this.raw_costv.const = arr_cost(1).v;
                        this.raw_costv.linear = arr_cost(2).v;
                        this.raw_costv.bInit = false;
                        return;
                    end
                    for i = 1:2
                        if length(this.raw_costv.(field_names{i})) < this.NUM_MEAN_BETA
                            this.raw_costv.(field_names{i})(end+1) = arr_cost(i).v;
                        else
                            this.raw_costv.(field_names{i}) = ...
                                [this.raw_costv.(field_names{i})(2:end), arr_cost(i).v];
                        end
                    end
                end                
            else
                if this.b_dim
                    field_names = {'x','z','v'};
                else
                    field_names = {'x','z'};
                end
                for i = 1:length(field_names)
                    if isfield(cost, field_names{i}) && isfield(linear_cost, field_names{i})
                        if cost.(field_names{i})<10^-6
                            %                         beta.(field_names{i}) = this.raw_beta.(field_names{i})(end)/2;
                            continue;
                        else
                            f = cost.(field_names{i})/linear_cost.(field_names{i});
                            if f > 2
                                f = 2;
                            elseif f<1/2
                                f = 1/2;
                            end
                            beta.(field_names{i}) = b.(field_names{i})*(1+f)/2;
                        end
                        if length(this.raw_beta.(field_names{i})) < this.NUM_MEAN_BETA
                            this.raw_beta.(field_names{i})(end+1) = beta.(field_names{i});
                        else
                            this.raw_beta.(field_names{i}) = ...
                                [this.raw_beta.(field_names{i})(2:end), beta.(field_names{i})];
                        end
                    end
                end
            end

            %             function [f,g] = fcn_optbeta(beta,cc,ca)
            %                 f = (norm(beta*ca-cc,2))^2;
            %                 g = sum(2*(beta*ca-cc).*ca);
            %             end
            %             function h = hess_optbeta(beta, lambda, ca) %#ok<INUSL>
            %                 h = sum(2*ca.^2);
            %             end
        end
    end
    
    
    methods (Access = {?Slice, ?CloudNetwork, ?SliceFlowEventDispatcher})
        [utility, node_load, link_load] = priceOptimalFlowRate(this, x0, new_opts);
    end
    methods (Access = {?Slice, ?CloudNetwork})
        %%
        % |new_opts|:
        % * *FixedCost*: used when slice's resource amount and resource prices are fixed,
        % resource cost is fixed. When this options is specified, the base method
        % <Slice.optimalFlowRate> returns profit without resource cost.
        % See also <DynamicSlice.fcnSocialWelfare>,<DynamicSlice.optimize>.
        function [profit, cost] = optimalFlowRate( this, new_opts )
            if nargin <= 1
                new_opts = struct;
            end
            theta0 = 0.9999;
            this.VirtualDataCenters.Capacity = theta0*this.VirtualDataCenters.Capacity;
            this.VirtualLinks.Capacity = theta0*this.VirtualLinks.Capacity;
            if ~isfield(new_opts, 'CostModel') || ~strcmpi(new_opts.CostModel, 'fixcost')
                [profit,cost] = optimalFlowRate@Slice( this, new_opts );
            else 
                if this.NumberFlows == 0
                    [profit, cost] = this.handle_zero_flow(new_opts);
                else
                    [profit, cost] = optimalFlowRate@Slice( this, new_opts );
                    if isfield(new_opts, 'CostModel') && strcmpi(new_opts.CostModel, 'fixcost')
                        profit = profit - cost;
                    end
                    if nargout >= 1
                        % When output argument specified, we finalize the VNF capacity.
                        % After reconfiguration VNF instance capcity has changed.
                        v = this.getVNFCapacity;
                        this.Variables.v = v(:);
                    end
                end   
                profit = profit - this.get_reconfig_cost('const');
            end
            this.VirtualDataCenters.Capacity = this.old_state.node_capacity;
            this.VirtualLinks.Capacity = this.old_state.link_capacity;
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
            this.max_flow_rate = [this.max_flow_rate; inf*ones(num_new_flows,1)];
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
            this.topts.old_variables_x = this.old_variables.x;
            this.topts.old_variables_z = this.old_variables.z;
            this.update_reconfig_cost('add');
            this.sh_options.action = 'add';
            [ef, ~] = this.executeMethod('add');
            %% Failure handling
            % Zero-rate: now we leave it unhandled, the zero-rate flow will stay in
            % the slice. It may be allocated with resource at later stage. On the
            % other hand, we can remove the zero-rate flows from the slice (reject).
            if ef < 0
                fidx = [];
                this.set_state();
            else
                this.max_flow_rate(this.max_flow_rate==inf) = 0;
                this.max_flow_rate = max(this.FlowTable.Rate, this.max_flow_rate);
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
            this.topts.old_variables_x = this.old_variables.x(~this.changed_index.x);
            this.topts.old_variables_z = this.old_variables.z(~this.changed_index.z);
            this.update_reconfig_cost('remove');
            % Before performing optimization, there is not situation that VNF instances
            % will be removed. So there is no need to copy the |old_variable.v| like x/z.
            
            %% Update flow table and local path id.
            % After removing flows, re-allocate local identifier to subsequent flows (flow
            % ID) and paths (path ID).
            pid = this.FlowTable{flow_index(1), 'Paths'}.paths{1}.local_id-1;
            fidt = flow_index(1):(this.NumberFlows-length(flow_index));
            this.FlowTable(fidx,:) = [];
            this.max_flow_rate(fidx) = [];
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
            
            this.sh_options.action = 'remove';
            [ef,~] = this.executeMethod('remove');
            if ef < 0
                this.set_state();
            else
                this.max_flow_rate = max(this.FlowTable.Rate, ...
                    this.max_flow_rate);
            end
        end
        
        
        % |parameters|: include fields x0, As, bs, Aeq, beq, lbs, ub, lb;
        % |options|: include fields fmincon_opt, CostModel, PricingPolicy (if
        %       CostModel='fixcost'), num_orig_vars;
        function [x, fval] = optimize(this, params, options)
            if isfield(options, 'CostModel') && strcmpi(options.CostModel, 'fixcost')
                [xs, fval, exitflag, output] = ...
                    fmincon(@(x)DynamicSlice.fcnSocialWelfare(x, this, ...
                    getstructfields(options, 'CostModel')), ...
                    params.x0, params.As, params.bs, params.Aeq, params.beq, ...
                    params.lb, params.ub, [], options.fmincon_opt);
                this.interpretExitflag(exitflag, output.message);
                if isfield(options, 'bCompact') && options.bCompact
                    x = zeros(options.num_orig_vars, 1);
                    x(this.I_active_variables) = xs;
                else
                    x = xs;
                end
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
        %% TODO postProcessing with reconfiguration examining
        % Fist check constraint violation, then post-process the reconfiguration.
        %   flow processing requirement;
        %   link capacity constraint;
        %   VNF capacity constraint;
        % Resolve constraint violation by down-scaling, which might induce more
        % reconfiguration. 
        %    
        % As a result, we currently ignore the violation, only focus on 
        % how to remove unacessary reconfigurations.
        function [b, vars] = postProcessing(this)
%             % [TODO] May need post processing for VNF capacity constraint;            
            %% post process reconfiguration
            postProcessing@Slice(this);
            options = getstructfields(this.Parent.options, ...
                {'DiffNonzeroTolerance', 'NonzeroTolerance', 'ConstraintTolerance'});
            if isfield(this.temp_vars, 'v')
                tol_zero = this.Parent.options.NonzeroTolerance;
                var_v = this.temp_vars.v;
                var_v(var_v<tol_zero*max(var_v)) = 0;
                this.Variables.v = var_v;
            end
            %% Another simple method
            % directly recover those little changes, which might cause constraint violation
            tol_vec = options.DiffNonzeroTolerance;
            if ~this.b_derive_vnf
                s = this.diffstate(true);
                if length(s.diff_x) > length(this.Variables.x)
                    s.diff_x_norm = s.diff_x_norm(~this.changed_index.x);
                    s.diff_z_norm = s.diff_z_norm(~this.changed_index.z);
                    old_x = this.old_variables.x(~this.changed_index.x);
                    old_z = this.old_variables.z(~this.changed_index.z);
                else
                    old_x = this.old_variables.x;
                    old_z = this.old_variables.z;
                end
                b_diss_x = s.diff_x_norm < tol_vec;
                this.Variables.x(b_diss_x) = old_x(b_diss_x);
                b_diss_z =  s.diff_z_norm < tol_vec;
                this.Variables.z(b_diss_z) = old_z(b_diss_z);
                if isfield(this.temp_vars, 'v')
                    if length(s.diff_v) > length(var_v)
                        s.diff_v_norm = s.diff_v_norm(~this.changed_index.v);
                        old_v = this.old_variables.v(~this.changed_index.v);
                    else
                        old_v = this.old_variables.v;
                    end
                    b_diss_v =  s.diff_v_norm < tol_vec;
                    var_v(b_diss_v) = old_v(b_diss_v);
                    this.Variables.v = var_v;
                end
                this.diff_state = struct([]);
                b = true; vars = this.Variables;
                return;
                %{
                find(this.As_res * [this.Variables.x; this.Variables.z]>0,1)
                res = this.I_edge_path*this.Variables.x-this.VirtualLinks.Capacity;
                find(res>0)
                disp(res)
                this.Hdiag*this.Variables.z-this.Variables.v>0
                %}                
            end
            if ~this.b_derive_vnf
                %% do processing to discard minor changes.
                % When performing slice dimensioning, this only be performed after the
                % final iteration. The intermediate iteration only use the linear
                % approximated cost.
                s = this.diffstate(true);
                if length(s.diff_x) > length(this.Variables.x)
                    s.diff_x = s.diff_x(~this.changed_index.x);
                    s.diff_x_norm = s.diff_x_norm(~this.changed_index.x);
                    s.diff_z = s.diff_z(~this.changed_index.z);
                    s.diff_z_norm = s.diff_z_norm(~this.changed_index.z);
                    old_x = this.old_variables.x(~this.changed_index.x);
                    old_z = this.old_variables.z(~this.changed_index.z);
                else
                    old_x = this.old_variables.x;
                    old_z = this.old_variables.z;
                end
                restore_x = this.Variables.x;
                restore_z = this.Variables.z;
                if isfield(this.temp_vars, 'v')
                    restore_v = this.temp_vars.v;
                    if length(s.diff_v) > length(restore_v)
                        s.diff_v = s.diff_v(~this.changed_index.v);
                        s.diff_v_norm = s.diff_v_norm(~this.changed_index.v);
                        s.mid_v = s.mid_v(~this.changed_index.v);
                        old_v = this.old_variables.v(~this.changed_index.v);
                    else
                        old_v = this.old_variables.v;
                    end
                    b_diss_v = (s.diff_v<0) & (s.diff_v_norm<tol_vec);
                    restore_v(b_diss_v) = old_v(b_diss_v);
                    s.diff_v_norm(b_diss_v) = 0;
                    s.diff_v(b_diss_v) = 0;
                else
                    restore_v = this.Variables.v;
                end
%                 b_diss_z = (s.diff_z_norm>0) && (s.diff_z_norm<tol_vec) && s.diff_z>0;
%                 if length(this.old_variables.x) <= length(this.temp_vars.x)
%                     this.Variables.x(b_diss_x) = this.old_variables.x(b_diss_x);
%                 else
%                     
%                 end
                % separately process the two part of flows, firstly process the path's
                % that might release some resource could improve the possibility to accept
                % the second part of flows, which need a little more resources.
                %%
                % first case: |x| does not change or the increase amount less than
                % |tol_vec|, try to restore these changes. 
                b_diss_x1 = s.diff_x>=0;   % x might not change while z might still change
                Ndc = this.NumberDataCenters;
                Np = this.NumberPaths;
                Nvnf = this.NumberVNFs;
                af = this.Parent.VNFTable{this.VNFList, 'ProcessEfficiency'};
                %% 1-1
                for p = (find(b_diss_x1))'
                    if s.diff_x_norm(p)<tol_vec
                        % If x has significant change, usually z also has significant
                        % change. But we cannot claim that all z components has
                        % significant changes.
                        restore_x(p) = old_x(p);
                        s.diff_x_norm(p) = 0;
                        s.diff_x(p) = 0;
                    end
                    check_z(1);
                end
                %% 1-2
                delta_vnf = restore_v - this.Hdiag*restore_z;  % residual VNF capacity
                for p = (find(b_diss_x1))'
                    check_z(2);
                end
                %% 2-1
                % second case: |x|'s decrease amount less than |tol_vec|, try to restore
                % these changes.  
                b_diss_x2 = s.diff_x<0;
                % VirtualLinks.Capacity has been updated in <finalize>.
                if isfield(this.temp_vars, 'c')
                    link_capacity = this.temp_vars.c;
                else
                    link_capacity = this.VirtualLinks.Capacity;
                end
                delta_link = link_capacity - this.I_edge_path*restore_x; % residual link capacity
                for p = (find(b_diss_x2))'
                    if s.diff_x_norm(p)<tol_vec
                        temp_delta_link = delta_link+this.I_edge_path(:,p).*s.diff_x(p);
                        if isempty(find(temp_delta_link<0,1)) % check link capacity constraint
                            temp_x = restore_x(p);
                            restore_x(p) = old_x(p); % restore: up-scaling
                            tf = check_z(1);
                            if tf == false % up-scaling restore failed, recover current value.
                                restore_x(p) = temp_x;
                            else
                                delta_link = temp_delta_link;
                                s.diff_x_norm(p) = 0;
                                s.diff_x(p) = 0;
                            end
                        end
                    else
                        check_z(1);
                    end
                end
                %% 2-2
                for p = (find(b_diss_x2))'
                    if s.diff_x_norm(p)<tol_vec
                        temp_delta_link = delta_link+this.I_edge_path(:,p).*s.diff_x(p);
                        if isempty(find(temp_delta_link<0,1)) % check link capacity constraint
                            temp_x = restore_x(p);
                            restore_x(p) = old_x(p); % restore: up-scaling
                            tf = check_z(2);
                            if tf == false % up-scaling restore failed, recover current value.
                                restore_x(p) = temp_x;
                            else
                                delta_link = temp_delta_link;
                                s.diff_x_norm(p) = 0;
                                s.diff_x(p) = 0;
                            end
                        end
                    else
                        check_z(2);
                    end
                end
                %% 1-3
                for p = (find(b_diss_x1))'
                    check_z(3);
                end
                %% 2-3
                for p = (find(b_diss_x2))'
                    if s.diff_x_norm(p)<tol_vec
                        temp_delta_link = delta_link+this.I_edge_path(:,p).*s.diff_x(p);
                        if isempty(find(temp_delta_link<0,1)) % check link capacity constraint
                            temp_x = restore_x(p);
                            restore_x(p) = old_x(p); % restore: up-scaling
                            tf = check_z(3);
                            if tf == false % up-scaling restore failed, recover current value.
                                restore_x(p) = temp_x;
                            else
                                delta_link = temp_delta_link;
                                s.diff_x_norm(p) = 0;
                                s.diff_x(p) = 0;
                            end
                        end
                    else
                        check_z(3);
                    end
                end
                %% reover VNF capacity 
                % must be performed after the flow assignement hast been done. Otherwise,
                % there is no chance to restore VNF capacity. Due to the reconfiguration
                % cost constraint, there is no redundant reconfiguration.
                % Instead, after the flow has been restored, there might be redundancy of
                % VNF capacity.
                if isfield(this.temp_vars, 'v')
                    b_diss_v = (s.diff_v>0) & (s.diff_v_norm<tol_vec);
                    idv = this.Hdiag*restore_z<=old_v & b_diss_v; % if the old value satisfy the capacity constraint
                    restore_v(idv) = old_v(idv);
                    s.diff_v_norm(idv) = 0;
                    s.diff_v(idv) = 0;
                    
%                     delta_dc = full(sum(reshape(this.Variables.v, Ndc, Nvnf),2)) -...
%                         full(sum(reshape(restore_v, Ndc, Nvnf),2)); % residual node capacity
%                     % all nodes can be check in parallel, with VNF varying
%                     idv = 1:Ndc;
%                     for v = 1:Nvnf
%                         tv = restore_v(idv);
%                         old_tv = old_v(idv);
%                         idtv = (s.diff_v_norm(idv)>0) & (s.diff_v_norm(idv)<tol_vec) & s.diff_v(idv)<0;
%                         tv(idtv) = old_tv(idtv);
%                         temp_delta_dc = delta_dc + restore_v(idv) - tv;
%                         fx = temp_delta_dc >= 0;
%                         restore_v(idv(fx)) = tv(idv(fx));
%                         delta_dc(fx) = temp_delta_dc(fx);
%                         s.diff_v_norm(idv(fx)) = 0;
%                         s.diff_v(idv(fx)) = 0;
%                         idv = idv + Ndc;
%                     end                    
                    this.Variables.v = restore_v;
                end
                this.Variables.x = restore_x;
                this.Variables.z = restore_z;
                this.diff_state = struct([]);
            end
            b = true; vars = this.Variables;
            
            %% nest function: check z
            function tf = check_z(t)
                tf = true;
                b_reset = false(Nvnf,1);
                temp_z = zeros(Ndc, Nvnf);
                temp_sz = zeros(Ndc, Nvnf);
                temp_szn = zeros(Ndc, Nvnf);
                delta_sum = zeros(Nvnf,1);
                if t > 1
                    temp_delta_vnf = zeros(Ndc, Nvnf);
                end
                
                idz = (1:Ndc) + (p-1)*Ndc;
                idn = 1:Ndc;                
                for k = 1:Nvnf
                    tz = restore_z(idz);
                    temp_z(:,k) = tz;
                    temp_sz(:,k) = s.diff_z(idz);
                    temp_szn(:,k) = s.diff_z_norm(idz);
                    if t > 1
                        temp_delta_vnf(:,k) = delta_vnf(idn);
                    end
                    old_tz = old_z(idz);
                    switch t
                        case 1
                            % try restore z(:,p,f) that should be restored by down-scaling;
                            % make space for other requests.
                            idtz = (s.diff_z(idz)>0) & (s.diff_z_norm(idz)<tol_vec);
                        case 2
                            %%% try restore z(:,p,f) that includes both up-sclaing and down-scaling;
                            idtz = (s.diff_z_norm(idz)>0) & (s.diff_z_norm(idz)<tol_vec);
                        case 3
                            %%% try restore z(:,p,f) that only includes up-scaling;
                            idtz = (s.diff_z(idz)<0) & (s.diff_z_norm(idz)<tol_vec);
                    end
                    if isempty(find(idtz,1))
                        if nargout == 1 && dot(this.I_node_path(:,p), tz) < af(k)*restore_x(p)
                            % when x increase, but z keeps unchange, so we need to check the processing
                            % constraint.
                            tf = false;
                        end
                        idz = idz + Np*Ndc;
                        idn = idn + Ndc;
                        continue;
                    end
                    tz(idtz) = old_tz(idtz);  % recover: down-scaling | down/up-scaling | up-scaling
                    if t == 2 || t == 3
                        delta_tz = restore_z(idz) - tz; % current - past
                        t_delta_vnf = delta_vnf(idn)+this.I_node_path(:,p).*delta_tz;
                    end
                    if t == 1
                        if dot(this.I_node_path(:,p), tz) >= af(k)*restore_x(p) % check processing constraint
                            restore_z(idz) = tz;
                            s.diff_z_norm(idz(idtz)) = 0;
                            s.diff_z(idz(idtz)) = 0;
                            delta_sum(k) = 0;
                        else
                            b_reset(k) = true;
                        end
                    end
                    if t == 2
                        if dot(this.I_node_path(:,p), tz) >= af(k)*restore_x(p) &&...
                                isempty(find(t_delta_vnf<0,1))
                            %% check both processing constraint and VNF capacity constraint
                            restore_z(idz) = tz; % accept old value or keep current value.
                            s.diff_z_norm(idz(idtz)) = 0;
                            s.diff_z(idz(idtz)) = 0;
                            delta_vnf(idn) = t_delta_vnf; % update residual capacity
                            delta_sum(k) = sum(delta_tz);
                        else
                            b_reset(k) = true;
                        end
                    end
                    if t == 3 
                        if (nargout == 0 && isempty(find(t_delta_vnf<0,1))) || (nargout == 1 && isempty(find(t_delta_vnf<0,1)) && dot(this.I_node_path(:,p), tz) >= af(k)*restore_x(p))
                            restore_z(idz) = tz; % accept old value or keep current value.
                            s.diff_z_norm(idz(idtz)) = 0;
                            s.diff_z(idz(idtz)) = 0;
                            delta_vnf(idn) = t_delta_vnf; % update residual capacity
                            delta_sum(k) = sum(delta_tz);
                        else
                            b_reset(k) = true;
                        end
                    end
                    idz = idz + Np*Ndc;
                    idn = idn + Ndc;
                end
                if ~isempty(find(b_reset,1))
                    tf = false;
                    idz = (1:Ndc) + (p-1)*Ndc;
                    idn = 1:Ndc;
                    for k = 1:Nvnf
                        if ~b_reset(k) && delta_sum(k) < 0
                            restore_z(idz) = temp_z(:,k);
                            s.diff_z(idz) = temp_sz(:,k);
                            s.diff_z_norm(idz) = temp_szn(:,k);
                            if t > 1
                                delta_vnf(idn) = temp_delta_vnf(:,k);
                            end
                        end
                        idz = idz + Np*Ndc;
                        idn = idn + Ndc;
                    end
                end
            end
        end
        
        % update reconfiguration cost with scaler.
        % If dimensioning is applied, the orginal cost is firstly updated.
        function update_reconfig_costv(this)
            this.vnf_reconfig_cost = this.Parent.options.VNFReconfigCoefficient * ...
                repmat(this.VirtualDataCenters.ReconfigCost, this.NumberVNFs, 1);
            if this.ENABLE_DYNAMIC_NORMALIZER
                beta = this.getbeta();
                this.topts.vnf_reconfig_cost = beta.v*this.vnf_reconfig_cost;
            else
                this.topts.vnf_reconfig_cost = ...
                    this.options.ReconfigScaler*this.vnf_reconfig_cost;
            end
        end
        
        % update reconfiguration cost with scaler.
        function update_reconfig_cost(this, action, bDim)
            %% Period re-dimensioing
            % The number of virtual nodes/links is not change in the slice, as well as the
            % number of VNF instances.
            % When performing re-dimensioning, the reconfiguration cost is larger than
            % that of fast reconfiguration, so we need to update the reconfiguration cost.
            if nargin >= 3 && bDim
                this.Parent.updateRedimensionCost(this);
                this.update_reconfig_costv();
            end
            if strcmpi(action, 'add')
                this.x_reconfig_cost = (this.I_edge_path)' * this.VirtualLinks.ReconfigCost;
                this.z_reconfig_cost = repmat(this.VirtualDataCenters.ReconfigCost, ...
                    this.NumberPaths*this.NumberVNFs, 1);
                %% Ignore the reconfiguration cost when adding flow
                % The newly added flow's reconfiguration cost set to zero.
                % This helps the new flow be admitted into the slice, when there is
                % resource while the reconfiguration cost of the new flow is too high.
                % This can be explained as: the reconfiguration cost of the first time can
                % be treated as that it is distrbuted in the whole life time of the flow.
                % As long as the flow's lifetime is long enough, the first time cost could
                % be ignored.
                %
                % In the output, we still count the reconfiguration cost for the new flow,
                % see also <diffstate> and <get_reconfig_cost>.
                this.x_reconfig_cost(this.changed_index.x) = 0;
                this.z_reconfig_cost(this.changed_index.z) = 0;
                if this.ENABLE_DYNAMIC_NORMALIZER
                    beta =this.getbeta();
                    this.topts.x_reconfig_cost = beta.x*this.x_reconfig_cost;
                    this.topts.z_reconfig_cost = beta.z*this.z_reconfig_cost;
                else
                    this.topts.x_reconfig_cost = this.options.ReconfigScaler*this.x_reconfig_cost;
                    this.topts.z_reconfig_cost = this.options.ReconfigScaler*this.z_reconfig_cost;
                end
            else
                %%
                % Since we need maintain information of the deleted variables, we should use
                % the old state to calculate the information.
                this.x_reconfig_cost = ...
                    (this.old_state.I_edge_path)' * this.VirtualLinks.ReconfigCost;
                old_num_paths = length(this.old_state.path_owner);
                this.z_reconfig_cost = repmat(this.VirtualDataCenters.ReconfigCost, ...
                    old_num_paths*this.NumberVNFs, 1);
                %%
                % The removed flow's reconfiguration cost is not considered in the
                % optimization, which is constant.
                if this.ENABLE_DYNAMIC_NORMALIZER
                    beta =this.getbeta();
                    this.topts.x_reconfig_cost = ...
                        beta.x*this.x_reconfig_cost(~this.changed_index.x);
                    this.topts.z_reconfig_cost = ...
                        beta.z*this.z_reconfig_cost(~this.changed_index.z);
                else
                    this.topts.x_reconfig_cost = ...
                        this.options.ReconfigScaler*this.x_reconfig_cost(~this.changed_index.x);
                    this.topts.z_reconfig_cost = ...
                        this.options.ReconfigScaler*this.z_reconfig_cost(~this.changed_index.z);
                end
            end
        end
    end
    
    methods(Static, Access = protected)
        [profit, grad]= fcnProfitReconfigureSlicing(vars, slice, options);
        [profit, grad]= fcnProfitReserveSlicing(vars, slice, options);
        [profit, grad] = fcnSocialWelfare(vars, slice, options);
        hs = hessSlicing(var_x, lambda, slice, options);
        hs = hessInitialSlicing(vars, lambda, slice, options);
        %%
        % The resource consumption cost is constant, since the slice will not
        % release/request resouces to/from the substrate network, and the resource prices
        % during the reconfiguration does not change. Therefore, the resource consumption
        % cost it not counted in the objective function. The objective function contains
        % the user utility and the reconfiguration cost.
        %
        % options:
        %   # |num_varx|,|num_varz|,|num_varv|: the number of main variables. For compact
        %     mode (bCompact=true), |num_varz| is less than the original number of variables |z|.
        %   # |bFinal|: set to 'true' to output the original objective function value.
        %
        % NOTE: hession matrix of the objective is only has non-zero values for x
        % (super-linear). See also <fcnHessian>.
        function [profit, grad] = fcnFastConfigProfit(vars, this, options)
            num_vars = length(vars);
            num_basic_vars = options.num_varx + options.num_varz;
            % |tx|,|tz| and |tv| are auxilliary variables to transform L1 norm to
            % linear constraints and objective part. The number of auxiliary variables
            % eqauls to the number of the main variables.
            num_var_tx = length(this.topts.x_reconfig_cost);
            num_var_tz = length(this.topts.z_reconfig_cost);
            if isfield(options, 'num_varv')  % for _fastReconfigure2_, including |v|, and |tv|
                var_v_index = num_basic_vars + (1:options.num_varv);
                num_basic_vars = num_basic_vars + options.num_varv;
                num_var_tv = length(this.topts.vnf_reconfig_cost);
                var_tv_index = num_basic_vars + num_var_tx + num_var_tz + (1:num_var_tv);
                var_tv = vars(var_tv_index);
            end
            var_x = vars(1:options.num_varx);
            var_tx_index = num_basic_vars+(1:num_var_tx); 
            var_tx = vars(var_tx_index);
            var_tz_index = num_basic_vars+num_var_tx+(1:num_var_tz);
            var_tz = vars(var_tz_index);
            flow_rate = this.getFlowRate(var_x);
            %% objective value
            profit = -this.weight*sum(fcnUtility(flow_rate));
            %%%
            % calculate reconfiguration cost by |tx| and |tz| (L1 approximation),
            % which is similar to calculating by |x-x0| and |z-z0|.
            profit = profit + dot(var_tx, this.topts.x_reconfig_cost) + ...
                dot(var_tz, this.topts.z_reconfig_cost);
            if isfield(options, 'num_varv') % for _fastConfigure2_
                profit = profit + dot(var_tv, this.topts.vnf_reconfig_cost);
            end
            % If there is only one output argument, return the real profit (positive)
            if isfield(options, 'bFinal') && options.bFinal
                profit = -profit;
            else
                %% Gradient
                % The partial derivatives are computed by dividing the variables into four
                % parts, i.e., $x,z,t_x,t_z$. Since |z| does not appear in objective
                % function, the corrsponding derivatives is zero.
                % The partial derivatives of x
                grad = spalloc(num_vars, 1, num_vars-options.num_varz);
                for p = 1:options.num_varx
                    i = this.path_owner(p);
                    grad(p) = -this.weight/(1+this.I_flow_path(i,:)*var_x); %#ok<SPRIX>
                end
                %%%
                % The partial derivatives of t_x is the vector |x_reconfig_cost|;
                % The partial derivatives of t_z is duplicaton of |z_reconfig_cost|,
                % since the node reconfiguration cost is only depend on node.
                grad(var_tx_index) = this.topts.x_reconfig_cost;
                grad(var_tz_index) = this.topts.z_reconfig_cost;
                if isfield(options, 'num_varv')
                    grad(var_v_index) = this.topts.vnf_reconfig_cost;
                end
            end
        end
        
        %% 
        % See also <Slice.fcnHessian>.
        function h = hessReconfigure(vars, lambda, this, options) %#ok<INUSL>
            var_x = vars(1:options.num_varx);
            num_vars = length(vars);
            h = spalloc(num_vars, num_vars, options.num_varx^2);   % non-zero elements less than | options.num_varx^2|
            for p = 1:options.num_varx
                i = this.path_owner(p);
                h(p,1:options.num_varx) = this.weight *...
                    this.I_flow_path(i,:)/(1+(this.I_flow_path(i,:)*var_x))^2; %#ok<SPRIX>
            end
        end
    end
end
