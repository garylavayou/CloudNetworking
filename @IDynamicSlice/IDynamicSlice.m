classdef (Abstract) IDynamicSlice < EventSender & EventReceiver
	%DynamicSlice Event-driven to dynamic configure slice.
	properties (Constant)
		GLOBAL_OPTIONS = StaticProperties;
	end

	properties(Constant, Access = protected)
		NUM_MEAN_BETA = 15;
		ENABLE_DYNAMIC_NORMALIZER = true;
		GET_BETA_METHOD = 'Average';    % 'Average', 'ExponetialMovingAverage';
	end
	
	properties
		FlowArrivalRate;
		FlowServiceInterval;
		time;           % .{Current, LastDimensioning, Interval}
		bOnDepart = false;
	end
	
	properties
		portion_adhoc_flows_threshold = 0.1;
		b_adhoc_flows = zeros(1000,1);  % only record a window of flows.
		num_total_arrive = 0;
		net_changes = struct();
		% the variable |b_derive_vnf| decide if update VNF instance capacity.
		%    |b_derive_vnf=true|: derive VNF instance capacity from the optimization
		%    results of |Variables.z|;
		%    |b_derive_vnf=false|: apply |Variables.v| as VNF instance capacity.
		% For first time slice dimensioning,
		%   VNF capacity is the sum of VNF instance assigment (sum of z_npf);
		% For later reallocation (FSR2), the VNF capacity should be set as the optimized
		% variables (|this.Variables.v|, which might be larger than sum of |z_npf|, due to
		% reconfiguration cost).
		b_derive_vnf = true;
		b_dim = true;      % reset each time before reconfigurtion.
		old_net_state;
		vnf_reconfig_cost;
	end

	properties (Access = protected)
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
		old_variables;      % last one configuration, last time's VNF instance capcity is the |v| field;
		old_state;
		changed_index;
		diff_state;         % reset each time before reconfiguration.
		lower_bounds = struct([]);
		upper_bounds = struct([]);
		topts;              % used in optimization, avoid passing extra arguments.
		reject_index;
		max_flow_rate;
		init_gamma_k;
		init_q_k;
		% Final results of VNF instance capacity on each node, this is configured by
		% inter/intra-slicing, the values are stored in |Variables.v|; This property is a
		% wrapper.
		VNFCapacity;
	end
	
	events
		AddFlowSucceed;
		AddFlowFailed;
		RemoveFlowSucceed;
		RemoveFlowFailed;       % NOT used.
		RequestDimensioning;    % request slice dimensioning at once
		DeferDimensioning;      % defer slice dimensioning until the network makes dicision.
	end

	methods
		function this = IDynamicSlice(slice_data)
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
				{'TimeInterval', 'EventInterval', 'Adhoc', 'Trigger'}, ...
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
				case {ReconfigMethod.Dimconfig,ReconfigMethod.Fastconfig}
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
			
		function c = get.VNFCapacity(this)
			if ~isfield(this.Variables, 'v') || isempty(this.Variables.v)
				warning('VNF capacity not set, set to VNF load.');
				c = this.getVNFCapacity;
				this.Variables.v = c;
			else
				c = this.Variables.v;
			end
		end
		
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
						if this.options.Adhoc
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
		
		function tf = isDynamicFlow(this)
			if isempty(this.FlowArrivalRate) || isempty(this.FlowServiceInterval)
				tf = false;
			else
				tf = true;
			end
		end
		
		function b = isAdhocFlow(this)
			num_total = min(this.num_total_arrive,length(this.b_adhoc_flows));
			if num_total == 0 || nnz(this.b_adhoc_flows)/num_total < this.portion_adhoc_flows_threshold
				b = true;
			else
				b = false;
			end
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
				reconfig_cost_model = 'none';
			end
			
			link_price = this.Links.Price;
			link_load = this.Links.Capacity;
			node_price = this.ServiceNodes.Price;
			node_load = this.ServiceNodes.Capacity;
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
				n_path = nnz(this.I_dc_path)/this.NumberServiceNodes;
				b.v = b.z/n_path;
			end
		end
	end
	
	methods (Abstract)
		fidx = OnAddingFlow(this, flows);
		
		ef = OnRemovingFlow(this, fid);
		
		ds = diffstate(this, isfinal);
		stat = get_reconfig_stat(this, stat_names);
		s = get_state(this);
		set_state(this);
		[profit, cost] = handle_zero_flow(this, new_opts);
		[profit,cost] = fastReconfigure(this, action, options);
		identify_change(this, changed_path_index);
		update_reconfig_cost(this, action, bDim);
		[exitflag,fidx] = executeMethod(this, action);
	end
	
end

