classdef (Abstract) DynamicSlice < EventSender & EventReceiver & Slice
  %DynamicSlice Event-driven to dynamic configure slice.
  properties
    FlowArrivalRate;
    FlowServiceInterval;
    % time;           % .{Current, LastDimensioning, Interval}
    bOnDepart = false;
  end
  
  properties
    portion_adhoc_flows_threshold = 0.1;
    b_adhoc_flows = zeros(1000,1);  % only record a window of flows.
    num_total_arrive = 0;
    net_changes Dictionary;
		old_state Dictionary;
    old_net_state Dictionary;
		reject_index;
  end
  
  events
    AddFlowSucceed;
    AddFlowFailed;
    RemoveFlowSucceed;
    RemoveFlowFailed;       % NOT used.
    RequestDimensioning;    % request slice dimensioning at once
    DeferDimensioning;      % defer slice dimensioning until the network makes dicision.
	end
  
  methods (Abstract, Access=protected)
		release_resource_description(this);
	end

	%% Constructor
  methods
    function this = DynamicSlice(slice_data)
			if nargin == 0
				args = {};
			else
				args = {slice_data};
			end
			this@Slice(args{:});
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
			this.net_changes = Dictionary();
			this.old_state = Dictionary();
			this.old_net_state = Dictionary();
      this.options = structmerge(this.options, ...
        getstructfields(slice_data, {'Adhoc'}, 'default', {false}));
    end
	end

	%% Public Methods
	methods
		function s = save_state(this)
			this.old_state.erase();
			this.old_state.flow_table = this.FlowTable;
			this.old_state.link_load = this.Links.Load;
			this.old_state.link_capacity = this.Links.Capacity;
			this.old_state.node_load = this.ServiceNodes.Load;
			this.old_state.node_capacity = this.ServiceNodes.Capacity;
      this.old_state.path_owner = this.path_owner;
			
			sop = this.op.save_state();
      
      this.old_state = structmerge(this.old_state, sop);
			if nargout >= 1
				s = this.old_state;
			end			
		end
		
		function restore_state(this, s)
			if nargin <=1
				s = this.old_state;
				this.op.restore_state();
			end
			this.FlowTable = s.flow_table;
			% As the removing/adding flow operations do not change the original flow data, we
			% shallow copy the flow table.
			% As path is shallow copyed, its local identifier might have been changed when removing
			% flows. So in the restore stage, the ID is recaculated.
			pid = 1;
			this.path_owner = zeros(size(s.path_owner));
			for fid = 1:height(this.FlowTable)
				for p = 1:this.FlowTable{fid, 'Paths'}.Width
					this.FlowTable{fid, 'Paths'}{p}.local_id = pid;
					this.path_owner(pid) = fid;
					pid = pid + 1;
				end
			end
			assert( isequal(this.path_owner,s.path_owner) );
			this.Links.Load = s.link_load;
			this.Links.Capacity = s.link_capacity;
			this.ServiceNodes.Load = s.node_load;
			this.ServiceNodes.Capacity = s.node_capacity;
		end

		function finalize(this, prices)
			finalize@Slice(this, prices); % post processing to variables has been done
			
			op = this.Optimizer;
			% Override the capacity setting in the super class.
			if isfield(op.temp_vars,'v') && ~isempty(op.temp_vars.v)
				this.ServiceNodes.Capacity = full(sum(reshape(op.getVNFCapacity(), ...
					this.NumberServiceNodes,this.NumberVNFs),2));
			elseif this.NumberFlows > 0
				warning('Variables.v not initialized.');
				this.Optimizer.Variables.v = this.op.getVNFLoad();
			else
				warning('zero flow, Variables.v not initialized.');
			end
			
			if isfield(op.temp_vars,'c') && ~isempty(op.temp_vars.c)
				this.Links.Capacity = full(op.temp_vars.c);
			end
			if this.NumberFlows > 0
				%% Reconfiguration Cost
				% intra-slice reconfiguration: the flow reassignment cost and the VNF
				% instance reconfiguration cost is denpendent on the resource
				% consummption in the slice;
				% inter-slice reoncfiguration: the VNF node resource allocation cost is
				% dependent on the total load of the substrat network.
				% if strcmp(this.options.PricingPolicy, 'quadratic-price')
				[~, this.Links.ReconfigCost] = this.Optimizer.fcnLinkPricing(...
					this.Links.Price, this.Links.Capacity);
        go = IDynamicSliceOptimizer.GLOBAL_OPTIONS;
        if ~isfield(go, 'eta')
          warning('the global option ''eta'' should be initialized, but disappeared.');
          go.eta = this.Parent.options.UnitReconfigureCost;
        end
				eta = this.op.GLOBAL_OPTIONS.eta;
				this.Links.ReconfigCost = ...
					(eta/this.time.ConfigureInterval) * this.Links.ReconfigCost;
				% here the |node_price| is the price of all data centers.
				[~, this.ServiceNodes.ReconfigCost] = this.Optimizer.fcnNodePricing(...
					this.ServiceNodes.Price, this.ServiceNodes.Capacity);
				this.ServiceNodes.ReconfigCost = ...
					(eta/this.time.ConfigureInterval) * this.ServiceNodes.ReconfigCost;
			end
			% end
			this.op.init_reconfig_info();
		end
    
		%% fast slice reconfiguration when flow arriving and depaturing
		% * *flow*: flow table entries.
		function fidx = OnAddingFlow(this, flows)
			% temporarily allocated flow id.
			global DEBUG total_iter_num; %#ok<NUSED>
			num_exist_flows = height(this.FlowTable);
			num_new_flows = height(flows);
			this.save_state; % if fail to adding flow, recover flow table
			fidx = (1:num_new_flows) + num_exist_flows;
			if num_exist_flows == 0
				flows{:,'Identifier'} = 1:num_new_flows;
			else
				flows{:,'Identifier'} = (1:num_new_flows) + this.FlowTable{end,'Identifier'}; % temporary identifier
			end
			this.FlowTable = [this.FlowTable; flows];
			
			this.op.onAddingFlow(fidx);			% invoked after the flow is added.
			
			[ef, ~] = this.op.executeMethod('add');
			%% Failure handling
			% Zero-rate: now we leave it unhandled, the zero-rate flow will stay in
			% the slice. It may be allocated with resource at later stage. On the
			% other hand, we can remove the zero-rate flows from the slice (reject).
			if ef < 0
				fidx = [];
				this.restore_state();
			end
		end

		function ef = OnRemovingFlow(this, fid)
			global DEBUG; %#ok<NUSED>
			flow_index = find(this.FlowTable.Identifier == fid);
			if isempty(flow_index)
				warning('flow (identifier %d) not in slice (identifer %d).',...
					fid, this.Identifier);
			end
			this.save_state; % if fail to adding flow, recover slice data
			
			changes = this.op.onRemovingFlow(flow_index); % invoked before the flow is removed.
			this.FlowTable(flow_index,:) = [];
			this.updateFlowTable('remove', flow_index, changes);
			
			[ef,~] = this.op.executeMethod('remove');
			if ef < 0
				this.restore_state();
			end
		end
		
    function eventhandler(this, source, eventData) %#ok<INUSL>
      global DEBUG; %#ok<NUSED>
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
          if isempty(fidx)
            notify(this, 'AddFlowFailed');          % notify the network
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
          
      end
    end
    
    % <public>
    function tf = isDynamicFlow(this)
      if isempty(this.FlowArrivalRate) || isempty(this.FlowServiceInterval)
        tf = false;
      else
        tf = true;
      end
    end
    
    % <public>
    function b = isNextFlowAdhoc(this)
      num_total = min(this.num_total_arrive,length(this.b_adhoc_flows));
      if num_total == 0 || nnz(this.b_adhoc_flows)/num_total < this.portion_adhoc_flows_threshold
        b = true;
      else
        b = false;
      end
    end
    
    %%%
    % * *getCost*
    % Get the cost of creating a slice, overriding <Slice.getCost>. The cost
    % including resource consumption cost and reconfiguration cost.
    %   cost = getCost(this, pricing_policy, reconfig_cost_model)
    % |reconfig_cost_model|: reconfigure cost model, including: |'const'|,
    % |'linear'|, and |'none'|. If |'none'| is adopted, then the reconfiguration cost
    % is not counted.
		function cost = getCost(this, reconfig_cost_model)
			if nargin <= 1 || isempty(reconfig_cost_model)
				reconfig_cost_model = 'none';
			end
			
			link_price = this.Links.Price;
			link_load = this.Links.Capacity;
			node_price = this.ServiceNodes.Price;
			node_load = this.ServiceNodes.Capacity;
			switch this.options.PricingPolicy
				case {'quadratic-price', 'quadratic'}
					link_payment = this.Optimizer.fcnLinkPricing(link_price, link_load);
					node_payment = this.Optimizer.fcnNodePricing(node_price, node_load);
					cost = link_payment + node_payment;
				case 'linear'
					cost = dot(link_price, link_load) + dot(node_price, node_load);
				otherwise
					error('%s: invalid pricing policy', calledby);
			end
			if ~strcmpi(reconfig_cost_model, 'none')
				cost = cost + this.op.get_reconfig_cost(reconfig_cost_model, true);
			end
		end
		
		function stat = get_reconfig_stat(this, stat_names)
			if nargin == 1
				stat_names = {'All'};
			elseif ischar(stat_names)
				stat_names = {stat_names};
			end
			
			stat = table;
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
					stat.Cost = this.op.get_reconfig_cost('const');
				end
				if contains(stat_names{i},{'All', 'LinearCost'},'IgnoreCase',true)
					stat.LinearCost = this.op.get_reconfig_cost('linear', true);
				end
				%% Number of Variables
				if contains(stat_names{i},{'All', 'NumVariables'},'IgnoreCase',true)
					stat.NumVariables = this.op.get_reconfig_stat('NumVariables');
				end
				if contains(stat_names{i},{'All', 'ReVariables'},'IgnoreCase',true)
					stat.ReVariables = this.op.get_reconfig_stat('ReVariables');
				end
				%% Number of Flows
				if contains(stat_names{i},{'All', 'Flows'},'IgnoreCase',true)
					% For convenience of comparison, we store the number of flows including the
					% removed one.
					stat.Flows = max(this.NumberFlows, height(this.old_state.flow_table));
				end
				if contains(stat_names{i},{'All', 'ReFlows'},'IgnoreCase',true)
					stat.ReFlows = this.op.get_reconfig_stat('ReFlows');
				end
				if contains(stat_names{i},{'All', 'ResourceCost'},'IgnoreCase',true)
					stat.ResourceCost = this.getCost('none');
				end
				if contains(stat_names{i},{'All', 'FairIndex'},'IgnoreCase',true)
					stat.FairIndex = (sum(this.FlowTable.Rate))^2/...
						(this.NumberFlows*sum(this.FlowTable.Rate.^2));
				end
				if contains(stat_names{i},{'All', 'Interval'},'IgnoreCase',true)
					if this.op.b_dim
						stat.Interval = this.time.DimensionIntervalModified;
					else
						stat.Interval = this.time.ConfigureInterval;
					end
				end
			end
		end
		
	end
  
	methods (Access = {?DynamicSlice, ?IDynamicSliceOptimizer})
		function updateFlowTable(this, action, varargin)
			switch action
				case 'remove'
					%% Update flow table and local path id.
					% After removing flows, re-allocate local identifier to subsequent flows (flow
					% ID) and paths (path ID).
					% old information of flow table required: number of flows,
					flow_index = varargin{1};
					changes = varargin{2};
					pid = changes.pid;
					fidt = flow_index(1):(changes.Nf-length(flow_index));
					this.path_owner(changes.pid) = [];
					for fid = fidt
						path_list = this.FlowTable{fid, 'Paths'};
						for j = 1:path_list.Width
							path_list{j}.local_id = pid;
							this.path_owner(pid) = fid;
							pid = pid + 1;
						end
					end
				case 'add'
					flow_index = varargin{1};
					changes = varargin{2};
					pid = changes.pid;
					for i = 1:length(flow_index)
						fid = flow_index(i);
						for p = 1:this.FlowTable{fid, 'Paths'}.Width
							pid = pid + 1;
							this.path_owner(pid) = fid;
							this.FlowTable{fid, 'Paths'}{p}.local_id = pid;
						end
					end
			end
		end
		
	end
	
	methods (Access = protected)
		function newobj = copyElement(this)
			if this.isShallowCopyable
				newobj = copyElement@Slice(this);
				newobj.isShallowCopyable = false;
				newobj = copyElement@EventSender(newobj);
				newobj.isShallowCopyable = true;
			else
				newobj = this;
			end
			%% Reset the listener of the new instance
			% We should reconfigure the listeners by using AddListeners
			% outside.
			% see <EventSender>, <RepeatSliceReconfiguration>.
			newobj.ClearListener();
		end
	end
end

