classdef (Abstract) DynamicSlice < EventSender & EventReceiver & Slice
  %DynamicSlice Event-driven to dynamic configure slice.
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
		old_state;
    old_net_state;
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
  
  methods (Abstract)
    fidx = OnAddingFlow(this, flows);
    
    ef = OnRemovingFlow(this, fid);
    
    ds = diffstate(this, isfinal);
    stat = get_reconfig_stat(this, stat_names);
		release_resource_description(this);
  end
  
  methods (Abstract, Access=protected)
    s = save_state(this);
    restore_state(this);
  end
  
  methods
    function this = DynamicSlice(slice_data)
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
      this.options = structmerge(this.options, ...
        getstructfields(slice_data, {'Adhoc'}, 'default', {'false'}));
    end
    
    function finalize(this, prices)
      finalize@Slice(this, prices);

      if this.NumberFlows == 0
        this.op.clear();
        this.ServiceNodes{:,'Load'} = 0;
        this.Links{:,'Load'} = 0;
        this.ServiceNodes{:,'Capacity'} = 0;
        this.Links{:,'Capacity'} = 0;
      else
        %% Reconfiguration Cost
        % intra-slice reconfiguration: the flow reassignment cost and the VNF
        % instance reconfiguration cost is denpendent on the resource
        % consummption in the slice;
        % inter-slice reoncfiguration: the VNF node resource allocation cost is
        % dependent on the total load of the substrat network.
        % if strcmp(this.options.PricingPolicy, 'quadratic-price')
        [~, this.Links.ReconfigCost] = this.op.fcnLinkPricing(...
          this.Links.Price, this.Links.Capacity);
        eta = this.op.GLOBAL_OPTIONS.get('eta');
        this.Links.ReconfigCost = ...
          (eta/this.time.ConfigureInterval) * this.Links.ReconfigCost;
        % here the |node_price| is the price of all data centers.
        [~, this.ServiceNodes.ReconfigCost] = this.op.fcnNodePricing(...
          this.ServiceNodes.Price, this.ServiceNodes.Capacity);
        this.ServiceNodes.ReconfigCost = ...
          (eta/this.time.ConfigureInterval) * this.ServiceNodes.ReconfigCost;
      end
      % end
      this.op.init_reconfig_info();
    end
    
  end
  
  methods
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
    function cost = getCost(this, pricing_policy, reconfig_cost_model)
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
          link_payment = this.op.fcnLinkPricing(link_price, link_load);
          node_payment = this.op.fcnNodePricing(node_price, node_load);
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
    
  end
  
end

