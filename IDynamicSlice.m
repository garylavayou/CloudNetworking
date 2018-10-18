classdef (Abstract) IDynamicSlice < EventSender & EventReceiver
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
		net_changes;
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
		
	end
	
	methods (Abstract)
		fidx = OnAddingFlow(this, flows);
		
		ef = OnRemovingFlow(this, fid);
	end
end

