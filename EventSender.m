classdef (Abstract) EventSender < handle
    % <EventSender> and <EventReceiver> are interfaces for inter-class communication.
    
    properties (Access = protected)
        listeners;  % We need to remove some listeners, when the listeners have been destoryed.
        targets;
    end
        
    methods
        function this = EventSender()
            this.listeners = ListArray('event.listener');
            this.targets = ListArray('EventReceiver');
        end
        
        %%%
        % |event_name|: cell array represents a set of events.
        function AddListener(this, target, event_name, callback)
            if ischar(event_name)
                event_name = {event_name};
            end
            idx = this.targets(:)==target;
            for i = 1:length(event_name)
                if ~ismember(event_name{i},this.listeners(idx).EventName)
                    this.targets.Add(target);
                    this.listeners.Add(addlistener(this, event_name{i}, callback));
                    idx = [idx true]; %#ok<AGROW>
                end
            end
            %%
            % When notify listener (callback function), the EventData includes the Sender
            % and the Target of the event. In the callback function, the Target objects
            % can handle the event form the Sender.
        end
        
        % When slice is leaving, we should remove the slice form the listener's set
        function RemoveListener(this, target, event)
            idx = find(this.targets(:)==target);
            b_removed = false(length(this.listeners));
            for i = idx
                if strcmp(this.listeners(i).EventName,event)
                    b_removed(i) = true;
                end
            end
            this.listeners.Remove(b_removed);
            [~] = this.targets.Removed(b_removed);
        end
        
        function t = findTargets(this, ev_name)
            b_targets = false(length(this.listeners),1);
            for i = 1:length(this.listeners)
                if strcmp(this.listeners(i).EventName, ev_name)
                    b_targets(i) = true;
                end
            end
            t = this.targets(b_targets);
        end
    end
end

