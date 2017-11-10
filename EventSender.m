classdef (Abstract) EventSender < matlab.mixin.Copyable
    % <EventSender> and <EventReceiver> are interfaces for inter-class communication.
    
    properties (Access = protected)
        listeners;  % We need to remove some listeners, when the listeners have been destoryed. [Deprecated].
        targets;
    end
    
    properties 
        Results;    % Used by <EventReceiver>, return the event handling results.
                    % The data structure is defined by <EventSender>, and should at least
                    % include a field named as 'Value'. <EventReceiver> should return
                    % infromation requested by <EventSender>.
    end
        
    methods
        function this = EventSender()
            this.listeners = ListArray('event.listener');
            this.targets = ListArray('EventReceiver');
        end
        
        function delete(this)
            delete(this.listeners);
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
            if nargin >= 3
                for i = idx
                    if strcmp(this.listeners(i).EventName,event)
                        b_removed(i) = true;
                    end
                end
            else
                b_removed(idx) = true;
            end
            this.listeners.Remove(b_removed);       % delete the listener objects;
            [~] = this.targets.Removed(b_removed);  % remove but not delete the targets;
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
    
    methods (Access = protected)
        %% Deep Copy
        function this = copyElement(es)
            % Make a shallow copy of all properties
            this = copyElement@matlab.mixin.Copyable(es);
            % Make a deep copy of the DeepCp object
            this.listeners = ListArray('event.listener');
            %% ISSUE
            % After copy, the new instance still send message to the old targets, since
            % we do not know whether the targets will be changed or not, inside the copy
            % method. 
            warning('Targets are not deep copyed.');
            for i=1:this.targets.Length
                this.listeners.Add(addlistener(this, ...
                    es.listeners(i).EventName, es.listeners(i).callback));
            end
        end
    end
end

