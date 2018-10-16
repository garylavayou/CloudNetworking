classdef (Abstract) EventSender < matlab.mixin.Copyable
    % <EventSender> and <EventReceiver> are interfaces for inter-class communication.
    
    properties (Access = protected)
        listeners;  % We need to remove some listeners, when the listeners have been destoryed. [Deprecated].
        targets;
				bCopyReady logical = true; 
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
				
				function ClearListener(this)
					delete(this.listeners);
					this.listeners = ListArray('event.listener');
					this.targets.Clear();
					this.bCopyReady = true;
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
			%% Copy the Listeners and Targets
			% After the default shallow copy, the new instance has the
			% same event targets as the current instance. 
			% However, the stored listeners in the 'newobj.listeners' instance is
			% the listeners of 'this'. Thus the listeners for the new instance
			% should be regenerated with the same set of targets as 'this'.
			%
			% After the copy is finished, the listeners might be changed outside.
			% Then we need to first <delete> the existing listeners and then add
			% new listeners.
        function newobj = copyElement(this)
            % Make a shallow copy of all properties
            newobj = copyElement@matlab.mixin.Copyable(this);
            % Make a deep copy of the DeepCp object
            newobj.listeners = ListArray('event.listener');
						cprintf('Comments', 'info:[%s] copy listeners to the new object.\n', calledby);
						for i=1:newobj.targets.Length
							newobj.listeners.Add(addlistener(newobj, ...
								this.listeners(i).EventName, this.listeners(i).callback));
						end
						
						this.bCopyReady = false;
        end
    end
end

