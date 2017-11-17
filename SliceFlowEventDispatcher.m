classdef SliceFlowEventDispatcher < RandomEventDispatcher & EventSender & EventReceiver
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here   
    properties (Constant)
        flow_depart_option = 'NaturalDepart';    % MandatoryDepart
    end
    %% Events and event handler
    events
        SliceArrive;
        SliceDepart;
        FlowArrive;
        FlowDepart;
        %PermanetArrival;        % permanent slice arrival
    end
    
    methods
        %%%
        % |source|: potential source might be <DynamicNetwork>.
        % |evenData|: <DispatchEventData> object with prperties including |event|,
        % |entity| and |userdata|.
        function eventhandler(this, source, eventData) %#ok<INUSL>
            global DEBUG; %#ok<NUSED>
            % target = eventData.targets; % where targets should be the <EventDispatcher>
            switch eventData.EventName
                case 'AddSliceSucceed'
                    % To add existing flows and the flow entity builder.
                    % |userdata| is the successfully added slice.
                    sl = eventData.userdata;
                    if sl.isDynamicFlow
                        flow_entity_builder = ...
                            FlowEntityBuilder(sl.FlowArrivalRate, sl.FlowServiceInterval, eventData.entity);
                        for i = 1:sl.NumberFlows
                            fe = flow_entity_builder.Build(this.CurrentTime, ...
                                this.randtime(flow_entity_builder.ServiceInterval),...
                                sl.FlowTable{i, 'Identifier'});
                            %%%
                            this.AddEntity(fe);
                        end
                        this.AddEntityBuilder(flow_entity_builder);
                    end
                case 'AddSliceFailed'
                    % Remove the failed slice entity;
                    warning(eventData.EventName);
                    % When the entity is removed, the related slice event is invalid,
                    % and removed when it is dequeuing.
                    this.entities.Remove(eventData.entity);
                case 'RemoveSliceSucceed'
                    if isa(eventData.slice, 'Slice')
                        this.RemoveListener(eventData.slice);
                    else
                        error('error: type error.');
                    end                    
                case 'RemoveSliceFailed'
                    pause;
                case 'AddFlowSucceed'
                    flow_entity = eventData.entity;
                    flow_entity.GlobalIdentifier = eventData.userdata;
                    % We can find the parent of flow(entity) by |flow_entity.Parent.Identifier|
                case 'AddFlowFailed'
                    pause;
                case 'RemoveFlowSucceed'
                    % The flow entity has been pop-out from the entity list.
                    % We should leave it alone, since it will be used after the event
                    % handler. See alse <nextEvent>.
                    % It will be cleaned automatically, when there is no reference to it.
                    %     eventData.userdata.delete;
                case 'RemoveFlowFailed'
                    pause;
                otherwise
                    error('error: can not handle %s', eventData.EventName);
            end
        end
        
        function sendoutevent(this, event)
%             rng(this.rand_state);
%             this.rand_state = rng;
            data = DispatchEventData;
            %             data.targets = this.findTargets(event.Name);
            % data.entity = entity;  % the event object stores the entity handle.
            data.event = event;
            notify(this, event.Name, data);
        end
    end
    
    %% Constructor
    methods
        %%%
        % *Parameters*
        %
        % |entity_source|:
        % |seed|: seed for the internal random number generator.
        % |start_time|: start time of the simulation.
        % TODO: pass the argument to initialize 'flow_depart_option' .
        function this = SliceFlowEventDispatcher(entity_source, varargin)
            this@RandomEventDispatcher(entity_source, varargin{:});
            this@EventSender();
        end
    end
    
    %% Methods
    methods
        %%%
        % # Add flow entities for initial flows in the slices.
        % # Add permanent slice entities
        function AddEntity(this, et)
            et = this.entities.Add(et);
            if et.isPermanent
                % Permanent slice or flow
                this.event_queue.PushBack(Event(this.CurrentTime,EventType.Arrive,et));
            elseif this.CurrentTime >= et.ArriveTime
                % the flow events with slice arrival, these flow arrival events have been
                % processed together with the slice arrival events.
                this.event_queue.PushBack(Event(et.DepartTime,EventType.Depart,et));
            else
                this.event_queue.PushBack(Event(et.ArriveTime,EventType.Arrive,et));
                this.event_queue.PushBack(Event(et.DepartTime,EventType.Depart,et));
            end
        end
        
        function ev = nextEvent(this)
            %             rng(this.rand_state);
            while true
                %%%
                % #Option 1:
                % When the event is dequeued form the queue, we can choose to discard the
                % event that is no longer valid (has no valid entity/source associated.)
                %this.rand_state = rng;
                ev = nextEvent@RandomEventDispatcher(this);
                assert(issorted(this.event_queue.Time), 'error: EventQueue.');
                if isvalidEvent(this, ev)
                    break;
                end
            end
            %             this.rand_state = rng;
            et = ev.userdata;       % equals to ev.Entity
            if isempty(et)
                % When |event=Depart|, the entity will be returned as userdata, if this is
                % not the truth, there is an error. Besides, the entity handle also be
                % stored in the event object.
                error('error: Entity unknown.');
            else
                %% Slice Arrival
                % We need to inform the network that a slice is arriving, and let the
                % network decide if this slice can be successfully admitted.
                global event_num; %#ok<TLEV>
                event_num = event_num + 1;
                if this.targets.Length > 0
                    this.sendoutevent(ev);
                end
                % After that, the network will notify the EventDispatcher the results.
                % Then, the EventDispatcher will dicide if more events will be added
                % to the queue.
                %%%
                % #Option 2:
                % we can choose to remove exsiting Flow Entity/Event associated with
                % the departing Slice from the entity list/event queue. After that, we
                % can remove flow entity source from the entity source list.
                %
                % However, remove elements is not supported by event queue, it is not
                % recommended to use this method. If we use this option, we should
                % implement the event queue as a List/Vector.
                switch ev.Name
                    case 'SliceDepart'
                        % A slice entity has been removed (a slice departed).
                        % Then we should remove the flow entity builder derived from the
                        % slice.
                        this.removeFlowEntityBuilder(et);
                    otherwise
                end
            end
        end
    end
    
    methods (Access = protected)
         function this = copyElement(ed)
            this = copyElement@RandomEventDispatcher(ed);
            %%
            % The copyed version may not have the same targets as the copy source. We can
            % mannually update the target/listener list using AddListener/RemoveListener.
            %{
              temp = copyElement@EventSender(ed);
              this.targets = temp.targets;
              this.listeners = temp.listeners;
            %}
            % To make the new object not influence the original one, we detach the link of targets
            % and listeners.
            this.listeners = ListArray('event.listener');
            this.targets = ListArray('EventReceiver');
         end
        function removeFlowEntityBuilder(this, parent)
            for i = 1:this.entity_builder.Length
                eb = this.entity_builder(i);
                if isa(eb,'FlowEntityBuilder') && eb.Parent == parent
                    switch this.flow_depart_option
                        case 'naturaldepart'
                        case 'mandatorydepart'
                            %%%
                            % Remove residual flow entities from the list.
                            % Those removed entities will not be used any more, it can be
                            % deleted.
                            % See also <RandomEventDispatcher.removeEntitySource>
                            b_remove = false(this.entities.Length,1);
                            for ei = 1:this.entities.Length
                                if this.entities(ei).Builder == entity_builder
                                    b_remove(ei) = true;
                                end
                            end
                            this.entities.Remove(b_remove);
                            %%%
                            % Since the last event of one entity is a Depart event, when
                            % we find the depart event, we can stop the remove procedure.
                            % Issue: cannot remove all events from queue. Due to this
                            % issue we do not clear the residual events.
                    end
                    %%%
                    % Call super class method: remove the flow entity builder
                    this.removeEntityBuilder(entity_builder, i);
                    break;
                end
            end
        end
        %% Deprecated
        %         function removeEntityBuilder(this, entity_builder)
        function removeEntity(this, et)
            this.entities.Remove(et);
        end
        
        function tf = isvalidEvent(~, ev)
            %%% If the entity has been removed, the event is invalid.
            if isvalid(ev.Entity)
                tf = true;
            else
                tf = false;
            end
        end
        
        function e = addnewentity(this)
            arrive_time = this.CurrentTime + exprnd(this.avg_arrive_interval);
            ei = this.nextEntityId;
            service_time = exprnd(this.entity_builder(ei).ServiceInterval);
            e = this.entity_builder(ei).Build(arrive_time, service_time);
            this.entities.Add(e);
            %             this.last_arrive_pos = ...
            this.event_queue.PushBack(...
                Event(arrive_time, EventType.Arrive, this.entities(end)), 'front');
            
            %             depart_pos = ...
            this.event_queue.PushBack(...
                Event(e.DepartTime, EventType.Depart, this.entities(end)));
            %             if depart_pos < this.next_depart_pos
            %                 this.next_depart_pos = depart_pos;
            %             end
%             if e.DepartTime < this.nextdeparttime
%                 this.nextdeparttime = e.DepartTime;
%             end
%             this.current_arrive_time = arrive_time;
        end
    end
    
end

