%% RandomEventDispatcher
% generate random events.
% This class has enhanced the functionality of <RequestEvnet>, specifically including
% following improvements:
%
% # Introduce the concept of <Entity> and <EntityBuilder> for better organizing events.
% <Entity> represents the entity that produces events (Arrival/Departure), and
% <EntityBuilder> is used to generate these entities.
% # Introduce <ListArray> to store entity builders and event entities, and <PriorityQueue>
% to store events. Therefore, the maintenance of entities and events can be much
% simplified. Besides, the storage capacity for <ListArray> and <PriorityQueue> can be
% flexibly adjusted.
% # The entities are added to <ListArray> by arrival time, while the events are added to
% <PriorityQueue> by departure time. On the other hand, in <RequestEvent> two lists are
% maintained for arrival and departure events, and we need to sort the departure events.
% # 

classdef RandomEventDispatcher < handle
    
    properties (SetAccess = protected)
        rand_state;
        seed;
        
        entity_builder;
        avg_arrive_interval;    % update when entity source changes
        entity_probability;     % update when entity source changes
                
        entities;
        event_queue;
    end
    
    properties (Access = private)
        current_time = 0;
    end
    properties(Dependent)
        CurrentTime;
    end
    
    methods
        % * |event_set|: specifies the parameters of comming events, including
        % _ArrivalRate_, _ServiceInterval_;
        %    this = RandomEventDispatcher(entity_builder, seed, cur_time)
        function this = RandomEventDispatcher(entity_builder, varargin)
            if nargin == 0
                error('error: RandomEvent must be initialized with event entity.');
            end
            if nargin >= 1
                %%%
                % Calculate the accumulate probability
                if isempty(entity_builder)
                    entity_builder = ListArray('EntityBuilder');
                end
                if isa(entity_builder, 'ListArray') && ...
                        iscompatible(entity_builder.TypeName, 'EntityBuilder')
                    this.entity_builder = entity_builder.copy;
                elseif isa(entity_builder, 'EntityBuilder')
                    this.entity_builder = ListArray([], entity_builder);
                else
                    error('error: Invalid data for entity source.');
                end
                this.update_entity_builder;
            end
            
            if length(varargin) >= 1
                rng(varargin{1});
                this.seed = varargin{1};
                rng(this.seed);
            else
                warning('random number seed is not specified (set as shuffle).');
                rng('shuffle');
                scurr = rng;
                this.seed = scurr.Seed;
            end
            if length(varargin) >= 2
                this.current_time = varargin{2};
                %                 this.current_arrive_time = this.current_time;
            end
            
            this.entities = ListArray('Entity');
            priority = struct('field', 'Time', 'sorttype', 'ascend');
            this.event_queue = PriorityQueue('Event', priority);

            this.rand_state = rng;
        end
        
        function reset(this, seed)
            if nargin >= 2
                this.seed = seed;
            end
            
            this.entities.Clear;
            this.event_queue.Clear;
            
            rng(this.seed);
                        
            this.rand_state = rng;
            
            %             this.sojourn_type = zeros(this.NumberEventType,1);
            %             this.stat.total_arrival = 0;
            %             this.stat.arrive_type = zeros(this.NumberEventType,1);
            %             this.depart_id(1) = 1;
            %             this.current_arrive_pos = 1;
            %             this.current_depart_pos = 1;
            %             this.accumulate_arrival = 1;
        end
        
    end
    
    methods
        function t = get.CurrentTime(this)
            t = this.current_time;
        end
    end
    
    methods        
        function ev = nextEvent(this)
            rng(this.rand_state);
            if this.event_queue.Length == 0 % && this.entities.Length == 0
                if this.entity_builder.Length > 0
                    this.addnewentity;
                else
                    warning('No more events will be generated.');
                    ev = Event.empty;
                    return;
                end
            end
            
            if this.event_queue(1).Type == EventType.Depart
                this.addnewentity;
            end
            ev = this.event_queue.PopFront;
            this.current_time = ev.Time;
            %%%
            % when removing events, we need to check if we also need to remove the
            % associated Entity.
            if ev.Type == EventType.Depart
                ev.userdata = this.entities.Remove(ev.Entity);
                %                 for i = 1:this.event_queue.Length
                %                     if this.event_queue(i).Type == EventType.Depart
                %                         this.nextdeparttime = this.event_queue(i).Time;
                %                         break;
                %                     end
                %                 end
            else
                ev.userdata = ev.Entity;
            end
            this.rand_state = rng;
        end
        
        function AddEntityBuilder(this, src)
            this.entity_builder.Add(src);
            this.update_entity_builder;
        end
        
    end
    
    methods (Access=protected)
        function t = nextEntityId(this)
            r = rand;
            for t = 1:this.entity_builder.Length
                if r<=this.entity_probability(t)
                    break;
                end
            end
        end
        
        function update_entity_builder(this)
            arrive_rate = this.entity_builder(:).ArrivalRate;
            if ~isempty(arrive_rate)
                this.avg_arrive_interval = 1/sum(arrive_rate);
                this.entity_probability = cumsum(arrive_rate/sum(arrive_rate));
            end
        end
        
        %%%
        % When Entity Source is removed, it will no longer be used.
        % Before remove it, we may need to remove the associated events/entities.
        %
        % NOTE: this class might not need this function, while it can be used by child
        % class. See also <SliceFlowEventDispatcher.removeEntityBuilder>
        function removeEntityBuilder(this, argin, index)
            if isinteger(argin)
                eb_id = argin;
                for i = 1:this.entity_builder.Length
                    if this.entity_builder(i).Identifer == eb_id
                        %%%
                        % with no return value, the object of entity source will be
                        % deleted.
                        this.entity_builder.Remove(i);
                        break;
                    end
                end
            elseif isa(argin, 'EntityBuilder')
                if nargin >= 3
                    this.entity_builder.Remove(index);
                else
                    this.entity_builder.Remove(argin);
                end
            else
                error('error: invalid input arguments.');
            end
            this.update_entity_builder;
        end

        function e = addnewentity(this)
            arrive_time = this.current_time + exprnd(this.avg_arrive_interval);
            ei = this.nextEntityId;
            service_time = exprnd(this.entity_builder(ei).ServiceInterval);
            e = Entity(arrive_time, service_time, this.entity_builder(ei));
            this.entities.Add(e);
            %%%
            % PushBack can return the storage index in the queue.
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
        
        function t = randtime(this, mean_time)
            rng(this.rand_state);
            t = exprnd(mean_time);
            this.rand_state = rng;
        end
    end
    
    properties
        %% Deprecated
        % |arrive_times|, |depart_times|, and |arrive_type| can be accessed from entity
        % list or event queue.
        %         arrive_times;
        %         depart_times;
        %         arrive_type;     
        % <EventType> should be used with <EntityType> to identify a unique event type.
        % The information is stored in |entity_builder|.
        %         NumberEventType;
        % |depart_id| in <RequestEvent> record the index of departure events by the order
        % of depart time. This is not need in <RandomEventDispatcher>, since we have
        % sorted all events by depart time.
        %         depart_id;
        % |sojourn_type| can be obtained, by access entity list or event queue. This is
        % equivalent to |arrive_type|, since this class will remove the departed entities
        % and events.
        %         sojourn_type;
        %         NumberSojourns;
        %% TO BE UPDATED
        % These fields can be maintained to accelarate the access of entity list and event
        % queue.
        %         nextdeparttime = inf;
        %         current_arrive_time = 0;
        %         next_depart_pos = inf;      % Next depart event's location in the event queue
        %         last_arrive_pos;            % The most recent arrival event's location in the event queue.
        %         lastarrivetime;

        %% TODO
        % Log the statistics of entities and events from the beginning.
        %         stat;
        
    end
    
    methods (Access = private)
        %% TO BE UPDATED
        % generate events in the given interval.
        % NOTE: in this interval, the event source (entity builder) will not change. Since
        % this may not applied to the subclass (<SliceFlowEventDispatcher>, which can
        % dynamic add/remove entity builder, this function is accessed privately.
        %         function generateEvents(this, start_time, interval)
        %             current_time = start_time;
        %             while current_time < interval
        %                 arrive_time = current_time + exprnd(this.avg_arrive_interval);
        %                 type = this.nextType;
        %                 service_time = exprnd(this.event_set.ServiceInterval(type));
        %                 e = Entity(type, arrive_time, service_time);
        %                 this.entities.Add(e);
        %                 current_time = arrive_time;
        %             end
        %             depart_time = zeros(this.entities.Length,1);
        %             [~,idx] = sort(depart_time, 'ascend');
        %             i = 1;
        %             j = 1;
        %             while i<this.entities.Length && j < this.entities.Length
        %                 ta = this.entities{i}.ArriveTime ;
        %                 ts = this.entities{idx(j)}.DepartTime;
        %                 if ta <= ts
        %                     this.event_sequence.Add(Event(ta, 'SliceArrive', i));
        %                     i = i + 1;
        %                 else
        %                     this.event_sequence.Add(Event(ts, 'SliceDepart', idx(j)));
        %                     j = j + 1;
        %                 end
        %             end
        %             while i<this.entities.Length
        %                 ta = this.entities{i}.ArriveTime ;
        %                 this.event_sequence.Add(Event(ta, 'SliceArrive', i));
        %                 i = i + 1;
        %             end
        %             while i<this.entities.Length && j < this.entities.Length
        %                 ts = this.entities{idx(j)}.DepartTime;
        %                 this.event_sequence.Add(Event(ts, 'SliceDepart', idx(j)));
        %                 j = j + 1;
        %             end
        %         end
        %
        %         function t = get.nextdeparttime(this)
        %             t = this.event_queue(this.next_depart_pos).Time;
        %         end
        %
        % Last arrival time can be calculated from the Entity Queue, or from the Event
        % Queue, which you need to record the location of the last arrival event.
        %         function t = get.lastarrivetime(this)
        %             t = this.entities(end).ArriveTime;
        %             %%%
        %             % t = this.event_queue(this.last_arrive_pos).Time;
        %         end
        %
        %         function n = get.NumberEventType(this)
        %             n = height(this.event_set);
        %         end
        %         function n = get.NumberSojourns(this)
        %             n = sum(this.sojourn_type);
        %         end
        %         function ev = nextEvent(this)
        % %             this.next_depart_pos = this.next_depart_pos - 1;
        % %             this.last_arrive_pos = this.last_arrive_pos - 1;
        % %             if this.next_depart_pos == 0
        % %                 for i = 1:this.event_queue.Length
        % %                     if this.event_queue(i).Type == EventType.Depart
        % %                         this.next_depart_pos = i;
        % %                         break;
        % %                     end
        % %                 end
        % %             end
        %         end

    end
end

