classdef RequestEvent < handle
    properties (SetAccess = private)
        rand_state;
        arrive_times;
        depart_times;
        depart_id;
        arrive_type;
        density;
        event_set;
        
        total_arrival;
        avg_arrive_interval;
        avg_serve_interval;
        current_times = 0;
        current_arrive_pos;     % to arrive
        current_depart_pos;     % to depart
        accumulate_arrival;
        num_sojourns;
        eid;
        seed;
    end
    properties (Dependent = true)
        EventId;
    end
    
    methods
        %% Constructor
        % * |arrive_args|: specifies the parameters of the arrival process, including
        % _Interval_, _Number_;
        % * |event_set|: specifies the parameters of comming events, including
        % _Probability_, _ServiceInterval_;
        % * |seed|:
        function this = RequestEvent(event_set, arrive_args, seed)
            this.event_set = struct2table(event_set);
            this.event_set.Probability = ...
                this.event_set.Probability/sum(this.event_set.Probability);
            for i =1:height(this.event_set)
                if i>1
                    this.event_set.Probability(i) = ...
                        this.event_set.Probability(i-1) + this.event_set.Probability(i);
                end
            end
            
            this.total_arrival = arrive_args.Number;
            this.arrive_times = zeros(this.total_arrival,1);
            this.depart_times = zeros(this.total_arrival,1);
            this.depart_id = zeros(this.total_arrival,1);
            this.arrive_type = zeros(this.total_arrival,1);
            this.avg_arrive_interval = arrive_args.Interval;
            if nargin >= 3
                rng(seed);
                this.seed = seed;
            else
                rng('shuffle');
                this.seed = 1;
            end
            this.arrive_times(1) = exprnd(this.avg_arrive_interval);
            this.arrive_type(1) = this.nextType;
            this.depart_times(1) = this.arrive_times(1) + ...
                exprnd(this.event_set.ServiceInterval(this.arrive_type(1)));
            this.depart_id(1) = 1;
            this.current_arrive_pos = 1;
            this.current_depart_pos = 1;
            this.accumulate_arrival = 1;
            this.num_sojourns = 0;
            this.rand_state = rng;
            this.eid = 0;
        end
        
        function e = nextEvent(this)
            if this.current_depart_pos > this.accumulate_arrival
                error('No more event to be processed');
            end
            rng(this.rand_state);
            next_depart_time = this.depart_times(this.depart_id(this.current_depart_pos));
            % this.current_arrive_pos == this.accumulate_arrival || ...
            if this.arrive_times(this.accumulate_arrival) < next_depart_time
                while this.accumulate_arrival< this.total_arrival && ...
                        this.arrive_times(this.accumulate_arrival) < next_depart_time
                    this.accumulate_arrival = this.accumulate_arrival + 1;
                    % Next arrival's arrive time is the sum of the previous arrival time and the
                    % arrival interval.
                    this.arrive_times(this.accumulate_arrival) = ...
                        this.arrive_times(this.accumulate_arrival-1) + ...
                        exprnd(this.avg_arrive_interval);
                    t = this.nextType;
                    this.arrive_type(this.accumulate_arrival) = t;
                    this.depart_times(this.accumulate_arrival) = ...
                        this.arrive_times(this.accumulate_arrival) + ...
                        exprnd(this.event_set.ServiceInterval(t));
                    idx = this.accumulate_arrival;
                    t = idx;
                    while idx > 1 && ...
                            this.depart_times(t) < this.depart_times(this.depart_id(idx-1))
                        this.depart_id(idx) = this.depart_id(idx-1);
                        idx = idx - 1;
                    end
                    this.depart_id(idx) = t;
                    next_depart_time = this.depart_times(this.depart_id(this.current_depart_pos));
                end
            end
            if this.current_arrive_pos <= this.accumulate_arrival && ...
                    this.arrive_times(this.current_arrive_pos) < next_depart_time
                this.num_sojourns = this.num_sojourns + 1;
                this.eid = this.eid + 1;
                fprintf('Event %d: New arrival at %.4f, service interval is %.4f, %d in service.\n', ...
                    this.eid, ...
                    this.arrive_times(this.current_arrive_pos), ...
                    this.depart_times(this.current_arrive_pos)-...
                    this.arrive_times(this.current_arrive_pos),...
                    this.num_sojourns);
                
                e.Description = 'arrival';
                e.Type = this.arrive_type(this.current_arrive_pos);
                e.RandomSeed = this.seed + this.current_arrive_pos;
                e.Identifier = this.current_arrive_pos;
                e.Time = this.arrive_times(this.current_arrive_pos);
                this.current_arrive_pos = this.current_arrive_pos + 1;
            else
                this.num_sojourns = this.num_sojourns - 1;
                this.eid = this.eid + 1;
                fprintf('Event %d: New departure (%d) at %.4f, %d in service.\n', ...
                    this.eid, ...
                    this.depart_id(this.current_depart_pos), ...
                    this.depart_times(this.depart_id(this.current_depart_pos)), ...
                    this.num_sojourns);
                e.Description = 'departue';
                e.Id = this.depart_id(this.current_depart_pos);
                e.Type = this.arrive_type(this.depart_id(this.current_depart_pos));
                e.Time = this.depart_times(this.depart_id(this.current_depart_pos));
                this.current_depart_pos = this.current_depart_pos + 1;
            end
            this.rand_state = rng;
        end
        
        function T = countArriveType(this)
            arrive_pos = this.current_arrive_pos-1;
            tbl = tabulate(this.arrive_type(1:arrive_pos));
            if isempty(tbl)
                T = table([],[],[], 'VariableNames', {'Value', 'Count', 'Percent'});
            else
                T = array2table(tbl, 'VariableNames', {'Value', 'Count', 'Percent'});
            end
        end
        
        function T = countCurrentType(this)
            arrive_pos = this.current_arrive_pos-1;
            current_id = this.depart_id(this.current_depart_pos:arrive_pos);
            tbl = tabulate(this.arrive_type(current_id));
            if isempty(tbl)
                T = table([],[],[], 'VariableNames', {'Value', 'Count', 'Percent'});
            else
                T = array2table(tbl, 'VariableNames', {'Value', 'Count', 'Percent'});
            end
        end
        
        function reset(this)
            rng(this.seed);
            this.arrive_times(1) = exprnd(this.avg_arrive_interval);
            this.arrive_type(1) = this.nextType;
            this.depart_times(1) = this.arrive_times(1) + ...
                exprnd(this.event_set.ServiceInterval(this.arrive_type(1)));
            this.depart_id(1) = 1;
            this.current_arrive_pos = 1;
            this.current_depart_pos = 1;
            this.accumulate_arrival = 1;
            this.num_sojourns = 0;
            this.rand_state = rng;
            this.eid = 0;
        end
        
        function eid = get.EventId(this)
            eid = this.eid;
        end
    end
    
    methods (Access = private)
        function t = nextType(this)
            r = rand;
            for t = 1:height(this.event_set)
                if r<=this.event_set.Probability(t)
                    break;
                end
            end
        end
    end
end

