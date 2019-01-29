%% Request Envent Generator
%% Arrival rate
% *Compund Poisson Process*:
% We assume each kind of events follows a Poissson process (with arrival rate
% $\lambda_s$), and they are independent. Therefore, arrival of requests can be seen as a
% compound Poission process, with arrival rate $\lambda = \sum_s{\lambda_s}$. On the other
% hand, if the event of a Poisson process (with arrival rate $\lambda$) can be classified
% as sub-types with the probability $p_i, \sum_i{p_i}=1$ of each type of events, the
% process can be decomposed as a group of Poisssion process with respect to each type of
% events. In this way, the arrival rate of each type of events is given by $\lambda_i =
% p_i \lambda$. 
%
% *Emulation of Poission Process*:In the simulation, we emulate the Poission process by
% genereating the interval between two arrivals, which follows an exponential
% distribution (with expectation of |T|, and thus the arrival rate is
% $\lambda=\frac{1}{T}$). Therefore we use exponential random number generator to produce
% the random interval. Once a event arrives, we determine its type according to the
% probability $p_i$ by using a uniform random generator. Thus the arrival rate of type |i|
% event is $\lambda_i = p_i \lambda = \frac{T}{p_i}$.    
%
%% Service interval 
% The average interval influences the average number of customers in the system (_the
% Little's Therom_). Given the average arrival rate $\lambda$ and average service interval
% $\frac{1}{\mu}$, we have the average number of customers $L = \frac{\lambda}{\mu}$.
% Specially, we set the service interval following exponential distribution, while any
% distribution with average service interval $\frac{1}{\mu}$ is acceptable. 
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
        % * |event_set|: specifies the parameters of comming events, including
        % _Probability_, _ServiceInterval_;
        % * |arrive_args|: specifies the parameters of the arrival process, including
        % _Interval_(average arrival interval), _Number_(number of arrivals);
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
                warning('random number seed is not specified (set as shuffle).');
                rng('shuffle');
                scurr = rng;
                this.seed = scurr.Seed;
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
        %%%
        % Public methods and constructor must be surrounded with \rng(this.rand_state)|
        % and |this.rand_state = rng|, so that the random generator's state is continuous
        % between two calls to the methods. On the other hand, the private methods do not
        % need to keep the states, since they will only be called by the public methods
        % and constructor.
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
                fprintf('(%s) Event %d: New arrival at %.4f, service interval is %.4f, %d in service.\n', ...
                    datestr(now),...
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
                fprintf('(%s) Event %d: New departure (%d) at %.4f, %d in service.\n', ...
                    datestr(now),...
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
          % current_arrive_pos is the position for next arrival.
          arrive_pos = this.current_arrive_pos-1;
          [types, counts, percents] = count(this.arrive_type(1:arrive_pos));
          T = table(types, counts, percents, 'VariableNames', {'Value', 'Count', 'Percent'});
        end
        
        function T = countCurrentType(this)
          arrive_ids = 1:(this.current_arrive_pos-1);
          depart_ids = this.depart_id(1:this.current_depart_pos - 1);
          % SETDIFF returns elements in |arrive_ids| but not in |depart_ids|.
          current_ids = setdiff(arrive_ids, depart_ids);
          [types, counts, percents] = count(this.arrive_type(current_ids));
          T = table(types, counts, percents, ...
            'VariableNames', {'Value', 'Count', 'Percent'});
        end
        
        function reset(this)
            rng(this.seed);
            this.arrive_times = zeros(this.total_arrival,1);
            this.arrive_times(1) = exprnd(this.avg_arrive_interval);
            this.arrive_type = zeros(this.total_arrival,1);
            this.arrive_type(1) = this.nextType;
            this.depart_times = zeros(this.total_arrival,1);
            this.depart_times(1) = this.arrive_times(1) + ...
                exprnd(this.event_set.ServiceInterval(this.arrive_type(1)));
            this.depart_id = zeros(this.total_arrival,1);
            this.depart_id(1) = 1;
            this.current_arrive_pos = 1;
            this.current_depart_pos = 1;
            this.accumulate_arrival = 1;
            this.num_sojourns = 0;
            this.eid = 0;
            this.rand_state = rng;
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

