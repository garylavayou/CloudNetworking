%% Slice Entity Builder
% In addition to functionality of <EntityBuilder>, *SliceEntityBuilder* can generate the
% global identifiers for slice entities.
classdef SliceEntityBuilder < EntityBuilder
        
    properties (Dependent)
        Type;
    end
    properties(Access = private)
        seed;         % used to generate random seed for enities
    end
    
    methods
        %%
        % See also <EntityBuilder>.
        % To create permanent entities, pass empty value of  ([]) to the constructor.
        %   this = SliceEntityBuilder(slice_options)
        %   this = SliceEntityBuilder(arrival_rate, service_interval, slice_options)
        % In the first prototype, |slice_options| include the |arrival_rate| and
        % |serivce_interval|.
        function this = SliceEntityBuilder(varargin)
            if isstruct(varargin{1})
                slice_opt = varargin{1};
                arrive_rate = slice_opt.ArrivalRate;
                service_interval = slice_opt.ServiceInterval;
%                 slice_opt = rmfield(slice_opt, {'ArrivalRate', 'ServiceInterval'});
            elseif nargin >= 3
                arrive_rate = varargin{1};
                service_interval = varargin{2};
                slice_opt = varargin{3};
            else
                error('error: arguments not enough.');
            end
            if isfield(slice_opt, 'RandomSeed')
                seed = slice_opt.RandomSeed;
                slice_opt = rmfield(slice_opt, 'RandomSeed');
            else
                warning('random seed not specified for slice entity builder.');
                seed = 0;
            end
            this@EntityBuilder(arrive_rate, service_interval, slice_opt);
            this.seed = seed + this.Identifier;
        end
        
        function t = get.Type(this) %#ok<MANU>
            t = EntityType.Slice;
        end
    end
    
    methods
        function entity = Build(this, time_arrive, time_serve)
            global slice_id;    % Global slice ID;
            if isempty(slice_id)
                slice_id = int64(0);
            end
            if slice_id == intmax('uint64')
                slice_id = int64(1);
                warning('reset slice id.');
            else
                slice_id = slice_id + 1;
            end

            % build the entity
            entity = SliceEntity(time_arrive, time_serve, this, this.seed);
            this.seed = this.seed + 1;
            entity.SliceIdentifier = slice_id;
        end        
    end
    
end

