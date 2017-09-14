%% Enity Builder
% *EntityBuilder* stores necessary information, which is used to generate information for
% creating entities. 
classdef EntityBuilder < matlab.mixin.Copyable & matlab.mixin.Heterogeneous
    
    properties
        ArrivalRate;
        ServiceInterval;
        Identifier;         % Global identifier of entity builder
        Options;            % Provide informations to build an object corresponding to the entity.
    end
    
    properties(Abstract, Dependent)
        Type;               % EntityType
    end
    
    methods
        %% Constructor
        % |arrival_rate| and |service_interval| are used by <EventDispatcher> to compute
        % the arriving events. The <EntityBuilder> itself only stores these parameters.
        %
        % NOTE: this class is abstract, no need to check the input arguments (perform the
        % check in concrete subclasses). 
        function this = EntityBuilder(arrival_rate, service_interval, options)
            this.ArrivalRate = arrival_rate;
            this.ServiceInterval = service_interval;
            this.Options = options;
            persistent builder_id;       % global entity ID;
            if isempty(builder_id)
                builder_id = int64(0);
            else 
                if builder_id == intmax('uint64')
                    builder_id = int64(1);
                    warning('reset id');
                else
                    builder_id = builder_id + 1;
                end
            end
            this.Identifier = builder_id;
        end
    end
    
    methods (Sealed)  % This method should have been abstract
        function entity = Build(this, time_arrive, time_serve, varargin)
            entity = this.buildentity(time_arrive, time_serve, varargin{:});
        end
    end
    
    methods (Access=protected, Abstract)
        entity = buildentity(this, time_arrive, time_serve, varargin);
    end
end

