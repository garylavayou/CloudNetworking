%% Enity Builder
% *EntityBuilder* stores necessary information, which is used to generate information for
% creating entities. 
%
% _copyElement_: no handle class member, no need to override the method.
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
            if nargin >= 3
                this.Options = options;
            end
            global builder_id;       % global entity ID;
            if isempty(builder_id)
                builder_id = int64(0);
            end
            if builder_id == intmax('uint64')
                builder_id = int64(1);
                warning('reset id');
            else
                builder_id = builder_id + 1;
            end
            this.Identifier = builder_id;
        end
    end
    
    methods (Abstract)
        entity = Build(this, time_arrive, time_serve, varargin);
    end
end

