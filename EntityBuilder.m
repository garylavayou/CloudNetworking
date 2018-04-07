%% Enity Builder
% *EntityBuilder* stores necessary information, which is used to generate information for
% creating entities. 
%
% _copyElement_: no handle class member, no need to override the method.
classdef EntityBuilder < matlab.mixin.Copyable & matlab.mixin.Heterogeneous
    
    properties
        ArrivalRate;
        ServiceInterval;
        Identifier uint64;  % Global identifier of entity builder
        Options;            % Provide informations to build an object corresponding to the entity.
    end
    
    properties (Constant, Access = private)
        sn = SerialNumber;  
        % NOTE: Each type of class need to define their own serial number;
        % If a serieral number is shared among multiple classes, we can define a common
        % super class for them.
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
            this.Identifier = EntityBuilder.sn.next();
        end
    end
    
    methods (Abstract)
        entity = Build(this, time_arrive, time_serve, varargin);
    end
    
    methods(Static)
        function n = getGlobalBuilderId()
            n = EntityBuilder.sn.ID;
        end
        function setGlobalBuilderId(value)
            h = EntityBuilder.sn;
            h.ID = value;
        end
    end
end

