classdef Entity < matlab.mixin.Copyable & matlab.mixin.Heterogeneous
    
    properties (SetAccess = protected)
        EntityIdentifier;   % Global identifer of entities, might not be used.
        ArriveTime;
        ServiceTime;      
        Builder;            % Entity builder;
    end
    properties (Access = {?PhysicalNetwork, ?Entity})
        Description;        % Detailed information about the type; 
    end
    properties (Dependent)
        DepartTime;
        Type;               % Type of entity, see also <EntityType>;
    end
    
    methods
        function this = Entity(time_arrive, time_serve, src, info)
            persistent entity_id;
            if isempty(entity_id)
                entity_id = int64(1);
            else 
                if entity_id == intmax('uint64')
                    entity_id = int64(1);
                    warning('reset id');
                else
                    entity_id = entity_id + 1;
                end
            end
            this.ArriveTime = time_arrive;
            this.ServiceTime = time_serve;
            this.Builder = src;
            if nargin >= 4
                this.Description = info;
            end
            this.EntityIdentifier = entity_id;
        end
        
        function t = get.DepartTime(this)
            t = this.ArriveTime + this.ServiceTime;
        end
        
        function t = get.Type(this)
            t = this.Builder.Type;
        end
        
    end

    methods (Sealed)
        function tf = isPermanent(this)
            if isempty(this.ArriveTime) || isempty(this.ServiceTime)
                tf = true;
            else
                tf = false;
            end
        end
    end
end

