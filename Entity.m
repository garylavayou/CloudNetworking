classdef Entity < matlab.mixin.Copyable & matlab.mixin.Heterogeneous
    
    properties (SetAccess = protected)
        EntityIdentifier;   % Global identifer of entities, might not be used.
        ArriveTime;
        ServiceTime;
    end
    properties (Access = {?PhysicalNetwork, ?Entity})
        Description;        % Detailed information about the type;
    end
    properties (SetAccess = {?Entity, ?RandomEventDispatcher})
        Builder;            % Entity builder;
    end
    properties (Dependent)
        DepartTime;
        Type;               % Type of entity, see also <EntityType>;
    end
    
    methods
        function this = Entity(time_arrive, time_serve, src, info)
            % replace <persistent>, thus we can backup and recover the identifier out of the class.
            global entity_id;
            if isempty(entity_id)
                entity_id = int64(0);
            end
            if entity_id == intmax('uint64')
                entity_id = int64(1);
                warning('reset id');
            else
                entity_id = entity_id + 1;
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
    
    methods (Access = protected)
        function this = copyElement(et)
            this = copyElement@matlab.mixin.Copyable(et);
            %% Deep Copy Issue
            % *Builder* is an exterior link. When performing copy, we should not make a copy of this
            % object. Instead, the link should be updated by the caller of the _copy_ function. To
            % secure the original data, we detach the link in the new copy from the original data.
            if ~isempty(et.Builder)
                this.Builder = et.Builder.empty();
            end
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

