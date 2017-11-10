classdef Event < matlab.mixin.Copyable
    
    properties (SetAccess = protected)
        Type;           % EventType
        Time = 0;
        Identifer;      
    end
    properties (Access = private)
        Description;
    end
    properties (SetAccess = {?Event, ?RandomEventDispatcher})
        Entity;         % Entity handle can be checked with isvalid().
    end
    properties (Dependent = true)
        EntityType;   
        Name;
    end
    properties  % public access properties
        userdata;
    end
    
    methods
        function this = Event(time, type, entity, desc)
            if nargin == 0
                return;
            end            
            global eid;
            if isempty(eid)
                eid = int64(0);
            end
            if eid == intmax('uint64')
                eid = int64(1);
                warning('reset id');
            else
                eid = eid + 1;
            end
            this.Time = time;
            this.Type = type;
            this.Entity = entity;
            this.Identifer = eid;
            if nargin == 4
                this.Description = desc;
            end
        end
        
        % One can directly access _Entity.Type_.
        function st = get.EntityType(this)
            st = this.Entity.Type;
        end
        
        function name_str = get.Name(this)
            name_str = strcat(this.Entity.Type.char, this.Type.char);
        end
    end
    
    methods (Access = protected)
        function this = copyElement(ev)
            this = copyElement@matlab.mixin.Copyable(ev);
            %% Deep Copy Issue
            % *Entity*: is an exterior link. When performing copy, we should not make a copy of this
            % object. Instead, the link should be updated by the caller of the _copy_ function. To
            % secure the original data, we detach the link in the new copy from the original data.
            % *userdata*: whether it is a handle class object is unknown. If yes, we detach the link
            % to the original data instead of make a copy.
            if ~isempty(ev.Entity)
                this.Entity = ev.Entity.empty();
            end
            if ishandle(ev.userdata)
                this.userdata = [];
            end
        end
    end
end

