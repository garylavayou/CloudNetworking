classdef Event < matlab.mixin.Copyable
    
    properties (SetAccess = protected)
        Type;           % EventType
        Time = 0;
        Identifer;      
        Entity;         % Entity handle can be checked with isvalid().
    end
    properties (Access = private)
        Description;
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
            persistent eid;
            if isempty(eid)
                eid = int64(1);
            else 
                if eid == intmax('uint64')
                    eid = int64(1);
                    warning('reset id');
                else
                    eid = eid + 1;
                end
            end
            if nargin == 0
                return;
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
    
end

