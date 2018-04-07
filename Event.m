classdef Event < matlab.mixin.Copyable
    
    properties (SetAccess = protected)
        Type EventType;           % EventType
        Time = 0;
        Identifer uint64;      
    end
    
    properties (Constant, Access = private)
        sn = SerialNumber;  
        % NOTE: Each type of class need to define their own serial number;
        % If a serieral number is shared among multiple classes, we can define a common
        % super class for them.
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
            len_t = length(time);
            len_ty = length(type);
            len_e = length(entity);
            len_args = [len_t, len_ty, len_e];
            len = max(len_args); 
            if find((len_args<len) & (len_args>1))
                error('%s: The length of arguments mismatch.', calledby);
            end
            id = Event.sn.next(len);
            for i = len:-1:1
                this(i,1).Identifer = id(i);
            end
            if ~isempty(time)
                for i = len:-1:1
                    this(i,1).Time = time(min(i,len_t));
                end
            end
            if ~isempty(type)
                for i = len:-1:1
                    this(i,1).Type = type(min(i, len_ty));
                end
            end
            if ~isempty(entity)
                for i = len:-1:1
                    this(i,1).Entity = entity(min(i, len_e));
                end
            end
            if nargin == 4
                if ischar(desc)
                    desc = {desc};
                end
                len_d = length(desc);
                if len_d == 1
                    for i = len:-1:1
                        this(i).Description = desc{1};
                    end
                elseif len_d > 1 
                    len_d = min(len_d, len);
                    for i = len_d:-1:1
                        if ischar(desc{i})
                            this(i).Description = desc{i};
                        else
                            warning('%s: the description is not a character string.', calledby);
                        end
                    end
                end
                
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
    
    methods(Static)
        function n = getGlobalEventId()
            n = Event.sn.ID;
        end
        function setGlobalEventId(value)
            h = Event.sn;
            h.ID = value;
        end
    end
end

