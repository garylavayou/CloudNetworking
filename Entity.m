classdef Entity < matlab.mixin.Heterogeneous & matlab.mixin.Copyable
    
    properties (SetAccess = {?Entity,?EntityBuilder})
        EntityIdentifier uint64;   % Global identifer of entities, might not be used.
        ArriveTime=-1;
        ServiceTime=-1;
    end

    properties (Constant, Access = private)
        sn = SerialNumber;  
        % NOTE: Each type of class need to define their own serial number;
        % If a serieral number is shared among multiple classes, we can define a common
        % super class for them.
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
        function this = Entity(time_arrive, time_serve, src, varargin)
            % Replace <persistent>, thus we can backup and recover the identifier out of the class.
            % TODO: Replace the 'entity_id' as a global SerialNumber instance.
            if nargin == 0 
                return;
            end
            if nargin >= 4 
                desc = varargin{1};
            end
            len_ta = length(time_arrive);
            len_ts = length(time_serve);
            len_s = length(src);
            len_args = [len_ta, len_ts, len_s];
            len = max(len_args); 
            if find((len_args<len) & (len_args>1))
                error('%s: The length of arguments mismatch.', calledby);
            end
            %% Handle Heterogeneous Array
            % Copy operation will not call the constructor again.
            % The results will be same type, not the hetero-array.
            % But for handle classes, the copy operation only duplicate the handle, not
            % the object it reference. So we need to call the copy method when perform
            % initialization.
            % NOTE: <repelem https://www.mathworks.com/help/matlab/ref/repelem.html> and
            % <repmat https://www.mathworks.com/help/matlab/ref/repmat.html>
            % have the same effect for scalar. 
            if len > 1
                this = repelem(this, len, 1);
                for i = 2:len
                    this(i,1) = copy(this(1));
                end
            end
            id = Entity.sn.next(len);
            for i = len:-1:1
                this(i).ArriveTime = time_arrive(min(i, len_ta));
                this(i,1).EntityIdentifier = id(i);
            end
            if ~isempty(src)
                for i = len:-1:1
                    this(i,1).Builder = src(min(i,len_s));
                end
            end
            if ~isempty(time_serve)
                for i = len:-1:1
                    this(i,1).ServiceTime = time_serve(min(i,len_ts));
                end
            end
            if exist('description', 'var')
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
            tf = false(size(this));
            for i = 1:numel(this)
                if isempty(this(i).ArriveTime) || isempty(this(i).ServiceTime) || ...
                        this(i).ArriveTime < 0 || this(i).ServiceTime < 0
                    tf = true;
                end
            end
        end
    end
    
    methods(Static)
        function n = getGlobalEntityId()
            n = Entity.sn.ID;
        end
        function setGlobalEntityId(value)
            h = Entity.sn;
            h.ID = value;
        end
    end
end

