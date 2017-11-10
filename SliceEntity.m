classdef SliceEntity < Entity
    %UNTITLED7 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(SetAccess = ?SliceEntityBuilder)
        SliceIdentifier;
    end
    properties (Access = protected)
        seed;
    end
    properties (SetAccess = {?SliceEntity,?RandomEventDispatcher})
        Child;          % see also <FlowEntity>.
    end
    properties(Dependent)
        Options;
    end
    
    methods (Access = ?SliceEntityBuilder)
        function this = SliceEntity(time_arrive, time_serve, src, varargin)
            this@Entity(time_arrive,time_serve,src,varargin{2:end});            
            if length(varargin) >= 1
                this.seed = varargin{1};
            end
            %
            % The identifier can be update in the constructor.
        end
    end
    
    methods
        % You can also directly access _Builder.Options_.
        function opt = get.Options(this)
            opt = this.Builder.Options;
            opt.RandomSeed = this.seed;
        end
    end
    
    methods (Access = protected)
        function this = copyElement(et)
            this = copyElement@Entity(et);
            %% Deep Copy Issue
            % *Child* is an exterior link. When performing copy, we should not make a copy of this
            % object. Instead, the link should be updated by the caller of the _copy_ function. To
            % secure the original data, we detach the link in the new copy from the original data.
            % See also <Entity>.
            if ~isempty(et.Child)
                this.Child = et.Child.empty();
            end
        end
    end
end

