classdef SliceEntity < Entity
    %UNTITLED7 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(SetAccess = ?SliceEntityBuilder)
        SliceIdentifier;
    end
    properties (Access = protected)
        Child;          % see also <FlowEntity>.
        seed;
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
    
end

