classdef (Abstract) EventReceiver < handle & matlab.mixin.Heterogeneous
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods(Static, Abstract)
        eventhandler(this, source, eventData);
    end
    
    methods(Sealed)
        function tf = eq(b, a)
            tf = eq@handle(b,a);
        end
    end
    
%     methods (Access = protected, Abstract)
%         tf = compare(a,b);
%     end
end

