classdef (Abstract, HandleCompatible) EventReceiver 
    
    properties (Access = protected)
        event = struct('Count', 0);     %statistics of events
    end
    
    methods(Abstract)
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

