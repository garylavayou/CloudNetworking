classdef (Abstract) EventReceiver < matlab.mixin.Copyable
    
    properties
        event = struct('Count', 0);     %statistics of events
				time = struct;
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

