classdef FlowEventData < DispatchEventData
    
    properties
        slice;
        flow;       % flow table entries.
    end
    
    methods
        function this = FlowEventData(ev, slice, flow, userdata)
            if nargin <= 3
                userdata = [];
            end
            this@DispatchEventData(ev, userdata);
            this.slice = slice;
            this.flow = flow;
        end
    end
    
end

