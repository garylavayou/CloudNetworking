%% Link delay specification
% How to determine the link latency if the latency data is not available.
% Options: 'Bandwidth', 'BandwidthInverse', 'Constant', 'Random'
%%
classdef LinkDelayOption < uint32
    
    properties
    end
    
    methods
    end
    
    enumeration
        BandwidthPropotion(1),
        BandwidthInverse(2),
        Constant(3),
        Random(4)
    end
    
end

