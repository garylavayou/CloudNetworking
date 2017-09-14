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
        BandwidthPropotion(1),  % Static delay
        BandwidthInverse(2),    % Dynamic delay
        Constant(3),            % Static
        Random(4)               % Static
    end
    
end

