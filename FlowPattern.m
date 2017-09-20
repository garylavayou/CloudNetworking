classdef FlowPattern < uint32
    
    properties
    end
    
    methods
    end
    
    enumeration
        Default(0);
        RandomSingleFlow(1);
        RandomMultiFlow(2);
        SingleSourceMultiFlow(3);
        MultiSourceMultiFlow(4);
        RandomInterDataCenter(5);
        RandomInterBaseStation(6);
        RandomDataCenter2BaseStation(7);
    end
end

%% Flow Pattern
% * *RandomSingleFlow*
% * *RandomMultiFlow*
% * *SingleSourceMultiFlow*
% * *MultiSourceMultiFlow*

