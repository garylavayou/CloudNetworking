classdef LinkCostOption <uint32
    %LinkCostOption Unit cost.
    
    enumeration
        Uniform(1);             % All links have the same unit cost 1.
        LengthDependent(2);     % Unit cost depends on the link length
        NetworkSpecified(3);    % Network oprator specifies the link cost
        CapacityInverse(4);
    end
    
    methods
    end
    
end

