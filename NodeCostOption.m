classdef NodeCostOption < uint32
    %UNTITLED8 Summary of this class goes here
    %   Detailed explanation goes here
    
    enumeration
        Uniform(1);            % All nodes have the same unit cost comparable to link cost.
        Weighted(2);           % The node cost is mapped from link cost, while the mapping weight is user specified.
        CapacityInverse(3);    % unit node cost inversely proportional to the capacity.
        NetworkSpecified(4);   % Network oprator specifies the node cost.
        Random(5);
    end
    
    methods
    end
    
end

