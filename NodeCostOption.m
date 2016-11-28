classdef NodeCostOption < uint32
    %UNTITLED8 Summary of this class goes here
    %   Detailed explanation goes here
    
    enumeration
        Uniform(1);            % All nodes have the same unit cost comparable to link cost.
        NetworkSpecified(2);   % Network oprator specifies the node cost.
    end
    
    methods
    end
    
end

