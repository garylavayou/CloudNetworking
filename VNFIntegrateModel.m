classdef VNFIntegrateModel < uint32
    
    properties
    end
    
    methods
    end
    
    enumeration
        AllInOne(1);   % All VNF instances at one node are encapsulated in one VM
        SameSliceInOne(2); % VNF instances in one node from the same slice are encapsulated in one VM
        SameTypeInOne(3); % The VNF instances of the same type on one node are encapsulated in one VM
        Separated(4); % Each VNF instance on one node is encapsulated as separated VM
    end
    
end

