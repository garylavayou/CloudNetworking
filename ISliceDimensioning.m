classdef (Abstract,HandleCompatible) ISliceDimensioning
	% Interface: slice dimensioning.
	%   Resource allocation for multiple coexisting slices.
		
	methods (Abstract)
		initialize(this);
		finalize(this);
		isFinal(this);
		r = getFlowRate(this);
		c = getLinkCapacity(this, isfinal);
		c = getNodeCapacity(this, isfinal);
		vc = getVNFCapacity(this, z);
		sc = getCost(this, load, model);
		r = getRevenue(this);
		b = checkFeasible(this, vars, opt_opts);
		
	end
end

