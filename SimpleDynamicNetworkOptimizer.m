classdef SimpleDynamicNetworkOptimizer < SimpleNetworkOptimizer & IDynamicNetworkOptimizer
	%UNTITLED5 Summary of this class goes here
	%   Detailed explanation goes here
	
	properties
		
		
	end
	
	%% Constructor
	methods
		function this = SimpleDynamicNetworkOptimizer(net)
			this@SimpleNetworkOptimizer(net);
			this@IDynamicNetworkOptimizer(net);
		end
		
		
	end
end

