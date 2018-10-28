classdef NormalDynamicNetworkOptimizer < NormalNetworkOptimizer & IDynamicNetworkOptimizer
	%UNTITLED5 Summary of this class goes here
	%   Detailed explanation goes here
	
	properties
		
		
	end
	
	%% Constructor
	methods
		function this = NormalDynamicNetworkOptimizer()
			this@NormalNetworkOptimizer();
			this@IDynamicNetworkOptimizer();
		end
		
		
	end

	methods
		function update_options(this, options)
			update_options@IDynamicSliceOptimizer(this, options);
			update_options@NormalSliceOpitmizer(this, options);
		end
	end

end

