classdef NormalDynamicNetworkOptimizer < NormalNetworkOptimizer & IDynamicNetworkOptimizer
	%UNTITLED5 Summary of this class goes here
	%   Detailed explanation goes here
	
	properties
		
		
	end
	
	%% Constructor
	methods
		function this = NormalDynamicNetworkOptimizer(varargin)
			this@NormalNetworkOptimizer(varargin{:});
			this@IDynamicNetworkOptimizer(varargin{:});
		end
		
		
	end

	methods
	end

end

