classdef SlimpleDynamicNetworkOptimizer < SimpleNetworkOptimizer & IDynamicNetworkOptimizer
	methods
		function this = SlimpleDynamicNetworkOptimizer(net, varargin)
			this@SimpleNetworkOptimizer(net, varargin{:});
			this@IDynamicNetworkOptimizer(net, varargin{:});
		end
	end
end