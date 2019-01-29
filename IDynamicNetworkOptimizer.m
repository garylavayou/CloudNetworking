classdef (HandleCompatible) IDynamicNetworkOptimizer 
	%UNTITLED5 Summary of this class goes here
	%   Detailed explanation goes here
	
	properties (Abstract, SetAccess = protected)
		options;
	end
	
	%% Constructor
	methods
		function this = IDynamicNetworkOptimizer(net, varargin)
			defaultopts = Dictionary(...
				'DiffNonzeroTolerance', 10^-3, ...  % 10^-4
				'ReconfigMethod', ReconfigMethod.Fastconfig ...
				);
			if length(varargin) >= 1
				defaultopts = structupdate(defaultopts, varargin{1});
			end
			this.options = setdefault(this.options, defaultopts);
		end
		
		
		
	end

	%%
	methods

	end
end

