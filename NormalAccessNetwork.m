classdef NormalAccessNetwork < IAccessNetwork & NormalCloudNetwork
	%UNTITLED6 Summary of this class goes here
	%   Detailed explanation goes here
	
	methods
		function this = NormalAccessNetwork(varargin)
			this@NormalCloudNetwork();
			this@IAccessNetwork(varargin{:});
		end
		
		
		function plot(this, b_undirect)
			if nargin <= 1
				b_undirect = false;
			end
			plot@NormalCloudNetwork(this, b_undirect);
		end
		
		function [output, varargout] = optimizeResourcePrice1(this, varargin)
			[output, varargout{1:nargout-1}] = optimizeResourcePrice1@NormalCloudNetwork(this, varargin{:});
		end
		
		%		[output, runtime] = optimizeResourcePrice(this, slices, options)
		function [output, varargout] = optimizeResourcePrice(this, varargin)
			[output, varargout{1:nargout-1}] = optimizeResourcePrice@NormalCloudNetwork(this, varargin{:});
		end
		
	end
	
	methods (Access = protected)
		function sl = createslice(this, slice_opt, varargin)
			sl = createslice@NormalCloudNetwork(this, slice_opt);
			sl.getOptimizer(slice_opt);
		end
	end
end

