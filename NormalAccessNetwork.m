classdef NormalAccessNetwork < IAccessNetwork & NormalCloudNetwork
	%UNTITLED6 Summary of this class goes here
	%   Detailed explanation goes here
	
	methods
		function this = NormalAccessNetwork(varargin)
			[node_opt, link_opt] = NormalAccessNetwork.loadNetworkData(varargin{1:2});
			varargin{1} = node_opt;
			varargin{2} = link_opt;
			this@NormalCloudNetwork(varargin{:});
			this@IAccessNetwork(node_opt);
		end
		
		function value = readNode(this, name, varargin)
			value = readNode@IAccessNetwork(this, name, varargin{:});
			if isempty(value)
				value = readNode@NormalCloudNetwork(this, name, varargin{:});
			end
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
	methods(Static)
		function [node_opt, link_opt] = loadNetworkData(node_opt, link_opt)
			[node_opt, link_opt] = loadNetworkData@IAccessNetwork(node_opt, link_opt);
		end
	end
	
	methods (Access = protected)
		function slice_opt = preAddingSlice(this, slice_opt)
			slice_opt = preAddingSlice@NormalCloudNetwork(this, slice_opt);
			slice_opt = preAddingSlice@IAccessNetwork(this, slice_opt);
		end

		function sl = createslice(this, slice_opt, varargin)
			sl = createslice@NormalCloudNetwork(this, slice_opt);
			sl.getOptimizer(slice_opt);
		end

		function end_points = generateEndPoints(this, slice_opt)
			end_points = generateEndPoints@IAccessNetwork(this, slice_opt);
			if isempty(end_points)
				end_points = generateEndPoints@NormalCloudNetwork(this, slice_opt);
			end
		end
	end
end

