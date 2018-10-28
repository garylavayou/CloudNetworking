classdef SubstrateAccessNetwork < CloudAccessNetwork & SubstrateNetwork
	%UNTITLED6 Summary of this class goes here
	%   Detailed explanation goes here
	
	methods
		function this = SubstrateAccessNetwork(varargin)
			this@CloudAccessNetwork(varargin{:});
			this@SubstrateNetwork();
		end
		
		
		function plot(this, b_undirect)
			if nargin <= 1
				b_undirect = false;
			end
			plot@CloudAccessNetwork(this, b_undirect);
			plot@SubstrateNetwork(this, b_undirect);
		end
		
		function [output, runtime] = optimizeResourcePrice1(this, slices)
			if nargin == 1
				arg_list = {this};
			else 
				arg_list = {this, slices};
			end
			if nargout <= 1
				output = optimizeResourcePrice1@SubstrateNetwork(arg_list{:});
			else
				[output, runtime] = optimizeResourcePrice1@SubstrateNetwork(arg_list{:});
			end
		end
		
		function [output, runtime] = optimizeResourcePriceNew(this, slices, options)
			switch nargin 
				case 1
					arg_list = {this};
				case 2
					arg_list = {this, slices};
				otherwise
					arg_list = {this, slices, options};
			end
			if nargout <= 1
				output = optimizeResourcePriceNew@SubstrateNetwork(arg_list{:});
			else
				[output, runtime] = optimizeResourcePriceNew@SubstrateNetwork(arg_list{:});
			end
		end
		
	end
	
	methods (Access = protected)
		function sl = createslice(this, slice_opt, varargin)
			sl = createslice@SubstrateNetwork(this, slice_opt);
		end
	end
end

