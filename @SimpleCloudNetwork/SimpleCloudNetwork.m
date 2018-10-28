%% Cloud Network
% Enable network resource allocation, by mechanisms such as pricing.
%
% In the network, forwarding nodes only take charge of packet processing. Forwarding nodes
% may connect to Data Center, and then VNF instances can be created in data center. The
% connection between forwarding nodes and data centers are assumed with infinity bandwidth
% and zero latency.
%
% The topology of CloudNetwork is predefined.
%%
classdef SimpleCloudNetwork < PhysicalNetwork
	%% Constructor
	methods
		%%
		% * *options*:
		%       _PricingFactor_ is used for <singleSliceOptimization> and <staticSlicing>.
		%       _Threshold_ is used for resource pricing.
		%       _Method_ is used for selecting method to solve sub-problem.
		%
		% NOTE: only put common options ('SlicingMethod', 'Form', etc.)in the constructor. Put
		% those method-specific options to the correspongding method.
		function this = SimpleCloudNetwork(varargin)
			this@PhysicalNetwork(varargin{:});
      this.op = SimpleNetworkOptimizer(this);
		end
	end
	
	%% Public Methods
	methods
    %     function update_options(this, opt_name, opt_value)
    %       rmidx = false(length(opt_name,1),1);
    %       for i = 1:length(opt_name)
    %         if contains(opt_name{i}, {'SlicingMethod', 'PricingFactor'})
    %           this.options.(opt_name{i})= opt_value{i};
    %           rmidx(i) = true;
    %         end
    %       end
    %       opt_name(rmidx) = [];
    %       opt_value(rmidx) = [];
    %       update_options@PhysicalNetwork(this, opt_name, opt_value);
    %     end
	end
	
	methods (Access=protected)	
		function sl = createslice(this, slice_opt, varargin)
			this.slices{end+1} = SimpleSlice(slice_opt);
			sl = this.slices{end};
    end
		
		% Now the same as <PhysicalNetwork>, subclasses might override it.
		%         function [flow_table, phy_adjacent, flag] = ...
		%                 generateFlowTable(this, graph, slice_opt)
		%             [flow_table, phy_adjacent, flag] = ...
		%                 this.generateFlowTable@PhysicalNetwork(graph, slice_opt);
		%         end
  end
	
	methods
		[output, runtime] = optimizeResourcePrice1(this, init_price);
		[output, runtime] = optimizeResourcePriceNew(this, init_price, sub_slices);
		function [output, runtime] = singleSliceOptimization(this, new_opts)
			if nargin <= 1
				argins = {this};
			else
				argins = {this, new_opts};
			end
			if nargout <= 1
				output = this.op.singleSliceOptimization(argins{:});
			else
				[output, runtime] = this.op.singleSliceOptimization(argins{:});
			end
		end
		output = StaticSlicing(this, slice);
		[tb, stbs] = saveStatTable(PN, output, rt, slice_types, method);
	end
	
	methods (Static)
		[stat, slice_stat] = createStatTable(num_point, num_type, type);
	end
end
