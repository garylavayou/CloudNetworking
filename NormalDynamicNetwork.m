classdef NormalDynamicNetwork < DynamicNetwork
	
	methods
		function this = NormalDynamicNetwork(node_opt, link_opt, VNF_opt, net_opt)
			this@DynamicNetwork(node_opt, link_opt, VNF_opt, net_opt);
		end
	end
	
	methods
		function op = getOptimizer(this, options)
			if nargin == 1
				this.op = NormalDynamicNetworkOptimizer(this);
			else
				this.op = NormalDynamicNetworkOptimizer(this, options);
			end
			op = this.op;
		end
	end
	
	methods (Access=protected)
		function sl = createslice(this, slice_opt, varargin)
			% examine flow arrival parameters.
			% usage of <Slice>: if a slice without |ArrivalRate| or |ServiceInterval| or
			% their values are invalid, the slice (<Slice> or <DynamicSlice>) is treated
			% as no dynamics of flow, and it will not handle flow events. Therefore
			% initilize it as class <Slice> is OK.
			if isfield(slice_opt, 'ClassName')
				sl = instantiateclass(slice_opt.ClassName, ...
					rmfield(slice_opt, 'ClassName'), varargin{:});
			else
				if ~isfield(slice_opt, 'ArrivalRate') || ~isfield(slice_opt, 'ServiceInterval')
					sl = NormalSlice(slice_opt);
				elseif isempty(slice_opt.ArrivalRate) || isempty(slice_opt.ServiceInterval)
					sl = NormalSlice(slice_opt);
					warning('slice created with type SimpleSlice.');
				else
					sl = NormalDynamicSlice(slice_opt);
				end
			end
			this.slices(end+1) = sl;
			sl.getOptimizer(slice_opt);
			sl.FlowTable{:, 'Type'} = FlowType.Normal;  % Specify flow type
		end
	end
end

