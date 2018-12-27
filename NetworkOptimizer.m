classdef NetworkOptimizer < handle
  %UNTITLED5 Summary of this class goes here
  %   Detailed explanation goes here
  
  properties (SetAccess = protected)
    hn;
    runtime;
		options Dictionary;
  end
  
  %% Constructor
  methods
    function this = NetworkOptimizer(net, options)
			defaultopts = Dictionary(...
				'NonzeroTolerance', 10^-3, ...
				'PostProcessing', 'round', ...
				'ConstraintTolerance', 10^-3, ...
				'Form', 'normal', ...
				'Threshold', 'min', ...
				'InterSlicePenalty', 1, ...
				'OptimizationTool', 'matlab' ... % {matlab|cvx}
				);
			if nargin >= 2
				this.options = structupdate(defaultopts, options);
			else
				this.options = defaultopts;
			end
			if isfield(options, 'IntraSlicePenalty')
				this.options.IntraSlicePenalty = options.IntraSlicePenalty;
			end
      this.hn = net;
    end

  end
	
	methods (Abstract)
	end
  
  methods (Access = {?NetworkOptimizer,?PhysicalNetwork})
    slice_data = updateSliceData(this, slice_data, options);
    
		function prices = initPrice(this, slices, unitcost, options)
			t1 = 1;           % {0.1|0.8|1}
			if nargin >= 4 && isfield(options, 'InitPrice')
				link_usage = false(this.hn.NumberLinks,1);
				node_usage = false(this.hn.NumberDataCenters,1);
				for i = 1:length(slices)
					link_usage(slices(i).Links.PhysicalLink) = true;
					node_usage(slices(i).getDCPI) = true;
				end
				if find([options.InitPrice.Link(link_usage)==0;options.InitPrice.Node(node_usage)==0],1)
					prices.Link = t1* unitcost.Link;
					prices.Node = t1* unitcost.Node;
					%         init_price.Link = prices.Link
					%         init_price.Node = prices.Node;
				else
					prices.Link = options.InitPrice.Link;
					prices.Node = options.InitPrice.Node;
				end
			else
				prices.Link = t1* unitcost.Link;
				prices.Node = t1* unitcost.Node;
			end
		end
  end
end

