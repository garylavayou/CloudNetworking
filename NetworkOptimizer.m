classdef NetworkOptimizer < handle
  %UNTITLED5 Summary of this class goes here
  %   Detailed explanation goes here
	properties (Abstract, SetAccess = protected)
		options Dictionary;
	end
	properties (Access = {?NetworkOptimizer, ?PhysicalNetwork})
    runtime = 0;
		iterations = 1;
	end
  properties (SetAccess = protected)
    hn;
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
				'InterSlicePenalty', 0.1, ...
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
			if isfield(options, 't1')
				t1 = options.t1;
			else
				t1 = 1;           % {0.1|0.8|1}
			end
			unitcost = this.hn.convertParameter(unitcost, 'vector');
			if nargin >= 4 && isfield(options, 'InitPrice')
				Ne = this.hn.NumberLinks;
				res_usage = false(Ne+this.hn.NumberDataCenters,1);
				for i = 1:length(slices)
					res_usage(slices(i).Links.PhysicalLink) = true;
					res_usage(Ne+slices(i).getDCPI()) = true;
				end
				init_price = convertParameter(options.InitPrice, 'vector');
				if find(init_price(res_usage)==0,1)
					prices = t1* unitcost;
				else
					prices = init_price;
				end
			else
				prices = t1* unitcost;
			end
		end
  end
end
