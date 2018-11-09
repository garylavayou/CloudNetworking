classdef NetworkOptimizer < handle
  %UNTITLED5 Summary of this class goes here
  %   Detailed explanation goes here
  
  properties
    hn;
    runtime;
		options;
  end
  
  %% Constructor
  methods
    function this = NetworkOptimizer(net, options)
			defaultopts = struct(...
				'NonzeroTolerance', 10^-3, ...
				'PostProcessing', 'round', ...
				'ConstraintTolerance', 10^-3, ...
				'Form', 'normal', ...
				'Threshold', 'min');
			if nargin >= 2
				this.options = structupdate(defaultopts, options);
			else
				this.options = defaultopts;
			end
			
      this.hn = net;
    end

  end
	
	methods (Abstract)
	end
  
  methods (Access = {?NetworkOptimizer,?PhysicalNetwork})
    slice_data = updateSliceData(this, slice_data, options);
    
    function prices = initPrice(this, unitcost) %#ok<INUSL>
      t1 = 1;           % {0.1|0.8|1}
      prices.Link = t1* unitcost.Link;
      prices.Node = t1* unitcost.Node;
    end
		
  end
end

