classdef NetworkOptimizer
  %UNTITLED5 Summary of this class goes here
  %   Detailed explanation goes here
  
  properties
    hn;
    runtime;
  end
  
  %% Constructor
  methods
    function this = NetworkOptimizer(net)
      this.hn = net;
    end

  end
	
	methods (Abstract)
	end
	
  methods (Abstract, Access = {?PhysicalNetwork})
		runtime = saveRuntime(this);
		
		[sp_profit, b_violate, violates] = SolveSCPCC(this, slices, prices, capacities, options);

  end
  
  methods (Access = {?PhysicalNetwork})
    slice_data = updateSliceData(this, slice_data, options);
    
    function prices = initPrice(this, unitcost)
      t1 = 1;           % {0.1|0.8|1}
      prices.Link = t1* unitcost.Link;
      prices.Node = t1* unitcost.Node;
    end
		
  end
end

