classdef NormalNetworkOptimizer < NetworkOptimizer
  %% Constructor
	methods
		function this = NormalNetworkOptimizer(net)
      this@NetworkOptimizer(net)
		end

  end
  
  methods (Access = {?PhysicalNetwork})
    function slice_data = updateSliceData(~, slice_data, options)
      defaultopts = struct('SlicingMethod', SlicingMethod.SingleNormal);
      if nargin <= 2
        options = defaultopts;
      else
        options = structupdate(defaultopts, options);
      end
      if options.SlicingMethod == SlicingMethod.SingleNormal
        warning("[%s] not supported for ordered service chain.",calledby(0));
      end
    end

  end
  
end

