classdef NormalNetworkOptimizer < NetworkOptimizer
	properties (SetAccess = protected)
		init_gamma_k double;
		init_q_k double;
	end	
  %% Constructor
	methods
		function this = NormalNetworkOptimizer(net, options)
			if nargin <= 1
				args = {};
			else
				args = {options};
			end
      this@NetworkOptimizer(net, args{:});
			
			this.options = setdefault(this.options, struct(...
				'RelativeTolerance', 10^-3, ...
				'AbsolueTolerance', 10^-6, ...
				'OptimizeOrder', 1 ...
				));
			this.options = structupdate(this.options, options);
		end

	end
  
	%% Public Methods
	methods
	end
	
  methods (Access = {?NetworkOptimizer,?PhysicalNetwork})
    function slice_data = updateSliceData(this, slice_data, options)
      defaultopts = Dictionary('SlicingMethod', this.hs.options.SlicingMethod);
      if nargin <= 2
        options = defaultopts;
      else
        options = setdefault(options, defaultopts);
      end
      if options.SlicingMethod == SlicingMethod.SingleNormal
        error("[%s] not supported for ordered service chain.",calledby(0));
      end
    end

  end
  
end

