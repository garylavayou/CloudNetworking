classdef SimpleNetworkOptimizer < NetworkOptimizer
	%UNTITLED5 Summary of this class goes here
	%   Detailed explanation goes here
	properties (SetAccess = protected)
		options Dictionary;
	end
	%% Constructor
	methods
		function this = SimpleNetworkOptimizer(net, varargin)
      this@NetworkOptimizer(net, varargin{:});
		end

  end
  
  methods
		[prices, results] = singleSliceOptimization(this, new_opts);
  end

  methods (Access = protected)
    %% This function use intermediate results.
		% called by <singleSliceOptimization>.
    function output = calculateOptimalOutput(this, ss)
      global DEBUG;
      if ~exist('DEBUG', 'var')
        DEBUG = false;
      end
      
      load = this.hn.getNetworkLoad(ss, struct('Stage', 'temp'));
      output.WelfareOptimal = sum(...
        ss.FlowTable.Weight.*fcnUtility(ss.Optimizer.getFlowRate(ss.Optimizer.temp_vars))) ...
        - this.hn.totalCost(load);
      if DEBUG
        cprintf('Comments','Info: [%s] The optimal net social welfare of the network: %G.\n', ...
          calledby, output.WelfareOptimal);
      end
		end
		
	end
  
  methods (Access = {?NetworkOptimizer, ?PhysicalNetwork})
    function slice_data = updateSliceData(this, slice_data, options)
      defaultopts = Dictionary(...
				'SlicingMethod', this.hn.options.SlicingMethod,...
				'Form', this.options.Form);
      if nargin <= 2
        options = defaultopts;
      else
        options = setdefault(options, defaultopts);
      end
      % options.Stage = 'temp';
      if options.SlicingMethod == SlicingMethod.SingleNormal
        switch options.Form
					case 'compact'
						b_vnf = false(this.hn.NumberVNFs, 1);
						for s = 1:this.hn.NumberSlices
							b_vnf(this.hn.slices(s).VNFList) = true;
							if isempty(find(b_vnf==false,1))
								break;
							end
						end
						slice_data.VNFList = find(b_vnf);
					case 'normal'
						slice_data.VNFList = 1:this.hn.NumberVNFs;
					otherwise
						error('error: un-recognized value for ''Form''.');
        end
      end
		end
		
	end
end

