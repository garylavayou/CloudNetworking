classdef SimpleNetworkOptimizer < NetworkOptimizer
	%UNTITLED5 Summary of this class goes here
	%   Detailed explanation goes here
	
	%% Constructor
	methods
		function this = SimpleNetworkOptimizer(net)
      this@NetworkOptimizer(net)
		end

  end
  
  methods
		[output, runtime] = singleSliceOptimization(this, new_opts);
  end

  methods (Access = protected)
    %% This function use intermediate results.
    function output = calculateOptimalOutput(this, ss)
      global DEBUG;
      if ~exist('DEBUG', 'var')
        DEBUG = false;
      end
      
      load = this.getNetworkLoad(ss, struct('Stage', 'temp'));
      output.WelfareOptimal = sum(...
        ss.FlowTable.Weight.*fcnUtility(ss.getFlowRate(ss.temp_vars.x))) ...
        - this.totalCost(load);
      if DEBUG
        cprintf('Comments','Info: [%s] The optimal net social welfare of the network: %G.\n', ...
          calledby, output.WelfareOptimal);
      end
    end
  end
  
  methods (Access = {?NetworkOptimizer, ?PhysicalNetwork})
    function slice_data = updateSliceData(this, slice_data, options)
      defaultopts = struct('SlicingMethod', SlicingMethod.SingleNormal);
      if nargin <= 2
        options = defaultopts;
      else
        options = structupdate(defaultopts, options);
      end
      % options.Stage = 'temp';
      if options.SlicingMethod == SlicingMethod.SingleNormal
        if options.bCompact
          b_vnf = false(this.hn.NumberVNFs, 1);
          for s = 1:this.hn.NumberSlices
            b_vnf(this.hn.slices{s}.VNFList) = true;
            if isempty(find(b_vnf==false,1))
              break;
            end
          end
          slice_data.VNFList = find(b_vnf);
        else
          slice_data.VNFList = 1:this.hn.NumberVNFs;
        end
      end
    end

    function prices = initPrice(this, slices, unitcost, options)
      t1 = 1;           % {0.1|0.8|1}
      if nargin >=4 && isfield(options, 'InitPrice')
        link_usage = false(this.hn.NumberLinks,1);
        node_usage = false(this.hn.NumberDataCenters,1);
        for i = 1:length(slices)
          link_usage(slices{i}.Links.PhysicalLink) = true;
          node_usage(slices{i}.getDCPI) = true;
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

    
		[output, runtime] = singleSliceOptimization( this, new_opts );
  end
end

