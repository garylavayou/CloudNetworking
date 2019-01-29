classdef NormalNetworkOptimizer < NetworkOptimizer	
	properties (SetAccess = protected)
		options;
	end
	
	properties (Access = protected)
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
				'AbsoluteTolerance', 10^-6, ...
				'OptimizeOrder', 1 ...
				));
			this.options = structupdate(this.options, options, ...
        {'RelativeTolerance', 'AbsoluteTolerance', 'OptimizeOrder'});
		end

	end
  
	%% Public Methods
	methods
		%% Optimization in Single Slice
		% Optimize the resource allocation in a single slice.
		%   All service flows are placed in one slice.
		%
		% options: [in-out]
		function [prices, results] = singleSliceOptimization( this, slices, options )
			% this.clearStates;
			net = this.hn;
			defaultopts = getstructfields(net.options, {'PricingPolicy'}, 'error');
			if nargin < 2
				options = defaultopts;
			else
				options = structmerge(defaultopts, options);
			end
			assert(net.options.SlicingMethod.IsSingle,...
				'error[%s]: unrecognized method (%s).', calledby(0), net.options.SlicingMethod.char);
			assert(isfield(net.options, 'PricingFactor'), ...
				'error[%s]: PricingFactor not specified, %s, %s',  calledby(0), ...
				'considerng provide it when creating the network',...
				'or specify it before calling this method.');
			if options.bCountTime
				t_start = tic; prt = 0; srt = 0;
			end
			
			%% Merge slices into one single big slice
			Nl = net.NumberLinks;
			Nn = net.NumberNodes;
			Ns = net.NumberSlices;
			slice_data.Adjacent = net.graph.Adjacent;
			slice_data.LinkMapS2P = (1:Nl)';
			slice_data.LinkMapP2S = (1:Nl)';
			slice_data.LinkCapacity = net.readLink('Capacity');
			slice_data.NodeMapS2P = (1:Nn)';
			slice_data.NodeMapP2S = (1:Nn)';
			slice_data.NodeCapacity = net.readDataCenter('Capacity');
			slice_data.FlowTable = table([],[],[],[],[],[],[],[], 'VariableNames',...
				{net.slices(1).FlowTable.Properties.VariableNames{:,:},'Weight', 'Type'});
			nf = 0;
			slice_data.VNFList = cell(Ns,1);
			slice_data.NumberPaths = 0;     
			for s = 1:Ns
				sl = net.slices(s);
				new_table = sl.FlowTable;
				% Map the virtual nodes to physical nodes.
				new_table.Source = sl.Nodes{new_table.Source, {'PhysicalNode'}};
				new_table.Target = sl.Nodes{new_table.Target, {'PhysicalNode'}};
				for f = 1:height(sl.FlowTable)
					% path_list is handle object, is should be copyed to the new table.
					path_list = PathList(sl.FlowTable{f,'Paths'});
					for p = 1:path_list.Width
						path_list{p}.node_list = sl.Nodes{path_list{p}.node_list,'PhysicalNode'};
					end
					new_table{f,'Paths'} = path_list;
				end
				new_table.Weight = sl.Weight*ones(height(new_table),1);
				new_table{:,'Type'} = s;   % replace 'flow_owner'
				slice_data.FlowTable = [slice_data.FlowTable; new_table];
				slice_data.VNFList{s} = sl.VNFList;
				slice_data.NumberPaths = max(slice_data.NumberPaths, sl.options.NumberPaths);     
				nf = nf + sl.NumberFlows;
			end
			slice_data.FlowPattern = FlowPattern.Default;
			slice_data.DelayConstraint = inf;
			
			slice_data = this.updateSliceData(slice_data, options);      % override by subclasses
			slice_data.Parent = net;
			slice_data.PricingPolicy = 'linear'; % the first step use the cost as price, so the policy is linear
			% the flow id and path id has been allocated in each slice already, no need to reallocate.
			ns = NormalSlice(slice_data);
			slice_data.Optimizer = 'SingleNormalSliceOptimizer';
			op = ns.getOptimizer(slice_data);
			op.setProblem('Price', options.UnitCost);
			options.slices = slices;
			op.optimalFlowRateSingleSlice(slice_data, options);
			if options.bCountTime
				t_stop = toc(t_start); prt = prt + t_stop; srt = srt + t_stop;
			end
			results = calculateOptimalOutput(this, ns);
			%% Compute the real resource demand with given prices
			% Individual slice adopt 'quadratic' pricing policy
			[prices, output] = pricingFactorAdjustment(net, options);
			if options.bCountTime
				output.runtime = output.runtime.Parallel + prt;
				output.runtime = output.runtime.Serial + srt;
				results.runtime = output.runtime;
			end
			%% Calculate the output
			results.SingleSlice = ns;
		end

	end
	
  methods (Access = {?NetworkOptimizer,?PhysicalNetwork})
    function slice_data = updateSliceData(this, slice_data, options)
      defaultopts = Dictionary('SlicingMethod', this.hn.options.SlicingMethod);
      if nargin <= 2
        options = defaultopts;
      else
        options = setdefault(options, defaultopts);
      end
      if options.SlicingMethod == SlicingMethod.SingleFunction
        error("[%s] not supported for ordered service chain.",calledby(0));
      end
    end

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
        ss.FlowTable.Weight.*fcnUtility(ss.Optimizer.getFlowRate(false))) ...
        - this.hn.totalCost(load);
      if DEBUG
        cprintf('Comments','Info: [%s] The optimal net social welfare of the network: %G.\n', ...
          calledby, output.WelfareOptimal);
      end
    end
  end
  
end

