classdef SliceOptimizer < handle & matlab.mixin.Heterogeneous
	properties (Abstract)
		temp_vars Dictionary;  % Temporary variables that directly resulted from optimization.
		problem Dictionary;	   % inlcuding 'numvar': the actual number of each component variable
		pardata Dictionary;    % pass data to parallel workers
	end
	properties (Abstract, SetAccess = protected)
		options Dictionary;
	end
	properties
    Variables Dictionary;  % Variables after post processing.
    prices Dictionary;
    capacities Dictionary;
	end
  	
  properties (SetAccess = protected)
    hs;									% reference to the host slice.
		num_vars;		% the raw number of each component variable
		I_active_variables logical;  % indicator of active variables
		I_edge_path logical;    % Edge-Path Incidence Matrix
		I_dc_path logical;      % Node-Path Incidence Matrix
		I_flow_path logical;    % Flow-Path Incidence Matrix
		x0;         % start point
		flow_rate;          % Temporary results of flow rate.
	end
  
	properties (Dependent)
		Host;
	end
	
  methods (Abstract)
		%% priceOptimalFlowRate
		% Find the optimal flow rate that maximizing the net profit of the network slice.
		%
		%      [rate, net_profit, node_load, link_load] = priceOptimalFlowRate(this, x0, options)
		[utility, load] = priceOptimalFlowRate(this, x0, options);
		[profit, cost, output] = optimalFlowRate(this, new_opts);
		r = getFlowRate(this, ~);
		% ps = getFlowPaths(this, ~);
		profit = getProfit(this, options);
		v_n = getNodeLoad(this, isfinal, vars);
		ye = getLinkLoad(this, isfinal, vars);
		vc = getVNFCapacity(this, z);
    [tf, vars] = postProcessing(this);
		setPathBandwidth(this, ~);
	end
	
	methods (Static, Abstract, Access = protected)
		[profit, grad] = fcnSocialWelfare(x_vars, op, options);
    [profit, grad] = fcnProfit(vars, op, options);		% Objective function and gradient
    hs = fcnHessian(var_x, ~, op, options);
	end
  
  methods (Access = protected, Abstract)
    [x, fval] = optimize(this, options);
  end
  
  %% Constructor
  methods 
    function this = SliceOptimizer(slice, options)
      this.hs = slice;
			op = slice.Parent.Optimizer;
			assert(~isempty(op), 'error: network optimizer not initialized.');
			defaultopts = getstructfields(op.options, {...
				'NonzeroTolerance', ...
				'PostProcessing', ...
				'ConstraintTolerance', ...
				'Form', ...
				'OptimizationTool', ... % {matlab|cvx}
				'InterSlicePenalty' ...
				}, 'error');
			if nargin >= 2
				defaultopts = structupdate(defaultopts, options);
			end
			this.options = Dictionary();
			this.options = setdefault(this.options, defaultopts);
			this.Variables = Dictionary();  % Variables after post processing.
			this.temp_vars = Dictionary(); % Temporary variables that directly resulted from optimization.
			this.prices = Dictionary();
			this.capacities = Dictionary();
			this.problem = Dictionary();	% inlcuding 'numvar': the actuall number of each component variable
			this.pardata = Dictionary();
		end
		
	end
  
	%% Property Access Methods
	methods
		%% Host
		function h = get.Host(this)
			h = this.hs;
		end
		function set.Host(this, hs)
			if isempty(hs) && isa(hs, 'Slice')
				this.hs = hs;
			elseif hs.Optimizer == this
				this.hs = hs;
			else
				error('error: host slice and optimizer do not match.');
			end
		end
		
	end

  %% Public Methods
  methods
		function initializeParallel(this, procedure, ~)
			slice = this.hs;
			this.pardata = Dictionary( ...
				'AggregateNodeUsage', slice.Parent.AggregateNodeUsage(slice.getSNPI()), ...
				'AggregateLinkUsage', slice.Parent.AggregateLinkUsage(slice.Links.PhysicalLink), ...
				'PhysicalLinkCapacity', slice.Parent.readLink('Capacity', slice.Links.PhysicalLink),...
				'PhysicalNodeCapacity', slice.Parent.readDataCenter('Capacity', slice.getDCPI()), ...
				'ProcessEfficiency', slice.Parent.VNFTable.ProcessEfficiency(slice.VNFList), ...
				'NumberLinks', slice.NumberLinks, ...
				'NumberPaths', slice.NumberPaths, ...
				'NumberVNFs', slice.NumberVNFs, ...
				'NumberFlows', slice.NumberFlows, ...
				'NumberServiceNodes', slice.NumberServiceNodes, ...
				'Weight', slice.Weight, ...
				'PathOwner', slice.path_owner, ...
				'SlicingMethod', slice.options.SlicingMethod, ...
				'PricingPolicy', slice.options.PricingPolicy ...
				);
			if strcmpi(procedure, 'optimalFlowRate')
				this.pardata.NodeCapacity = slice.ServiceNodes.Capacity;
				this.pardata.LinkCapacity = slice.Links.Capacity;
				this.pardata.ParentNumberVNFs = slice.Parent.NumberVNFs;
				if isempty(this.pardata.Weight)
					this.pardata.FlowWeight = slice.FlowTable.Weight;
				end
			end
		end
			
		function clear(this)      % clear variables.
      this.Variables.erase();
			this.temp_vars.erase();
		end
		
    function update_options(this, options) %#ok<INUSD>

    end
    %%
    % save the results from parfor, the original results in the slice cannot be
    % retrived.
    function saveTempResults(this, result) 
      this.temp_vars = Dictionary(result.temp_vars);
      this.x0 = result.x0;
			% FOR DEBUG
			% this.setPathBandwidth(this.temp_vars.x);
    end

    function setProblem(this, varargin) 
			for i = 1:2:(length(varargin)-1)
				switch varargin{i}
					case 'Price'
						p = varargin{i+1};
						if isempty(p)
							this.prices.Link = [];
							this.prices.Node = [];
						else
							% Raw price: corresponding price of links and nodes should be selected out.
							this.prices.Link = p.Link(this.hs.Links.PhysicalLink);
							this.prices.Node = p.Node(this.hs.getDCPI());
						end
					case 'Capacity'
						c= varargin{i+1};
						if isempty(c)
							this.capacities.Link = [];
							this.capacities.Node = [];
						else
							this.capacities.Link = c.Link(this.hs.Links.PhysicalLink);
							this.capacities.Node = c.Node(this.hs.getDCPI());
						end
					case 'NodePrice'
						this.prices.Node = varargin{i+1};
					case 'LinkPrice'
						this.prices.Link = varargin{i+1};
					case 'Array'
					case 'Problem'
						if isempty(varargin{i+1})
							this.problem.erase();
						else 
							this.problem = varargin{i+1};
							prbm = this.problem(1);
							if isfield(prbm, 'num_full_vars')  % only valid in parallel mode.
								this.num_vars = prbm.num_full_vars;
								prbm = rmfield(prbm, 'num_full_vars'); 
							end
							if isfield(prbm, 'I_active_variables') % only valid in parallel mode.
								this.I_active_variables = prbm.I_active_variables;
								prbm = rmfield(prbm, 'I_active_variables');  %#ok<NASGU>
							end
						end
					case 'Index'
					case 'ParallelData'
						this.pardata = Dictionary(varargin{i+1});
				end
			end
		end
    
		%% initializeState
		function initializeState(this)
			slice = this.hs;
			Nsn = slice.NumberServiceNodes;
			Np = slice.NumberPaths;
			Nve = slice.NumberLinks;
			Nf = slice.NumberFlows;
			
			this.I_dc_path = sparse(Nsn, Np);  % MOVE to SimpleSliceOptimizer
			this.I_edge_path = sparse(Nve, Np);
			this.I_flow_path = sparse(Nf, Np);
			slice.path_owner = zeros(Np,1);
			% this.local_path_id = zeros(this.NumberPaths, 1);
			pid = 0;
			for fid=1:Nf
				path_list = slice.FlowTable{fid, 'Paths'};
				for j = 1:path_list.Width
					pid = pid + 1;
					this.I_flow_path(fid,pid) = 1;
					slice.path_owner(pid) = fid;
					path = path_list{j};
					path.local_id = pid;    % record the local path in the slice.
					for k = 1:path.Length
						e = path.Link(k);
						eid = slice.graph.IndexEdge(e(1),e(2));
						this.I_edge_path(eid, pid) = 1;
						dc_index = slice.Nodes{e(1),'ServiceNode'};
						if dc_index~=0
							this.I_dc_path(dc_index, pid) = 1;
						end
					end
					dc_index = slice.Nodes{e(2),'ServiceNode'}; % last node
					if dc_index~=0
						this.I_dc_path(dc_index, pid) = 1;
					end
				end
			end
		end
	
	end
  
	methods
		[payment, grad, pseudo_hess] = fcnLinkPricing(this, link_price, link_load, unit);
		[payment, grad, pseudo_hess] = fcnNodePricing(this, node_price, node_load, unit);
	end

	
	methods (Static, Access = protected)
		function interpretExitflag(exitflag, foutput)
			global DEBUG INFO;
			if nargin <= 1
				message = '';
			elseif ischar(foutput)
				message = strtok(foutput, newline);
			else
				message = strtok(foutput.message, newline);
			end
			switch exitflag
				case 0
					if ~isempty(DEBUG) && DEBUG
						warning(message);    % max-iteration number exceed.
					elseif ~isempty(INFO) && INFO
						cprintf('SystemCommands', '%s\n', message);
					end
				case 1
					if ~isempty(INFO) && INFO
						fprintf('(%d) %s\n', exitflag, message);
					end
				case 2
					if ~isempty(DEBUG) && DEBUG
						warning('%s(%d)', message, exitflag);
					elseif ~isempty(INFO) && INFO
						cprintf('Comment', '%s(%d)\n', message, exitflag);
					end
				case -3
					error('error: Objective function unbounded below (%d). %s', exitflag, message);
				otherwise
					if nargin >= 2 && isstruct(foutput)
						fprintf('Constraint violation: %f.\n', foutput.constrviolation);
					end
					error('error: Abnormal exit (%d). %s', exitflag, message);
			end
		end
		
	end
end