classdef SliceOptimizer < handle
	properties
    Variables;  % Variables after post processing.
		temp_vars = struct; % Temporary variables that directly resulted from optimization.
    prices;
    capacities;
		problem;
  end
  
  properties (SetAccess = protected)
		I_active_variables logical;  % indicator of active variables
		I_edge_path logical;    % Edge-Path Incidence Matrix
		I_dc_path logical;      % Node-Path Incidence Matrix
		I_flow_path logical;    % Flow-Path Incidence Matrix
    hs;         % reference to the host slice.
		x0;         % start point
		flow_rate;          % Temporary results of flow rate.
    options;
  end
  
  properties(Abstract)
    % Final results of VNF instance capacity on each node, this is configured by
		% inter/intra-slicing, the values are stored in |Variables.v|; This property is a
		% wrapper.
		VNFCapacity;
  end
  
  properties (Dependent)
    NumberVariables int32;	% Number of optimization variables.
    NumberLinearConstraints int32; % Number of linear constraints
  end
  
  methods (Abstract)
    clear(this);      % clear variables.
		%% priceOptimalFlowRate
		% Find the optimal flow rate that maximizing the net profit of the network slice.
		%
		%      [rate, net_profit, node_load, link_load] = priceOptimalFlowRate(this, x0, options)
		[utility, load] = priceOptimalFlowRate(this, x0, options);
		[profit,cost] = optimalFlowRate(this, new_opts);
		r = getFlowRate(this, ~);
		% ps = getFlowPaths(this, ~);
		setPathBandwidth(this, ~);
		profit = getProfit(this, options);
		v_n = getNodeLoad(this, isfinal, vars);
		ye = getLinkLoad(this, isfinal, vars);
    [tf, vars] = postProcessing(this);
	end
	
	methods (Static, Abstract)
		[profit, grad] = fcnSocialWelfare(x_vars, op, options);
    [profit, grad] = fcnProfit(vars, op, options);		% Objective function and gradient
    hs = fcnHessian(var_x, ~, op, options);
	end
  
  methods (Access = protected, Abstract)
    %% Override the method to provide property overriding
    n = get_number_variables(this);
    n = get_number_linear_constraints(this);
    % temp_vars = get_temp_variables(this);
    [x, fval] = optimize(this, options);
    
  end
  
  %% Constructor
  methods 
    function this = SliceOptimizer(slice, options)
      this.hs = slice;
			defaultopts = struct(...
				'NonzeroTolerance', 10^-3, ...
				'PostProcessing', 'round', ...
				'ConstraintTolerance', 10^-3, ...
				'Form', 'normal');
			if nargin >= 2
				this.options = structupdate(defaultopts, options);
			else
				this.options = defaultopts;
			end
			
			this.initializeState();
    end
  end
  
  %% Property Get Methods
  methods
		function n = get.NumberVariables(this)
			n = this.get_number_variables();
    end
    
    function n = get.NumberLinearConstraints(this)
			n = this.get_number_linear_constraints();
		end    
  end
  
  %% Public Methods
  methods				
    function update_options(this, options) %#ok<INUSD>

    end
    %%
    % save the results from parfor, the original results in the slice cannot be
    % retrived.
    function saveTempResults(this, result) 
      this.temp_vars = result.temp_vars;
      this.flow_rate = result.flow_rate;
      this.prices.Link = [];
      this.prices.Node = [];
      this.capacities = [];
      this.x0 = result.x0;
    end

    function setProblem(this, varargin) 
			for i = 1:2:(length(varargin)-1)
				switch varargin{i}
					case 'Price'
						p = varargin{i+1};
						if isempty(p)
							this.prices = struct('Link', [], 'Node', []);
						else
							% Raw price: corresponding price of links and nodes should be selected out.
							this.prices.Link = p.Link(this.hs.Links.PhysicalLink);
							this.prices.Node = p.Node(this.hs.getDCPI);
						end
					case 'Capacity'
						c= varargin{i+1};
						this.capacities.Link = c.Link(this.hs.Links.PhysicalLink);
						this.capacities.Node = c.Node(this.hs.getDCPI);
					case 'NodePrice'
						this.prices.Node = varargin{i+1};
					case 'LinkPrice'
						this.prices.Link = varargin{i+1};
					case 'Array'
					case 'Problem'
					case 'Index'
				end
			end
		end
    
	end
  
	methods (Access = protected)
		%% initializeState
		function initializeState(this)
			slice = this.hs;
			Nsn = slice.NumberServiceNodes;
			Np = slice.NumberPaths;
			Nve = slice.NumberLinks;
			Nf = slice.NumberFlows;
			
			this.I_dc_path = sparse(Nsn, Np);
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
					for k = 1:(path.Length-1)
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
	
	methods (Static, Access = protected)
		function interpretExitflag(exitflag, foutput)
			global DEBUG INFO;
			if nargin <= 1
				message = '';
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
					fprintf('Constraint violation: %f.\n', foutput.constrviolation);
					error('error: Abnormal exit (%d). %s', exitflag, message);
			end
		end
		
	end
end