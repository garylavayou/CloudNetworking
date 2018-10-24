classdef SliceOptimizer
	properties
    Variables;  % Variables after post processing.
		prices;
  end
  
  properties (Access = protected)
    hs;         % reference to the host slice.
		x0;         % start point
		temp_vars = struct; % Temporary variables that directly resulted from optimization.
		flow_rate;          % Temporary results of flow rate.
		I_active_variables logical;  % indicator of active variables
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
  end
  
  methods (Access = protected, Abstract)
    %% Override the method to provide property overriding
    n = get_number_variables(this);
    n = get_number_linear_constraints(this);
    % temp_vars = get_temp_variables(this);
    [x, fval] = optimize(this, params, options);
    [tf, vars] = postProcessing(this);
  end
  
  %% Constructor
  methods 
    function this = SliceOptimizer(slice, options)
      this.hs = slice;
      this.options = getstructfields(options, ...
        {'NonzeroTolerance', 'PostProcessing', 'ConstraintTolerance'}, 'ignore');
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
    
    %% priceOptimalFlowRate
    % Find the optimal flow rate that maximizing the net profit of the network slice.
    %
    %      [rate, net_profit, node_load, link_load] = priceOptimalFlowRate(this, x0, options)
    [utility, load] = priceOptimalFlowRate(this, x0, options);
    [profit,cost] = optimalFlowRate(this, new_opts);
  end
  
end