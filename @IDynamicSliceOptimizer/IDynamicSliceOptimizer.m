classdef (Abstract) IDynamicSliceOptimizer < handle
  properties (Constant)
    GLOBAL_OPTIONS Dictionary = Dictionary('theta0', 0.999);
  end

	properties(Constant)
		NUM_MEAN_BETA = 15;
		ENABLE_DYNAMIC_NORMALIZER = true;
		GET_BETA_METHOD = 'ExponetialMovingAverage';    % 'Average', 'ExponetialMovingAverage'; 'OverallAverage'
  end 
 
  properties
    % the variable |b_derive_vnf| decide if update VNF instance capacity.
    %    |b_derive_vnf=true|: derive VNF instance capacity from the optimization
    %    results of |Variables.z|;
    %    |b_derive_vnf=false|: apply |Variables.v| as VNF instance capacity.
    % For first time slice dimensioning,
    %   VNF capacity is the sum of VNF instance assigment (sum of z_npf);
    % For later reallocation (FSR2), the VNF capacity should be set as the optimized
    % variables (|this.Variables.v|, which might be larger than sum of |z_npf|, due to
    % reconfiguration cost).
    b_derive_vnf = true;
    b_dim = true;      % reset each time before reconfigurtion.
    %% define reconfiguraton cost for link and node variables.
    z_reconfig_cost;    % for z_nvf, used in optimization, the value updates during each fast reconfiguration
    x_reconfig_cost;    % for x_ef, used in optimization, the value updates during each fast reconfiguration
		vnf_reconfig_cost;
    old_state;
		changed_index;
		invoke_method = 0;
	end
	
	properties (Abstract)
		temp_vars;
		problem;	  % inlcuding 'numvar': the actuall number of each component variable
		pardata;    % pass data to parallel workers
	end
	
  properties (Abstract, SetAccess = protected)
		hs;
		num_vars;
    options;        % Options on performing optimizarion.
  end
  
	properties (Access = protected)
		a = 0.8;        % a for the history, should have a larger weight (EMA)
		%% options for slice dimensioning schedule algorithm
		% 'omega_upper':
		% 'omega_lower':
		% 'alpha':
		% 'series_length'
		% 'trend_length'
		sh_options Dictionary;
		%% data for slice dimensioning schdule algorithm
		% 'omegas':
		% 'profits'
		% 'omega_trend':
		% 'profit_trend':
		sh_data Dictionary;
		old_variables;      % last one configuration, last time's VNF instance capcity is the |v| field;
		lower_bounds = struct([]);
		upper_bounds = struct([]);
		topts;              % used in optimization, avoid passing extra arguments.
		max_flow_rate;
		diff_state;         % reset each time before reconfiguration.

    raw_beta;
		raw_cost;
		raw_costv;
  end

  methods (Abstract)
    [exitflag,fidx] = executeMethod(this, action);
    s = save_state(this);
    restore_state(this);
		onAddingFlow(this);
		onRemovingFlow(this);
    stat = get_reconfig_stat(this, stat_names);
	end
  
  methods (Abstract, Access = protected)
    identify_change(this, ~);
    [profit,cost] = fastReconfigure(this, action, options);
	end
	
	methods (Abstract, Access = {?IDynamicSliceOptimizer, ?DynamicNetwork})
		update_reconfig_costinfo(this, action, bDim);
	end
	
  methods
    function this = IDynamicSliceOptimizer(slice, slice_data)
      % Interval for performing dimensioning should be configurable.
      this.options.bDistributed = false;
			defaultopts = getstructfields(slice.Parent.Optimizer.options, { ...
				'DiffNonzeroTolerance', ...
				'ReconfigMethod' ...
				}, 'error');
			defaultopts = structmerge(defaultopts, getstructfields(slice.Parent.Optimizer.options, { ...
				'IntraSlicePenalty' ...  % for distributed algorithm
				}, 'ignore'));
			this.options = setdefault(this.options, defaultopts);
			this.options = structmerge(this.options, getstructfields(slice_data, {...
				'TimeInterval', ...
				'EventInterval', ...
				'Trigger', ...
				'bDistributed', ...
				}, 'ignore'));
			this.sh_data = Dictionary('omegas', [] , 'profits', [], ...
				'omega_trend', struct('ascend', [], 'descend', []),...
				'profit_trend', struct('ascend', [], 'descend', []),...
				'index', 0, 'reserve_dev', 0);
			this.sh_options = Dictionary('omega_upper', 0.97, 'omega_lower', 0.85, ...
				'alpha', 0.3, 'series_length', 40, 'trend_length', 20);
			
      switch this.options.ReconfigMethod
        case {ReconfigMethod.DimconfigReserve, ReconfigMethod.DimconfigReserve0,...
            ReconfigMethod.FastconfigReserve}
          this.options = structmerge(this.options, ...
            getstructfields(slice_data, 'bReserve', 'default', 2));
        case ReconfigMethod.DimBaseline
          this.options = structmerge(this.options, ...
            getstructfields(slice_data, 'bReserve', 'default', 0));
          this.sh_options.omega_upper = 1;
          this.sh_options.omega_lower = 0.9;
        case {ReconfigMethod.Dimconfig,ReconfigMethod.Fastconfig}
          this.options = structmerge(this.options, ...
            getstructfields(slice_data, 'bReserve', 'default', 0));
          this.sh_options.omega_upper = 1;
        otherwise
          this.options = structmerge(this.options, ...
            getstructfields(slice_data, 'bReserve', 'default', 1));
      end
      if this.options.ReconfigMethod == ReconfigMethod.DimconfigReserve
        this.sh_options = structmerge(this.sh_options, ...
          getstructfields(slice_data, 'UtilizationVariance', 'default', 0.05));
			else
				switch this.options.ReconfigMethod
					case {ReconfigMethod.DimBaseline, ReconfigMethod.Dimconfig}
						u = 0.2;
					case ReconfigMethod.DimconfigReserve0
						u = 0.1;
					otherwise
						u = 0;
				end
        this.sh_options.UtilizationVariance = u;
      end
      if this.options.bReserve
        if this.options.ReconfigMethod == ReconfigMethod.DimconfigReserve0
          this.options.Reserve = 0;
        else
          this.options = structmerge(this.options, ...
            getstructfields(slice_data, 'Reserve','default', 0.9));
        end
      end
      
      if isfield(this.options, 'EventInterval')
        slice.time.DimensionInterval = this.options.EventInterval/(2*slice.FlowArrivalRate); % both arrival and departure events
        slice.time.MinDimensionInterval = 1/5*slice.time.DimensionInterval;
      elseif isfield(this.options, 'TimeInterval')
        slice.time.DimensionInterval = slice.options.TimeInterval;
      elseif isfield(this.options, 'Trigger')
        slice.time.DimensionInterval = 10/slice.FlowArrivalRate;
      end
      slice.time.DimensionIntervalModified = slice.time.DimensionInterval;
			if ~IDynamicSliceOptimizer.ENABLE_DYNAMIC_NORMALIZER
				if ~isfield(slice_data, 'ReconfigScaler')
					this.options.ReconfigScaler = 1;
				else
					this.options.ReconfigScaler = slice_data.ReconfigScaler;
				end
			end    
			
    end
  
	end
  
	%% Public Methods
  methods		
		function initializeParallel(this, procedure, options) %#ok<INUSD>
			if isfield(this.hs.old_state, 'link_capacity')
				this.old_variables.c = this.hs.old_state.link_capacity;
			end
			if strcmpi(procedure, 'optimalFlowRate')
				this.pardata.FlowRate = this.hs.FlowTable.Rate;
			end
		end
		
		function setProblem(this, varargin)
			for i = 1:2:(length(varargin)-1)
				switch varargin{i}
					case 'Problem'
						if ~isempty(varargin{i+1})
							prbm = this.problem;
							if isfield(prbm, 'invoke_method') % only valid in parallel mode.
								this.invoke_method = prbm.invoke_method;
								this.problem = rmfield(prbm, 'invoke_method');
							end
						end
				end
			end
		end
		
    function update_options(this, options)
      if isfield(options, 'ResidualCapacity')
				slice = this.hs;
				res_caps.Link = options.ResidualCapacity.Link(slice.Links.PhysicalLink);
				res_caps.Node = options.ResidualCapacity.Node(slice.getDCPI());
        this.options.ResidualCapacity = res_caps;
      end
		end
		
		%% Get reconfiguration cost
		% model = {'linear'|'const'|'none'}
		% Return: reconfiguration cost, number of reconfiguration, ration of
		% reconfigurations, and total number of variables.
		function [total_cost, reconfig_cost] = get_reconfig_cost(this, model, isfinal)
			if nargin <= 1
				model = 'const';
			end
			if strcmpi(model, 'const')
				isfinal = true;
			elseif strcmpi(model, 'linear') && nargin <= 2
				isfinal = false;
			end
			tol_vec = this.options.DiffNonzeroTolerance;
			s = this.diffstate(isfinal);
			
      if strcmpi(model, 'const')
        % logical array cannot be used as the first argument of dot.
        reconfig_cost.x = dot(this.x_reconfig_cost, s.diff_x_norm>tol_vec);
        reconfig_cost.z = dot(this.z_reconfig_cost, s.diff_z_norm>tol_vec);
        if isfield(s, 'diff_v_norm')
          reconfig_cost.v = dot(this.vnf_reconfig_cost, s.diff_v_norm>tol_vec);
        end
      elseif strcmpi(model, 'linear')
        reconfig_cost.x = dot(abs(s.diff_x), this.x_reconfig_cost);
        reconfig_cost.z = dot(abs(s.diff_z), this.z_reconfig_cost);
        if isfield(s, 'diff_v')
          reconfig_cost.v = dot(abs(s.diff_v), this.vnf_reconfig_cost);
        end
        % The linear cost coefficient (*_reconfig_cost) is the original
        % one, which is not scaled. In orde to compute the scaled cost as the
        % problem formulation, we additionally mutiply the scaler.
        if this.ENABLE_DYNAMIC_NORMALIZER
          beta = this.getbeta();
          reconfig_cost.x = beta.x * reconfig_cost.x;
          reconfig_cost.z = beta.z * reconfig_cost.z;
          if isfield(s, 'diff_v')
            reconfig_cost.v = beta.v * reconfig_cost.v;
          end
        else
          reconfig_cost.x = this.options.ReconfigScaler * reconfig_cost.x;
          reconfig_cost.z = this.options.ReconfigScaler * reconfig_cost.z;
          if isfield(s, 'diff_v')
            reconfig_cost.v = this.options.ReconfigScaler * reconfig_cost.v;
          end
        end
      else
        error('%s: invalid cost model <%s>.', calledby, model);
      end
      vars = fieldnames(reconfig_cost);
      total_cost = 0;
      for i = 1:length(vars)
        total_cost = total_cost + reconfig_cost.(vars{i});
      end
		end
		
		function init_reconfig_info(this)
			global DEBUG g_results;
			slice = this.hs;
			if this.b_derive_vnf
				% We initialize the reconfiguration cost coefficient when the slice is
				% added. The coefficients is required to compute the cost after the
				% optimizaton at <DynamicCloudNetwork.onAddingSlice>.
				% After redimensioning, the former cost coefficients are still needed to
				% calculate the reconfiguration cost. So in that case, we update it before
				% performing reconfigure/redimensioning.
				this.b_derive_vnf = false;
				this.update_reconfig_costinfo('add', false, 1); % the fourth argument can be any value.
				%         this.x_reconfig_cost = (this.I_edge_path)' * slice.Links.ReconfigCost;
				%         this.z_reconfig_cost = repmat(slice.ServiceNodes.ReconfigCost, ...
				%           slice.NumberPaths*slice.NumberVNFs, 1);
				this.update_reconfig_costvinfo(1);
				%         this.vnf_reconfig_cost = slice.Parent.options.VNFReconfigCoefficient * ...
				%           repmat(slice.ServiceNodes.ReconfigCost, slice.NumberVNFs, 1);
				this.old_variables = Dictionary(this.Variables);
				slice.save_state;
				if this.ENABLE_DYNAMIC_NORMALIZER
					%% Initialize L1-approximation normalizer
					% change of magnitude:
					%   x -> 0; 0 -> x:    |x|      : small portion
					%   x -> x1:           |x-x1|   : large portion
					% we choose beta to be (1/2) of the magnitude.
					field_names = {'x', 'z', 'v'};
					cost_entries = {'x_reconfig_cost', 'z_reconfig_cost', 'vnf_reconfig_cost'} ;
					if strcmpi(this.GET_BETA_METHOD, 'LeastSquare')
						for i = 1:2
							this.raw_cost.const.(field_names{i}) = ...
								dot(this.(cost_entries{i}), this.Variables.(field_names{i})~=0);
							this.raw_cost.linear.(field_names{i}) = ...
								dot(this.(cost_entries{i}), this.Variables.(field_names{i}))/2;
						end
						this.raw_costv.const = ...
							dot(this.vnf_reconfig_cost, this.Variables.v~=0);
						this.raw_costv.linear = ...
							dot(this.vnf_reconfig_cost, this.Variables.v)/4;
						this.raw_costv.bInit = true;
					else
						for i =1:3
							cost = dot(this.(cost_entries{i}), this.Variables.(field_names{i}))/2;
							cost0 = dot(this.(cost_entries{i}), this.Variables.(field_names{i})~=0);
							this.raw_beta.(field_names{i}) = cost0/cost;
						end
					end
				end
				
				stat = slice.get_reconfig_stat();
				% options.bFinal = true;
				stat.Profit = this.getProfit();
				stat.ReconfigType = ReconfigType.Dimensioning;
				stat.ResourceCost = slice.getCost();   % _optimizeResourcePriceNew_ use 'quadratic-price'
				stat.FairIndex = (sum(slice.FlowTable.Rate))^2/(slice.NumberFlows*sum(slice.FlowTable.Rate.^2));
				if this.options.ReconfigMethod >= ReconfigMethod.Dimconfig
					this.sh_data.profits = stat.Profit;
					this.sh_data.omegas = stat.Utilization;
				end
				if ~isempty(DEBUG) && DEBUG
					disp(stat);
				end
				g_results = stat;       % The first event.
				this.b_dim = 0;
			else
				%% Reach here when redimensioning the slice
				% Reduce reconfigurations by perform one more FSR after DSR.
				% Post processing has been done to variables {x,z,v}.
				%
				if this.options.ReconfigMethod >= ReconfigMethod.Dimconfig && ...
						this.options.ReconfigMethod < ReconfigMethod.DimBaseline% && false
					% not necessary: the variables {x,z} will be updated in FSR, while others
					% including {v} is unchanged. And {v} will be used later, see
					% <SimpleDynamicSlice>.<finalize>.
					old_num_vars = this.num_vars; % keep the DSR variables for later use.
					old_pardata = this.pardata;
					old_prbm = this.problem;
					this.temp_vars = rmstructfields(this.temp_vars, {'x', 'z', 'tx', 'tz'});
					this.update_reconfig_costinfo(this.sh_options.action);
					new_opt = struct;
					if this.options.ReconfigMethod == ReconfigMethod.DimconfigReserve
						new_opt.bEnforceReserve = true;
					end
					this.max_flow_rate = slice.FlowTable.Rate;
					this.fastReconfigure(this.sh_options.action, new_opt);
					this.num_vars = old_num_vars;
					this.pardata = old_pardata;
					this.problem = old_prbm;
				end
			end
			[omega, sigma_o] = slice.utilizationRatio();
			if this.options.ReconfigMethod >= ReconfigMethod.Dimconfig
				this.sh_options.sigma = sigma_o + var_sigma(omega);
			end
			%% Nest function: fiiting the variation of utilization
			% fitting y = ax^2+bx+c
			function s = var_sigma(omega, b) %#ok<INUSD>
				%                 A = [1 1 1;
				%                     0.9^2 0.9 1;
				%                     0.8^2 0.8 1;
				%                     0.5^2 0.5 1;
				%                     0 0 1];
				%                 if nargin == 1
				%                     if this.NumberFlows <= 60
				%                         b = [0.3; 0.2; 0.15; 0.1; 0.3];
				%                     elseif this.NumberFlows <= 120
				%                         b = [0.2; 0.14; 0.1; 0.05; 0.2];
				%                     elseif this.NumberFlows <= 180
				%                         b = [0.18; 0.12; 0.09; 0.04; 0.18];
				%                     elseif this.NumberFlows <= 240
				%                         b = [0.15; 0.1; 0.08; 0.03; 0.15];
				%                     else
				%                         b = [0.12; 0.09; 0.075; 0.025; 0.12];
				%                     end
				%                 end
				%                 c = A\b;
				%                 f = @(x) c(1)*x^2+c(2)*x+c(3);
				%                 s = f(omega);
				%%
				s = 0.2*sqrt(100/slice.NumberFlows);
			end
		end
		
	end
	
	methods (Access = protected)
    [profit, exitflag, fidx] = DimensioningReconfigure( this, action, new_opts );
		%% handle_zero_flow
		% called by fast slice reconfigure procedure.
		function [profit, cost, output] = handle_zero_flow(this, ~)
			cost = this.hs.getCost();
			profit = -cost;
			this.clear();
			this.hs.Links{:,'Load'} = 0;
			this.hs.ServiceNodes{:,'Load'} = 0;
			%% Capacity/Prices Not Change
			% When the number of flows reduced to 0, we do not change the VNF instance
			% capcity, so that the reconfiguration cost is reduced. However, at the
			% next time when a flow arrive, the reconfiguration cost might be larger.
			% Another method is to set the VNF instance capacity to 0.
			%                 if nargout >= 1
			%                     this.VNFCapacity = 0;
			%                 end
			if nargout >= 3
				output = struct();
			end
		end
		
		% update reconfiguration cost with scaler.
		% If dimensioning is applied, the orginal cost is firstly updated.
		function update_reconfig_costvinfo(this, ~)
			slice = this.hs;
			this.vnf_reconfig_cost = slice.Parent.options.VNFReconfigCoefficient * ...
				repmat(slice.ServiceNodes.ReconfigCost, slice.NumberVNFs, 1);
			if nargin >= 2
				return;
			end
			if this.ENABLE_DYNAMIC_NORMALIZER
				beta = this.getbeta();
				this.topts.vnf_reconfig_cost = beta.v*this.vnf_reconfig_cost;
			else
				this.topts.vnf_reconfig_cost = ...
					this.options.ReconfigScaler*this.vnf_reconfig_cost;
			end
		end
		
		%%
		% dynamic adjust the normalizer for L1-approximation of the reconfiguration cost.
		% This method is used when 'ENABLE_DYNAMIC_NORMALIZER' is set to true.
		function postl1normalizer(this)
			[~, cost] = this.get_reconfig_cost('const');
			[~, linear_cost] = this.get_reconfig_cost('linear', true);
			field_names = {'const', 'linear'};
			b = this.getbeta();
			if strcmpi(this.GET_BETA_METHOD, 'LeastSquare')
				arr_cost = [cost, linear_cost];
				fn = fieldnames(arr_cost(2));
				for j = 1:length(fn)
					arr_cost(2).(fn{j}) = arr_cost(2).(fn{j})/b.(fn{j});
				end
				for i = 1:2
					if length(this.raw_cost.(field_names{i})) < this.NUM_MEAN_BETA
						this.raw_cost.(field_names{i})(end+1) = ...
							getstructfields(arr_cost(i), {'x','z'});
					else
						this.raw_cost.(field_names{i}) = ...
							[this.raw_cost.(field_names{i})(2:end), ...
							getstructfields(arr_cost(i), {'x','z'})];
					end
				end
				if isfield(arr_cost, 'v')
					if arr_cost(1).v < 10^-3 || arr_cost(2).v < 10^-3
						return;
					end
					if this.raw_costv.bInit
						this.raw_costv.const = arr_cost(1).v;
						this.raw_costv.linear = arr_cost(2).v;
						this.raw_costv.bInit = false;
						return;
					end
					for i = 1:2
						if length(this.raw_costv.(field_names{i})) < this.NUM_MEAN_BETA
							this.raw_costv.(field_names{i})(end+1) = arr_cost(i).v;
						else
							this.raw_costv.(field_names{i}) = ...
								[this.raw_costv.(field_names{i})(2:end), arr_cost(i).v];
						end
					end
				end
			else
				if this.b_dim
					field_names = {'x','z','v'};
				else
					field_names = {'x','z'};
				end
				for i = 1:length(field_names)
					if isfield(cost, field_names{i}) && isfield(linear_cost, field_names{i})
						if cost.(field_names{i})<10^-6
							% beta.(field_names{i}) = this.raw_beta.(field_names{i})(end)/2;
							continue;
						else
							f = cost.(field_names{i})/linear_cost.(field_names{i});
							if f > 2
								f = 2;
							elseif f<1/2
								f = 1/2;
							end
							beta.(field_names{i}) = b.(field_names{i})*(1+f)/2;
						end
						if length(this.raw_beta.(field_names{i})) < this.NUM_MEAN_BETA
							this.raw_beta.(field_names{i})(end+1) = beta.(field_names{i});
						else
							this.raw_beta.(field_names{i}) = ...
								[this.raw_beta.(field_names{i})(2:end), beta.(field_names{i})];
						end
					end
				end
			end
			
			%             function [f,g] = fcn_optbeta(beta,cc,ca)
			%                 f = (norm(beta*ca-cc,2))^2;
			%                 g = sum(2*(beta*ca-cc).*ca);
			%             end
			%             function h = hess_optbeta(beta, lambda, ca) %#ok<INUSL>
			%                 h = sum(2*ca.^2);
			%             end
		end
		
		function b = getbeta(this)
			field_names = {'x','z','v'};
			if strcmpi(this.GET_BETA_METHOD, 'LeastSquare')
				for i=1:2
					b.(field_names{i}) = sum([this.raw_cost.const.(field_names{i})].*[this.raw_cost.linear.(field_names{i})])./ ...
						sum([this.raw_cost.linear.(field_names{i})].^2);
				end
				b.v = sum(this.raw_costv.const.*this.raw_costv.linear)./ ...
					sum(this.raw_costv.linear.^2);
			else
				if strcmpi(this.GET_BETA_METHOD, 'ExponentiaLMovingAverage')
					a0 = 0.25;
					for i=1:2
						%                 b.(field_names{i}) = mean(this.raw_beta.(field_names{i}));
						n = length(this.raw_beta.(field_names{i}));
						alpha = zeros(n,1);
						alpha(2:n) = a0 * (1-a0).^(n-2:-1:0);
						alpha(1) = (1-a0)^(n-1);
						b.(field_names{i}) = dot(alpha, this.raw_beta.(field_names{i}));
					end
				else
					for i=1:2
						b.(field_names{i}) = mean(this.raw_beta.(field_names{i}));
					end
				end
				% v_nf >= sum_{p}{z_npf}
				n_path = nnz(this.I_dc_path)/this.hs.NumberServiceNodes;
				b.v = b.z/n_path;
			end
		end

		
		%% calculate the difference of slice state
		% Arguments:
		%   # isfinal: If the slice is not in the final stage, the difference should be
		%     update with each call. Otherwise, the difference can be saved for later use.
		%     see also <Slice.isFinal>;
		function ds = diffstate(this, isfinal)
			if isfinal && ~isempty(this.diff_state)
				if nargout == 1
					ds = this.diff_state;
				end
				return;
			end
			
			% At the beginging, |old_variables| equals to |Variables|.
			% |old_variables| is still valid after adding/removing flows.
			% |old_variables| have more elements than |this.Variables.x| when removing
			% flows.
			if length(this.old_variables.x) <= length(this.temp_vars.x)  % only ==
				if isfinal
					new_x = this.Variables.x;
				else
					new_x = this.temp_vars.x;
				end
				%%
				% In the optimization, we ignore the reconfiguration cost of the new flow,
				% while in the results, we should consider the new flow's reconfiguration
				% cost. See also <update_reconfig_costinfo>.
			else
				new_x = sparse(numel(this.old_variables.x),1);
				if isfinal
					new_x(~this.changed_index.x) = this.Variables.x;
				else
					new_x(~this.changed_index.x) = this.temp_vars.x;
				end
				%% The removed flow has been considered with the reconfguration cost
				% The removed flow's reconfiguration cost is constant and does not appear
				% in the optimization, see also <update_reconfig_costinfo>.
			end
			ds.diff_x = sparse(new_x-this.old_variables.x);
			mid_x = 1/2*(abs(new_x)+abs(this.old_variables.x));
			nz_index_x = mid_x~=0;
			ds.diff_x_norm = sparse(length(ds.diff_x),1);
			% Here, we assume that all tiny variables (smaller than NonzeroTolerance) have been
			% eleminated. So that the following operation can identify changes.
			ds.diff_x_norm(nz_index_x) = abs(ds.diff_x(nz_index_x)./mid_x(nz_index_x));
			
			if length(this.old_variables.z) <= length(this.temp_vars.z)
				if isfinal
					new_z = this.Variables.z;
				else
					new_z = this.temp_vars.z;
				end
			else
				new_z = sparse(numel(this.old_variables.z),1);
				if isfinal
					new_z(~this.changed_index.z) = this.Variables.z;
				else
					new_z(~this.changed_index.z) = this.temp_vars.z;
				end
			end
			ds.diff_z = sparse(new_z-this.old_variables.z);
			mid_z = 1/2*(abs(new_z)+abs(this.old_variables.z));
			nz_index_z = mid_z~=0;
			ds.diff_z_norm = sparse(length(ds.diff_z),1);
			ds.diff_z_norm(nz_index_z) = abs(ds.diff_z(nz_index_z)./mid_z(nz_index_z));
			%%%
			% Reconfiguration cost of VNF capacity.
			% No reconfiguration of VNF instance for 'fastconfig'.
			if isfield(this.temp_vars, 'v')
				if length(this.old_variables.v) <= length(this.temp_vars.v)
					if isfinal
						new_vnf_capacity = this.getVNFCapacity();
					else
						new_vnf_capacity = this.temp_vars.v;
					end
				else
					new_vnf_capacity = zeros(length(this.old_variables.v),1);
					if isfinal
						new_vnf_capacity(~this.changed_index.v) = this.getVNFCapacity();
					else
						new_vnf_capacity(~this.changed_index.v) = this.temp_vars.v;
					end
				end
				ds.diff_v = new_vnf_capacity-this.old_variables.v(:);
				ds.mid_v = 1/2*(abs(new_vnf_capacity)+abs(this.old_variables.v(:)));
				nz_index_v = ds.mid_v~=0;
				ds.diff_v_norm = zeros(length(ds.diff_v),1);
				ds.diff_v_norm(nz_index_v) = abs(ds.diff_v(nz_index_v)./ds.mid_v(nz_index_v));
			end
			if isfinal
				this.diff_state = ds;
			end
		end
		
		
	end % end protected methods
end