classdef SimpleDynamicSliceOptimizer < SimpleSliceOptimizer & IDynamicSliceOptimizer
  
  properties (Access=protected)
    %% define reconfiguraton cost for link and node variables.
    z_reconfig_cost;    % for z_npf, used in optimization, the value updates during each fast reconfiguration
    x_reconfig_cost;    % for x_p, used in optimization, the value updates during each fast reconfiguration
  end
	properties(Dependent, Access=protected)
		num_varv;
  end
  
  %% Constructor
  methods
    function this = SimpleDynamicSliceOptimizer(slice, options)
      this@SimpleSliceOptimizer(slice);
      this@IDynamicSliceOptimizer(slice, options);
    end
  end
  
  %% Property Get Methods
  methods
		function n = get.num_varv(this)
			n = this.NumberServiceNodes*this.NumberVNFs;
    end
  end
  
  %% Public Methods
  methods
		[exitflag,fidx] = executeMethod(this, action);
    
    function clear(this)
      clear@SimpleSliceOptimizer(this);
      this.Variables.v = [];
      this.temp_vars.v = [];
    end
		%% TODO: only update the colums arriving/removing
		% since the resource changes after each slice configuration, the matrix needs be
		% computed each time. Therefore, we define the access function to evalue it.
		% This decoupt the computation of incident matrix and As_res, in
		% <initializeState>.
		%         function As = getAs_res(this)
		%         end
 
    function update_options(this, options)
      update_options@IDynamicSliceOptimizer(this, options);
      update_options@SimpleSliceOpitmizer(this, options);
    end
    
    function s = save_state(this)
      this.old_state.vnf_capacity = this.VNFCapacity;
      this.old_state.I_dc_path = this.I_dc_path;
      this.old_state.I_edge_path = this.I_edge_path;
      this.old_state.I_flow_path = this.I_flow_path;
      this.old_state.path_owner = this.path_owner;
      this.old_state.variables = this.Variables;
      this.old_state.x_reconfig_cost = this.x_reconfig_cost;
      this.old_state.z_reconfig_cost = this.z_reconfig_cost;
      this.old_state.vnf_reconfig_cost = this.vnf_reconfig_cost;
      this.old_state.As_res = this.As_res;
      if nargout >= 1
        s = this.old_state;
      end
    end
    
    function restore_state(this, s)
      if nargin <= 1
        s = this.old_state;
      end
      this.op.I_dc_path = s.I_dc_path;
      this.op.I_edge_path = s.I_edge_path;
      this.op.I_flow_path = s.I_flow_path;
      this.op.path_owner = s.path_owner;
      this.op.Variables = s.Variables;
      this.op.x_reconfig_cost = s.x_reconfig_cost;
      this.op.op.z_reconfig_cost = s.z_reconfig_cost;
      this.op.vnf_reconfig_cost = s.vnf_reconfig_cost;
      this.op.As_res = s.As_res;
    end
    
    function init_reconfig_info(this)
			global DEBUG g_results;
      if this.b_derive_vnf
        % We initialize the reconfiguration cost coefficient when the slice is
        % added. The coefficients is required to compute the cost after the
        % optimizaton at <DynamicCloudNetwork.onAddingSlice>.
        % After redimensioning, the former cost coefficients are still needed to
        % calculate the reconfiguration cost. So in that case, we update it before
        % performing reconfigure/redimensioning.
        this.b_derive_vnf = false;
        this.x_reconfig_cost = (this.I_edge_path)' * this.hs.Links.ReconfigCost;
        this.z_reconfig_cost = repmat(this.hs.ServiceNodes.ReconfigCost, ...
          this.hs.NumberPaths*this.hs.NumberVNFs, 1);
        this.vnf_reconfig_cost = this.hs.Parent.options.VNFReconfigCoefficient * ...
          repmat(this.hs.ServiceNodes.ReconfigCost, this.hs.NumberVNFs, 1);
        this.old_variables = this.Variables;
        this.hs.get_state;
        this.max_flow_rate = this.hs.FlowTable.Rate;
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
        
        stat = this.hs.get_reconfig_stat();
        options = getstructfields(this.options, 'PricingPolicy', 'default', {'quadratic'});
        options.bFinal = true;
        stat.Profit = this.hs.getProfit(options);
        stat.ReconfigType = ReconfigType.Dimensioning;
        stat.ResourceCost = this.hs.getCost(options.PricingPolicy);   % _optimizeResourcePriceNew_ use 'quadratic-price'
        stat.FairIndex = (sum(this.FlowTable.Rate))^2/(this.NumberFlows*sum(this.FlowTable.Rate.^2));
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
        if this.options.ReconfigMethod >= ReconfigMethod.Dimconfig && ...
            this.options.ReconfigMethod < ReconfigMethod.DimBaseline% && false
          this.temp_vars = rmstructfields(this.temp_vars, {'x','z','tx','tz'});
          this.update_reconfig_costinfo(this.sh_options.action);
          new_opt = struct;
          if this.options.ReconfigMethod == ReconfigMethod.DimconfigReserve
            new_opt.bEnforceReserve = true;
          end
          this.fastReconfigure(this.sh_options.action, new_opt);
        end
      end
      [omega, sigma_o] = this.utilizationRatio();
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
        s = 0.2*sqrt(100/this.NumberFlows);
      end
    end
    
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
      options = getstructfields(this.hs.Parent.options, ...
        {'DiffNonzeroTolerance', 'NonzeroTolerance'});
      tol_vec = options.DiffNonzeroTolerance;
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
      total_cost = reconfig_cost.x + reconfig_cost.z;
      if isfield(reconfig_cost, 'v')
        total_cost = total_cost + reconfig_cost.v;
      end
    end
    
    [utility, node_load, link_load] = priceOptimalFlowRate(this, x0, new_opts);
    
    %% optimalFlowRate
    % |new_opts|:
    % * *FixedCost*: used when slice's resource amount and resource prices are fixed,
    % resource cost is fixed. When this options is specified, the base method
    % <Slice.optimalFlowRate> returns profit without resource cost.
    % See also <DynamicSlice.fcnSocialWelfare>,<DynamicSlice.optimize>.
    function [profit, cost] = optimalFlowRate( this, new_opts )
      if nargin <= 1
        new_opts = struct;
      end
      theta0 = 0.9999;
      this.hs.ServiceNodes.Capacity = theta0*this.hs.ServiceNodes.Capacity;
      this.hs.Links.Capacity = theta0*this.hs.Links.Capacity;
      if ~isfield(new_opts, 'CostModel') || ~strcmpi(new_opts.CostModel, 'fixcost')
        [profit,cost] = optimalFlowRate@SimpleSliceOptimizer( this, new_opts );
      else
        if this.NumberFlows == 0
          [profit, cost] = this.handle_zero_flow(new_opts);
        else
          [profit, cost] = optimalFlowRate@SimpleSliceOptimizer( this, new_opts );
          if isfield(new_opts, 'CostModel') && strcmpi(new_opts.CostModel, 'fixcost')
            profit = profit - cost;
          end
          if nargout >= 1
            % When output argument specified, we finalize the VNF capacity.
            % After reconfiguration VNF instance capcity has changed.
            v = this.getVNFCapacity;
            this.Variables.v = v(:);
          end
        end
        profit = profit - this.get_reconfig_cost('const');
      end
      this.ServiceNodes.Capacity = this.old_state.node_capacity;
      this.Links.Capacity = this.old_state.link_capacity;
    end
    
  end
  
  methods (Access=protected)
		% |parameters|: include fields x0, As, bs, Aeq, beq, lbs, ub, lb;
		% |options|: include fields fmincon_opt, CostModel, PricingPolicy (if
		%       CostModel='fixcost'), num_orig_vars;
		function [x, fval] = optimize(this, params, options)
			if isfield(options, 'CostModel') && strcmpi(options.CostModel, 'fixcost')
				[xs, fval, exitflag, output] = ...
					fmincon(@(x)DynamicSlice.fcnSocialWelfare(x, this.hs, ...
					getstructfields(options, 'CostModel')), ...
					params.x0, params.As, params.bs, params.Aeq, params.beq, ...
					params.lb, params.ub, [], options.fmincon_opt);
				this.interpretExitflag(exitflag, output.message);
				if isfield(options, 'bCompact') && options.bCompact
					x = zeros(options.num_orig_vars, 1);
					x(this.I_active_variables) = xs;
				else
					x = xs;
				end
				assert(this.checkFeasible(x, ...
					struct('ConstraintTolerance', options.fmincon_opt.ConstraintTolerance)), ...
					'error: infeasible solution.');
			else
				[x, fval] = optimize@SimpleSliceOpitmizer(this, params, rmstructfields(options, 'CostModel'));
			end
		end
    
    %% handle_zero_flow
    % called by fast slice reconfigure procedure.
    function [profit, cost] = handle_zero_flow(this, new_opts)
      cost = this.hs.getCost(new_opts.PricingPolicy);
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
    end
 
    [profit,cost] = fastReconfigure2(this, action, options);
    [profit,cost] = fastReconfigure(this, action, options);

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
    
    %%
    % dynamic adjust the normalizer for L1-approximation of the reconfiguration cost.
    % This method is used when 'ENABLE_DYNAMIC_NORMALIZER' is set to true.
    function postl1normalizer(this)
      [~, cost] = this.hs.get_reconfig_cost('const');
      [~, linear_cost] = this.hs.get_reconfig_cost('linear', true);
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
              %                         beta.(field_names{i}) = this.raw_beta.(field_names{i})(end)/2;
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
    
		%% TODO postProcessing with reconfiguration examining
		% Fist check constraint violation, then post-process the reconfiguration.
		%   flow processing requirement;
		%   link capacity constraint;
		%   VNF capacity constraint;
		% Resolve constraint violation by down-scaling, which might induce more
		% reconfiguration.
		%
		% As a result, we currently ignore the violation, only focus on
		% how to remove unacessary reconfigurations.
		function [b, vars] = postProcessing(this)
			%             % [TODO] May need post processing for VNF capacity constraint;
			%% post process reconfiguration
			postProcessing@SimpleSlice(this);
			options = getstructfields(this.Parent.options, ...
				{'DiffNonzeroTolerance', 'NonzeroTolerance', 'ConstraintTolerance'});
			if isfield(this.temp_vars, 'v')
				tol_zero = this.Parent.options.NonzeroTolerance;
				var_v = this.temp_vars.v;
				var_v(var_v<tol_zero*max(var_v)) = 0;
				this.Variables.v = var_v;
			end
			%% Another simple method
			% directly recover those little changes, which might cause constraint violation
			tol_vec = options.DiffNonzeroTolerance;
			if ~this.b_derive_vnf
				s = this.diffstate(true);
				if length(s.diff_x) > length(this.Variables.x)
					s.diff_x_norm = s.diff_x_norm(~this.changed_index.x);
					s.diff_z_norm = s.diff_z_norm(~this.changed_index.z);
					old_x = this.old_variables.x(~this.changed_index.x);
					old_z = this.old_variables.z(~this.changed_index.z);
				else
					old_x = this.old_variables.x;
					old_z = this.old_variables.z;
				end
				b_diss_x = s.diff_x_norm < tol_vec;
				this.Variables.x(b_diss_x) = old_x(b_diss_x);
				b_diss_z =  s.diff_z_norm < tol_vec;
				this.Variables.z(b_diss_z) = old_z(b_diss_z);
				if isfield(this.temp_vars, 'v')
					if length(s.diff_v) > length(var_v)
						s.diff_v_norm = s.diff_v_norm(~this.changed_index.v);
						old_v = this.old_variables.v(~this.changed_index.v);
					else
						old_v = this.old_variables.v;
					end
					b_diss_v =  s.diff_v_norm < tol_vec;
					var_v(b_diss_v) = old_v(b_diss_v);
					this.Variables.v = var_v;
				end
				this.diff_state = struct([]);
				b = true; vars = this.Variables;
				return;
				%{
                find(this.As_res * [this.Variables.x; this.Variables.z]>0,1)
                res = this.I_edge_path*this.Variables.x-this.Links.Capacity;
                find(res>0)
                disp(res)
                this.Hdiag*this.Variables.z-this.Variables.v>0
				%}
			end
			if ~this.b_derive_vnf
				%% do processing to discard minor changes.
				% When performing slice dimensioning, this only be performed after the
				% final iteration. The intermediate iteration only use the linear
				% approximated cost.
				s = this.diffstate(true);
				if length(s.diff_x) > length(this.Variables.x)
					s.diff_x = s.diff_x(~this.changed_index.x);
					s.diff_x_norm = s.diff_x_norm(~this.changed_index.x);
					s.diff_z = s.diff_z(~this.changed_index.z);
					s.diff_z_norm = s.diff_z_norm(~this.changed_index.z);
					old_x = this.old_variables.x(~this.changed_index.x);
					old_z = this.old_variables.z(~this.changed_index.z);
				else
					old_x = this.old_variables.x;
					old_z = this.old_variables.z;
				end
				restore_x = this.Variables.x;
				restore_z = this.Variables.z;
				if isfield(this.temp_vars, 'v')
					restore_v = this.temp_vars.v;
					if length(s.diff_v) > length(restore_v)
						s.diff_v = s.diff_v(~this.changed_index.v);
						s.diff_v_norm = s.diff_v_norm(~this.changed_index.v);
						s.mid_v = s.mid_v(~this.changed_index.v);
						old_v = this.old_variables.v(~this.changed_index.v);
					else
						old_v = this.old_variables.v;
					end
					b_diss_v = (s.diff_v<0) & (s.diff_v_norm<tol_vec);
					restore_v(b_diss_v) = old_v(b_diss_v);
					s.diff_v_norm(b_diss_v) = 0;
					s.diff_v(b_diss_v) = 0;
				else
					restore_v = this.Variables.v;
				end
				%                 b_diss_z = (s.diff_z_norm>0) && (s.diff_z_norm<tol_vec) && s.diff_z>0;
				%                 if length(this.old_variables.x) <= length(this.temp_vars.x)
				%                     this.Variables.x(b_diss_x) = this.old_variables.x(b_diss_x);
				%                 else
				%
				%                 end
				% separately process the two part of flows, firstly process the path's
				% that might release some resource could improve the possibility to accept
				% the second part of flows, which need a little more resources.
				%%
				% first case: |x| does not change or the increase amount less than
				% |tol_vec|, try to restore these changes.
				b_diss_x1 = s.diff_x>=0;   % x might not change while z might still change
				Ndc = this.NumberServiceNodes;
				Np = this.NumberPaths;
				Nvnf = this.NumberVNFs;
				af = this.Parent.VNFTable{this.VNFList, 'ProcessEfficiency'};
				%% 1-1
				for p = (find(b_diss_x1))'
					if s.diff_x_norm(p)<tol_vec
						% If x has significant change, usually z also has significant
						% change. But we cannot claim that all z components has
						% significant changes.
						restore_x(p) = old_x(p);
						s.diff_x_norm(p) = 0;
						s.diff_x(p) = 0;
					end
					check_z(1);
				end
				%% 1-2
				delta_vnf = restore_v - this.Hdiag*restore_z;  % residual VNF capacity
				for p = (find(b_diss_x1))'
					check_z(2);
				end
				%% 2-1
				% second case: |x|'s decrease amount less than |tol_vec|, try to restore
				% these changes.
				b_diss_x2 = s.diff_x<0;
				% Links.Capacity has been updated in <finalize>.
				if isfield(this.temp_vars, 'c')
					link_capacity = this.temp_vars.c;
				else
					link_capacity = this.Links.Capacity;
				end
				delta_link = link_capacity - this.I_edge_path*restore_x; % residual link capacity
				for p = (find(b_diss_x2))'
					if s.diff_x_norm(p)<tol_vec
						temp_delta_link = delta_link+this.I_edge_path(:,p).*s.diff_x(p);
						if isempty(find(temp_delta_link<0,1)) % check link capacity constraint
							temp_x = restore_x(p);
							restore_x(p) = old_x(p); % restore: up-scaling
							tf = check_z(1);
							if tf == false % up-scaling restore failed, recover current value.
								restore_x(p) = temp_x;
							else
								delta_link = temp_delta_link;
								s.diff_x_norm(p) = 0;
								s.diff_x(p) = 0;
							end
						end
					else
						check_z(1);
					end
				end
				%% 2-2
				for p = (find(b_diss_x2))'
					if s.diff_x_norm(p)<tol_vec
						temp_delta_link = delta_link+this.I_edge_path(:,p).*s.diff_x(p);
						if isempty(find(temp_delta_link<0,1)) % check link capacity constraint
							temp_x = restore_x(p);
							restore_x(p) = old_x(p); % restore: up-scaling
							tf = check_z(2);
							if tf == false % up-scaling restore failed, recover current value.
								restore_x(p) = temp_x;
							else
								delta_link = temp_delta_link;
								s.diff_x_norm(p) = 0;
								s.diff_x(p) = 0;
							end
						end
					else
						check_z(2);
					end
				end
				%% 1-3
				for p = (find(b_diss_x1))'
					check_z(3);
				end
				%% 2-3
				for p = (find(b_diss_x2))'
					if s.diff_x_norm(p)<tol_vec
						temp_delta_link = delta_link+this.I_edge_path(:,p).*s.diff_x(p);
						if isempty(find(temp_delta_link<0,1)) % check link capacity constraint
							temp_x = restore_x(p);
							restore_x(p) = old_x(p); % restore: up-scaling
							tf = check_z(3);
							if tf == false % up-scaling restore failed, recover current value.
								restore_x(p) = temp_x;
							else
								delta_link = temp_delta_link;
								s.diff_x_norm(p) = 0;
								s.diff_x(p) = 0;
							end
						end
					else
						check_z(3);
					end
				end
				%% reover VNF capacity
				% must be performed after the flow assignement hast been done. Otherwise,
				% there is no chance to restore VNF capacity. Due to the reconfiguration
				% cost constraint, there is no redundant reconfiguration.
				% Instead, after the flow has been restored, there might be redundancy of
				% VNF capacity.
				if isfield(this.temp_vars, 'v')
					b_diss_v = (s.diff_v>0) & (s.diff_v_norm<tol_vec);
					idv = this.Hdiag*restore_z<=old_v & b_diss_v; % if the old value satisfy the capacity constraint
					restore_v(idv) = old_v(idv);
					s.diff_v_norm(idv) = 0;
					s.diff_v(idv) = 0;
					
					%                     delta_dc = full(sum(reshape(this.Variables.v, Ndc, Nvnf),2)) -...
					%                         full(sum(reshape(restore_v, Ndc, Nvnf),2)); % residual node capacity
					%                     % all nodes can be check in parallel, with VNF varying
					%                     idv = 1:Ndc;
					%                     for v = 1:Nvnf
					%                         tv = restore_v(idv);
					%                         old_tv = old_v(idv);
					%                         idtv = (s.diff_v_norm(idv)>0) & (s.diff_v_norm(idv)<tol_vec) & s.diff_v(idv)<0;
					%                         tv(idtv) = old_tv(idtv);
					%                         temp_delta_dc = delta_dc + restore_v(idv) - tv;
					%                         fx = temp_delta_dc >= 0;
					%                         restore_v(idv(fx)) = tv(idv(fx));
					%                         delta_dc(fx) = temp_delta_dc(fx);
					%                         s.diff_v_norm(idv(fx)) = 0;
					%                         s.diff_v(idv(fx)) = 0;
					%                         idv = idv + Ndc;
					%                     end
					this.Variables.v = restore_v;
				end
				this.Variables.x = restore_x;
				this.Variables.z = restore_z;
				this.diff_state = struct([]);
			end
			b = true; vars = this.Variables;
			
			%% nest function: check z
			function tf = check_z(t)
				tf = true;
				b_reset = false(Nvnf,1);
				temp_z = zeros(Ndc, Nvnf);
				temp_sz = zeros(Ndc, Nvnf);
				temp_szn = zeros(Ndc, Nvnf);
				delta_sum = zeros(Nvnf,1);
				if t > 1
					temp_delta_vnf = zeros(Ndc, Nvnf);
				end
				
				idz = (1:Ndc) + (p-1)*Ndc;
				idn = 1:Ndc;
				for k = 1:Nvnf
					tz = restore_z(idz);
					temp_z(:,k) = tz;
					temp_sz(:,k) = s.diff_z(idz);
					temp_szn(:,k) = s.diff_z_norm(idz);
					if t > 1
						temp_delta_vnf(:,k) = delta_vnf(idn);
					end
					old_tz = old_z(idz);
					switch t
						case 1
							% try restore z(:,p,f) that should be restored by down-scaling;
							% make space for other requests.
							idtz = (s.diff_z(idz)>0) & (s.diff_z_norm(idz)<tol_vec);
						case 2
							%%% try restore z(:,p,f) that includes both up-sclaing and down-scaling;
							idtz = (s.diff_z_norm(idz)>0) & (s.diff_z_norm(idz)<tol_vec);
						case 3
							%%% try restore z(:,p,f) that only includes up-scaling;
							idtz = (s.diff_z(idz)<0) & (s.diff_z_norm(idz)<tol_vec);
					end
					if isempty(find(idtz,1))
						if nargout == 1 && dot(this.I_dc_path(:,p), tz) < af(k)*restore_x(p)
							% when x increase, but z keeps unchange, so we need to check the processing
							% constraint.
							tf = false;
						end
						idz = idz + Np*Ndc;
						idn = idn + Ndc;
						continue;
					end
					tz(idtz) = old_tz(idtz);  % recover: down-scaling | down/up-scaling | up-scaling
					if t == 2 || t == 3
						delta_tz = restore_z(idz) - tz; % current - past
						t_delta_vnf = delta_vnf(idn)+this.I_dc_path(:,p).*delta_tz;
					end
					if t == 1
						if dot(this.I_dc_path(:,p), tz) >= af(k)*restore_x(p) % check processing constraint
							restore_z(idz) = tz;
							s.diff_z_norm(idz(idtz)) = 0;
							s.diff_z(idz(idtz)) = 0;
							delta_sum(k) = 0;
						else
							b_reset(k) = true;
						end
					end
					if t == 2
						if dot(this.I_dc_path(:,p), tz) >= af(k)*restore_x(p) &&...
								isempty(find(t_delta_vnf<0,1))
							%% check both processing constraint and VNF capacity constraint
							restore_z(idz) = tz; % accept old value or keep current value.
							s.diff_z_norm(idz(idtz)) = 0;
							s.diff_z(idz(idtz)) = 0;
							delta_vnf(idn) = t_delta_vnf; % update residual capacity
							delta_sum(k) = sum(delta_tz);
						else
							b_reset(k) = true;
						end
					end
					if t == 3
						if (nargout == 0 && isempty(find(t_delta_vnf<0,1))) || (nargout == 1 && isempty(find(t_delta_vnf<0,1)) && dot(this.I_dc_path(:,p), tz) >= af(k)*restore_x(p))
							restore_z(idz) = tz; % accept old value or keep current value.
							s.diff_z_norm(idz(idtz)) = 0;
							s.diff_z(idz(idtz)) = 0;
							delta_vnf(idn) = t_delta_vnf; % update residual capacity
							delta_sum(k) = sum(delta_tz);
						else
							b_reset(k) = true;
						end
					end
					idz = idz + Np*Ndc;
					idn = idn + Ndc;
				end
				if ~isempty(find(b_reset,1))
					tf = false;
					idz = (1:Ndc) + (p-1)*Ndc;
					idn = 1:Ndc;
					for k = 1:Nvnf
						if ~b_reset(k) && delta_sum(k) < 0
							restore_z(idz) = temp_z(:,k);
							s.diff_z(idz) = temp_sz(:,k);
							s.diff_z_norm(idz) = temp_szn(:,k);
							if t > 1
								delta_vnf(idn) = temp_delta_vnf(:,k);
							end
						end
						idz = idz + Np*Ndc;
						idn = idn + Ndc;
					end
				end
			end
    end
		
    %% Convert the optimizatio results to temporary variables
    % temp variables might include: x,z,v,tx,tz,tv,c (see also
    % <Dynamic.priceOptimalFlowRate>);
    function convert(this, x, ~)
      Np = this.hs.NumberPaths;
      this.temp_vars.x = x(1:Np);
      this.temp_vars.z = x((Np+1):this.NumberVariables);
      offset = this.NumberVariables;
      this.temp_vars.v = x(offset+(1:this.num_varv));
      offset = offset + this.num_varv;
      if this.invoke_method >= 2
        this.temp_vars.tx = x(offset+(1:Np));
        this.temp_vars.tz = x(offset+((Np+1):this.NumberVariables));
        offset = offset + this.NumberVariables;
        this.temp_vars.tv = x(offset+(1:this.num_varv));
        offset = offset + this.num_varv;
      end
      this.temp_vars.c = x(offset+(1:this.NumberLinks));
      %             if nargin >= 3
      %                 this.Variables = getstructfields(this.temp_vars, {'x','z','v','c'}, 'ignore');
      %             end
		end
    
		%% TODO
		% optimize only on the subgraph, which carries the flows that intersection with
		% the arriving/departuring flow. this can significantly reduce the problem scale
		% when there is lots of flow, and the numbe of hops of the flow is small.
		% 		function temp_vars = get_temp_variables(this, bfull)
		% 			if nargin >= 2 && bfull
		% 				temp_vars = [this.temp_vars.x; this.temp_vars.z; this.temp_vars.v; ...
		% 					this.temp_vars.tx; this.temp_vars.tz; this.temp_vars.tv];
		% 			else
		% 				temp_vars = [this.temp_vars.x; this.temp_vars.z; this.temp_vars.v];
		% 			end
		% 		end
		
    %% Find the index of variables for the newly added/removed flow
    % assume the newly added/removed flow index is u, the the index of
    % corrsponding flow and VNF allocation variables is
    % [x_u, z_npf], for all p in p_u.
    %
    % NOTE: This function should be called after adding new flow, or before removing
    % departing flow.
    function identify_change(this, changed_path_index)
      global DEBUG; %#ok<NUSED>
      Nsn = this.hs.NumberServiceNodes;
      Nvf = this.hs.NumberVNFs;
      this.changed_index.x = logical(sparse(changed_path_index));
      base_z_index = false(Nsn, this.hs.NumberPaths);
      base_z_index(:,changed_path_index) = true;
      if ~isempty(fieldnames(this.hs.net_changes))
        base_z_index(this.hs.net_changes.DCIndex,:) = true;
      end
      base_z_index = sparse(base_z_index);
      this.changed_index.z = repmat(base_z_index(:), Nvf, 1);
      if ~isempty(fieldnames(this.hs.net_changes))
        this.changed_index.v = sparse(Nsn, Nvf);
        this.changed_index.v(this.hs.net_changes.DCIndex,:) = true;
        this.changed_index.v = logical(this.changed_index.v(:));
      end
    end
    
    % update reconfiguration cost with scaler.
    function update_reconfig_costinfo(this, action, bDim)
      %% Period re-dimensioing
      % The number of virtual nodes/links is not change in the slice, as well as the
      % number of VNF instances.
      % When performing re-dimensioning, the reconfiguration cost is larger than
      % that of fast reconfiguration, so we need to update the reconfiguration cost.
      if nargin >= 3 && bDim
        this.Parent.updateRedimensionCost(this);
        this.update_reconfig_costvinfo();
      end
      if strcmpi(action, 'add')
        this.x_reconfig_cost = (this.I_edge_path)' * this.Links.ReconfigCost;
        this.z_reconfig_cost = repmat(this.ServiceNodes.ReconfigCost, ...
          this.NumberPaths*this.NumberVNFs, 1);
        %% Ignore the reconfiguration cost when adding flow
        % The newly added flow's reconfiguration cost set to zero.
        % This helps the new flow be admitted into the slice, when there is
        % resource while the reconfiguration cost of the new flow is too high.
        % This can be explained as: the reconfiguration cost of the first time can
        % be treated as that it is distrbuted in the whole life time of the flow.
        % As long as the flow's lifetime is long enough, the first time cost could
        % be ignored.
        %
        % In the output, we still count the reconfiguration cost for the new flow,
        % see also <diffstate> and <get_reconfig_cost>.
        this.x_reconfig_cost(this.changed_index.x) = 0;
        this.z_reconfig_cost(this.changed_index.z) = 0;
        if this.ENABLE_DYNAMIC_NORMALIZER
          beta =this.getbeta();
          this.topts.x_reconfig_cost = beta.x*this.x_reconfig_cost;
          this.topts.z_reconfig_cost = beta.z*this.z_reconfig_cost;
        else
          this.topts.x_reconfig_cost = this.options.ReconfigScaler*this.x_reconfig_cost;
          this.topts.z_reconfig_cost = this.options.ReconfigScaler*this.z_reconfig_cost;
        end
      else
        %%
        % Since we need maintain information of the deleted variables, we should use
        % the old state to calculate the information.
        this.x_reconfig_cost = ...
          (this.old_state.I_edge_path)' * this.Links.ReconfigCost;
        old_num_paths = length(this.old_state.path_owner);
        this.z_reconfig_cost = repmat(this.ServiceNodes.ReconfigCost, ...
          old_num_paths*this.NumberVNFs, 1);
        %%
        % The removed flow's reconfiguration cost is not considered in the
        % optimization, which is constant.
        if this.ENABLE_DYNAMIC_NORMALIZER
          beta =this.getbeta();
          this.topts.x_reconfig_cost = ...
            beta.x*this.x_reconfig_cost(~this.changed_index.x);
          this.topts.z_reconfig_cost = ...
            beta.z*this.z_reconfig_cost(~this.changed_index.z);
        else
          this.topts.x_reconfig_cost = ...
            this.options.ReconfigScaler*this.x_reconfig_cost(~this.changed_index.x);
          this.topts.z_reconfig_cost = ...
            this.options.ReconfigScaler*this.z_reconfig_cost(~this.changed_index.z);
        end
      end
    end
    
    % update reconfiguration cost with scaler.
    % If dimensioning is applied, the orginal cost is firstly updated.
    function update_reconfig_costvinfo(this)
      this.vnf_reconfig_cost = this.hs.Parent.options.VNFReconfigCoefficient * ...
        repmat(this.hs.ServiceNodes.ReconfigCost, this.hs.NumberVNFs, 1);
      if this.ENABLE_DYNAMIC_NORMALIZER
        beta = this.getbeta();
        this.topts.vnf_reconfig_cost = beta.v*this.vnf_reconfig_cost;
      else
        this.topts.vnf_reconfig_cost = ...
          this.options.ReconfigScaler*this.vnf_reconfig_cost;
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
				ds.tI_dc_path = this.I_dc_path;
				% adding flow, |path| will increase, and |edge| might(might not) increase;
				ds.tI_edge_path = this.I_edge_path;
			else
				new_x = zeros(size(this.old_variables.x));
				if isfinal
					new_x(~this.changed_index.x) = this.Variables.x;
				else
					new_x(~this.changed_index.x) = this.temp_vars.x;
				end
				%% The removed flow has been considered with the reconfguration cost
				% The removed flow's reconfiguration cost is constant and does not appear
				% in the optimization, see also <update_reconfig_costinfo>.
				%
				% removing flow: |path| will decrease, and |edge| might(might not) decrease;
				ds.tI_dc_path = this.old_state.I_dc_path;
				ds.tI_edge_path = this.old_state.I_edge_path;
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
				new_z = zeros(size(this.old_variables.z));
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
						new_vnf_capacity = this.VNFCapacity(:);
					else
						new_vnf_capacity = this.temp_vars.v;
					end
				else
					new_vnf_capacity = zeros(length(this.old_variables.v),1);
					if isfinal
						new_vnf_capacity(~this.changed_index.v) = this.VNFCapacity(:);
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

  end
end