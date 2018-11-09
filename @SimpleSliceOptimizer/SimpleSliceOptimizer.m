classdef SimpleSliceOptimizer < SliceOptimizer
  properties
		As_res;   % Coefficient matrix of the processing constraints.
		Hrep;     % Coefficient used to compute node resource consumption.
		Hdiag;		% Coefficient used to compute VNF intance capacity.
    num_varz; % number of variables in vector |z_npf|.
  end

  properties (Dependent)
		VNFCapacity;
	end
	
	methods
		profit = getProfit(this, options);
	end
	methods(Static)
		[profit, grad] = fcnSocialWelfare(x_vars, op, options);
		%% fcnProfit
    % Evalute the objective function and gradient.
    %
    %      [profit, grad] = fcnProfit(vars, op, options)
    %
    % |grad|: the gradient value of the objective function.
    % The upper bound number of non-zero elements in the gradient vector: the gradient on path
    % variable is nonzeros, so there is |P| components; whether the gradient on node variable
    % is zeros is decided by the node-path incidence matrix, i.e. |nnz(I_dc_path)*F|.
    [profit, grad] = fcnProfit(vars, op, options);
		%% fcnHessian
    % Hessian matrix of the Largrangian.
    %
    %      hs = fcnHessian(var_x, ~, slice, options)
    %
    % Since the problem only contains linear constraint, the hessian matrix of the
    % Largrangian is equal to the second derivatives of the objective function, and the
    % Largrangian multipliers $\lambda$ takes no effect.
    % The Hessian matrix contains only $P^2$ nonzeros elements on the diagonal,
    % which is the second derviatives on path variables.
    hs = fcnHessian(var_x, ~, slice, options);
	end
  
  methods
    function this = SimpleSliceOptimizer(slice, options)
			if nargin >= 2
				args = {slice, options};
			elseif nargin == 1
				args = {slice};
			else
				args = {};
			end
      this@SliceOptimizer(args{:});

      this.getAs_res;
      this.getHrep;
      this.getHdiag;
    end
  end
  
  %% Property Get Methods
  methods
    function c = get.VNFCapacity(this)
      if ~isfield(this.Variables, 'v') || isempty(this.Variables.v)
        warning('VNF capacity not set, set to VNF load.');
        c = this.getVNFCapacity;
        this.Variables.v = c;
      else
        c = this.Variables.v;
      end
    end
    
		function n = get.num_varz(this)
			n = this.hs.NumberVNFs*this.hs.NumberServiceNodes*this.hs.NumberPaths;
		end
    
  end
  
  %% Public Methods
  methods
		[profit,cost] = optimalFlowRate(this, new_opts);		% see also <SliceOptimizer>
		[utility, load] = priceOptimalFlowRate(this, x0, options); % see also <SliceOptimizer>

		
    function clear(this)
      this.Variables.x = [];
      this.Variables.z = [];
      this.temp_vars.x = [];
      this.temp_vars.z = [];
		end
    
		%%
		% Implement <SliceOptimizer.getFlowRate>
    function r = getFlowRate(this, vars)
      if nargin == 1
        r = this.I_flow_path * this.Variables.x;
			elseif isnumeric(vars)
        r = this.I_flow_path * vars;
			else
				r = this.I_flow_path * vars.x;
      end
      r = full(r);
		end
		
		%%
		% Override <SliceOptimizer.setPathBandiwdth>
		function setPathBandwidth(this, x)
			if nargin == 1
				x = this.Variables.x;
			end
			p = 1;
			for i = 1:this.hs.NumberFlows
				pathlist = this.hs.FlowTable{i,'Paths'};
				for l = 1:length(pathlist)
					pathlist{l}.bandwidth = x(p);
					p = p + 1;
				end
			end
		end
		
    function b = checkFeasible(this, vars, opt_opts)
      if nargin <= 1 || isempty(vars)
        vars = [this.Variables.x; this.Variables.z];
      else
        vars = vars(1:this.NumberVariables);
      end
      if nargin >=3 && isfield(opt_opts, 'ConstraintTolerance')
        b = isempty(find(this.As_res*vars>opt_opts.ConstraintTolerance,1));
      else
        b = isempty(find(this.As_res*vars>1e-10,1));
      end
    end
    
    function initializeProblem(this)
      
    end
    
		%     function setProblem(this, varargin)
		% 			idx = true(length(varargin,1));  % filter varargins to
		% 			setProblem@SliceOptimizer
		% 			for i = 1:2:(length(varargin)-1)
		% 				switch varargin{i}
		% 					case ''
		% 					otherwise
		% 						error('un-recoginzed input for the problem.');
		% 				end
		% 			end
		% 		end
		
    %% Get Load
    % <getLinkLoad> and <getNodeLoad> are only used inside the <Slice> class. The
    % substrate network cares how much resources are ocuppied by the slice, that is
    % the capacity of resources. See also <getLinkCapacity> and <getNodeCpacity>.
    function ye = getLinkLoad(this, isfinal, vars)
      % retrive the link load of the slice, given the path variables.
			if nargin >= 3
				if isnumeric(vars)
					edge_vars = vars;
				else
					edge_vars = vars.x;
				end
			else
				if nargin == 1 || isfinal
					edge_vars = this.Variables.x;
				elseif nargin >= 2
					edge_vars = this.temp_vars.x;
				end
			end
			ye = this.I_edge_path * edge_vars;
      ye = full(ye);
		end
		
		%%%
		% retrive the node load of the slice, given the node variables. The node variables
		% represent the resource allocation of data center nodes.
		%
		% |v_n|: data center's resource consumption.
		function v_n = getNodeLoad(this, isfinal, vars)
			if nargin >= 3
				if isnumeric(vars)
					node_vars = vars;
				else
					node_vars = vars.z;
				end
			else
				if nargin == 1 || isfinal
					node_vars = this.Variables.z;
				elseif nargin >= 2
					node_vars = this.temp_vars.z;
				end
			end
			%%
			% |node_vars| is index by |(node,path,function)|.
			v_n = full(this.Hrep*node_vars);
			% assert(isempty(find(abs(this.Hrep*node_vars-v_n)>10^-6,1)), 'error: unequal node load');
			
			%% Alternative way to compute the node load.
			% _Hrep_ replace the following procedure.
			%             node_load = sum(f, node_vars(:,:,f).*I_dc_path).
			%             Nsn = this.op.NumberServiceNodes;
			%             Np = this.op.NumberPaths;
			%             v_n = zeros(Nsn,1);
			%             np = Nsn * Np;
			%             z_index = 1:np;
			%             for i = 1:this.NumberVNFs
			%                 node_vars_fi = reshape(node_vars(z_index), Nsn, NP);
			%                 v_n = v_n + sum(this.op.I_dc_path.*node_vars_fi,2);
			%                 z_index = z_index + np;
			%             end
		end
    %% priceOptimalFlowRateCompact
    % Find the optimal flow rate that maximizing the net profit of the network slice.
    %
    %      x = priceOptimalFlowRateCompact(this, x0)
    %
    % The optimizition procedure in this method remove the unnecessary components from the
    % independent variable |x|, so that the problem scale is cut down.
		
		
		function runtime = optimalFlowRateSingleSlice(this, slice_data, options)
			slice = this.hs;
			Nf = slice.NumberFlows;
			if options.SlicingMethod == SlicingMethod.SingleFunction
				this.getAs_res(slice_data.flow_owner, slice_data.Alpha_f);
				options.Alpha_f = slice_data.Alpha_f;
			elseif options.SlicingMethod == SlicingMethod.SingleNormal
				%% Coefficient for global optimization
				% When all slices are combined into one slice, a VNF might not be used by all paths
				% (_i.e._ all flows). If a VNF |f| is not used by a path |p|, there is no
				% processing-rate constraints on $f \times p$(|delete_items|). To form the constraint
				% coefficient matrix, the related items in |As| should be removed. See also <Slice
				% file://E:/workspace/MATLAB/Projects/Documents/CloudNetworking/Slice.html>.
				%
				% By the way, $z_{n,p,f}=0, \forall n$, if |p| does not use NFV |f|.
				I_flow_function = zeros(Nf, slice.NumberVNFs);
				for f = 1:Nf
					[~, vid] = ismember(slice.Parent.slices{slice_data.flow_owner(f)}.VNFList, slice.VNFList);
					I_flow_function(f, vid) = 1;
				end
				I_path_function = this.I_flow_path'*I_flow_function;
				this.As_res = this.As_res(logical(I_path_function(:)),:);
			end
			
			if nargout == 1
				tic;
			end
			% Only return intermediate results, so no return value provided.
			this.optimalFlowRate(options);
			if nargout == 1
				runtime.Serial = toc;
				runtime.Parallel = runtime.Serial;
			end
			if options.SlicingMethod == SlicingMethod.SingleNormal
				nz = slice.NumberServiceNodes*slice.NumberPaths;
				z_index = 1:nz;
				for v = 1:slice.NumberVNFs
					mask_npf = this.I_dc_path.*I_path_function(:,v)'; % compatible arithmetic operation
					this.temp_vars.z(z_index) = mask_npf(:).*this.temp_vars.z(z_index);
					z_index = z_index + nz;
				end
			end
			
			%% Partition the network resources according to the global optimization
			pid_offset = 0;
			z_npf = reshape(full(this.temp_vars.z), slice.NumberServiceNodes, slice.NumberPaths, slice.NumberVNFs);
			% node_load = zeros(this.NumberNodes, 1);
			% link_load = zeros(this.NumberLinks, 1);
			fmincon_opt = optimoptions('fmincon');
			for s = 1:slice.Parent.NumberSlices
				sl = slice.Parent.slices{s};
				op = sl.Optimizer;
				pid = 1:sl.NumberPaths;
				op.temp_vars.x = this.temp_vars.x(pid_offset+pid);
				nid = sl.getDCPI;       % here is the DC index, not the node index.
				[~, vid] = ismember(sl.VNFList, slice.VNFList);
				op.temp_vars.z = ...
					reshape(z_npf(nid,pid+pid_offset,vid), op.NumberVariables-sl.NumberPaths, 1);
				pid_offset = pid_offset + sl.NumberPaths;
				assert(op.checkFeasible([op.temp_vars.x; op.temp_vars.z], ...
					struct('ConstraintTolerance', fmincon_opt.ConstraintTolerance)), ...
					'error: infeasible solution.');
				sl.ServiceNodes.Capacity = sl.getNodeCapacity(false);
				sl.Links.Capacity = sl.getLinkCapacity(false);
				% DEBUG
				%     eid = sl.Links.PhysicalLink;
				%     node_load(nid) = node_load(nid) + sl.Nodes.Capacity;
				%     link_load(eid) = link_load(eid) + sl.Links.Capacity;
			end
			% disp(max(node_load-this.readNode('Capacity')));
			% disp(max(link_load-this.readLink('Capacity')));
		end
		
		%% Post Processing
		% *Flow processing constraint* might be violated, due to rounding of x,z. Three
		% candidate methods can be applied for post processing:
		%       (1) round the x components, so that flow processing demand is meet;
		%       (2) drop the items in x that violates the constraint;
		%       (3) recover those z components corresponding to violated constraints. This
		%           would result in more small components.
		% NOTE 1: (2) and (3) cannot retain the feasibility of the solution, since it
		% ignore small violations.
		%
		% NOTE 2: due to rounding to zero, the capacity constraints will not be violated
		% (assuming the original solution not violate the capacity constraints).
		%
		% NOTE 3: the optimization variables is continuous, discrete variables (due to the
		% scheduling granularity) are not considered here.
		function [tf, vars] = postProcessing(this)
			global DEBUG INFO;
			Np = this.hs.NumberPaths;
			Nvf = this.hs.NumberVNFs;
			Nsn = this.hs.NumberServiceNodes;
			var_x = this.temp_vars.x;
			var_z = this.temp_vars.z;
			tol_zero = this.options.NonzeroTolerance;
			var_x(var_x<tol_zero*max(var_x)) = 0;
			var_z(var_z<tol_zero*max(var_z)) = 0;       % [min(var_z(var_z~=0)) max(var_z)]
			A1 = this.As_res(:,1:Np);
			A2 = this.As_res(:,(Np+1):end);
			v1 = A1*var_x;
			v2 = A2*var_z;
			index_violate = find(v1+v2>0);      % re = v1+v2
			if ~isempty(index_violate)
				message = sprintf('[%s] Maximal violation is %.4f before processing.', ...
					calledby, full(max(v1+v2)));
				if ~isempty(DEBUG) && DEBUG
					warning(message); %#ok<SPWRN>
				elseif ~isempty(INFO) && INFO
					cprintf('SystemCommands', 'Warning: %s\n', message);
				end
				
				post_process = this.options.PostProcessing; % TODO: replace with the optimizer.options.
				switch post_process
					case 'round'
						b_violate = false(size(v1));
						b_violate(index_violate) = true;
						pid_offset = 0:Np:((Nvf-1)*Np);
						for i = 1:Np
							pid = pid_offset + i;
							p_violate = find(b_violate(pid));
							if ~isempty(p_violate)
								t = v1(pid(p_violate))./abs(v2(pid(p_violate)));
								var_x(i) = var_x(i) ./ max(t);
							end
						end
						% The rounding error may still leads to postive error, therefore
						% we set a relatively small tolerance, i.e. 1e-10.
						assert(this.checkFeasible([var_x; var_z], ...
							struct('ConstraintTolerance', 1e-10)), ...
							'[%s] failed infeasible solution.', calledby);
					case {'drop','recover'}
						%% Ignore the tiny components
						% Use the differentiation of |v1-v2|/|v1+v2| is more flexible than use only
						% |v1-v2|.
						%   # if |v1-v2| is remarkable, both form can eaily identify the difference;
						%   # If both v1 and v2 is large while |v1-v2| is small, the former can
						%     identify the small value regardless of the tolerance, while the later
						%     relies on the tolerance;
						%   # Otherwise, if both v1 and v2 is small (most small value should have been
						%     rounding in the former process), then both form rely on tolerance.
						%
						%% ISSUE
						% |tol_zero| should not be too large compared to the absolute
						% value of variables. It leads to inaccurate solution (the
						% violation still exists), which cannot pass the feasible assertion.
						% Therfore, *'drop' is not recommend*.  So |tol_zero| should be
						% smaller, so that all violated ones can be
						% identified.
						%   # If we use the 'drop' option, some large components of the
						%     solution might be discarded (even if the absolute value of
						%     residual is large),  which degenerate the solution to the
						%     sub-optmal.
						%   # If we use the 'recover' option, then some trivial components
						%     will be recovered. However the quality of solution declines.
						%
						% Since the |tol_zero| might lead to miss some violated components, we
						% do not perform feasible assertion after processing. We Since the
						% constraint tolerance in optimization is 10^-6 by default, the
						% residual error's magnitude may be higher than 10^-6.
						% Therefore, we have to set |tol_zero| to 10^-3~10^-4 (tuned according
						% to the real data).
						re = v1 + v2;
						re(re<tol_zero) = 0;
						mid_v = (abs(v1)+abs(v2))/2;
						nz_index = find(mid_v~=0);
						tol_con = this.options.ConstraintTolerance; % TODO: replace with the optimizer.options.
						b_violate = (re(nz_index)./mean(mid_v(nz_index))) > tol_con;
						% => mean(mid_v(nz_index)) | max(mid_v(nz_index))
						vi_idx = nz_index(b_violate);
						pidx = mod(vi_idx-1, Np)+1;
						if strcmpi(post_process, 'recover')
							fidx = ceil(vi_idx/Nvf);
							for i = 1:length(vi_idx)
								zidx = (1:Nsn)+ ((fidx(i)-1)*Nsn*Np + (pidx(i)-1)*Nsn);
								var_z(zidx) = this.temp_vars.z(zidx);  % recover
							end
						else
							pidx = unique(pidx);
							var_x(pidx) = 0;
						end
					otherwise
						error('error: [%s] invalid processing option (%s).', calledby, post_process);
				end
			end
			if nargout >= 2
				vars = [var_x; var_z];
			end
			this.Variables.x = var_x;
			this.Variables.z = var_z;
			tf = true;
		end
		
	end
	methods (Static)
		function [profit, grad] = fcnSocialWelfareCompact(act_vars, op)
			vars = zeros(options.num_orig_vars,1);
			vars(op.I_active_variables) = act_vars;
			
			if nargout <= 1
				profit = SimpleSliceOptimizer.fcnSocialWelfare(vars, op);
			else
				[profit, grad] = SimpleSliceOptimizer.fcnSocialWelfare(vars, op);
				grad = grad(op.I_active_variables);
			end
		end
		
		function [profit, grad] = fcnProfitCompact(act_vars, slice, options)
			vars = zeros(options.num_orig_vars,1);
			vars(slice.I_active_variables) = act_vars;
			
			% we extend the active variables by adding zeros to the inactive ones.
			if nargout <= 1
				profit = SimpleSlice.fcnProfit(vars, slice, options);
			else
				[profit, grad] = SimpleSlice.fcnProfit(vars, slice, options);
				% eliminate the inactive variable's derivatives.
				grad = grad(slice.I_active_variables);
			end
		end
	end

  
  %% Protected Methods
  methods (Access = protected)
		function n = get_number_variables(this)
			n = (this.hs.NumberVNFs*this.hs.NumberServiceNodes+1)*this.hs.NumberPaths;
    end
    
    function n = get_number_linear_constraints(this)
			n = size(this.As_res,1);
    end

    %     function temp_vars = get_temp_variables(this)
    %       temp_vars = [this.temp_vars.x; this.temp_vars.z];
    %     end
				
    %% Linear constraint without consdiering the bound constraint.
    % |As| takes the following form
    %
    % $$\left[ \begin{array}{ccccc}
    %   I_1     & H_s &      &          &       \\
    %   I_2     &     & H_s  &          &       \\
    %   \vdots  &     &      &  \ddots  &       \\
    %   I_F     &     &      &          & H_s
    % \end{array} \right]
    % \left[ \begin{array}{c}x\\z_1\\z_2\\\vdots\\z_f\\\vdots\\z_F\end{array}\right] $$
    %
    % where
    %
    % $$I_f = {\left[ \begin{array}{cccc}
    % \alpha_f &          &        &         \\
    %          & \alpha_f &        &         \\
    %          &          & \ddots &         \\
    %          &          &        & \alpha_f
    % \end{array} \right]}_{P\times P},$$
    % $$H_s = {\left[ \begin{array}{cccccccccc}
    % -h_{1,1} & \cdots & -h_{N,1} &          &        &          &        &         &        & \\
    %          &        &          & -h_{1,2} & \cdots & -h_{N,2} &        &         &        & \\
    %          &        &          &          &        &          & \ddots &         &        & \\
    %          &        &          &          &        &          &        &-h_{1,P} & \cdots & -h_{N,P}
    % \end{array} \right]}_{P\times NP},$$
    % $$ x = \left[\begin{array}{c}x_1\\x_2\\ \vdots\\x_P\end{array}\right]$$
    % $$ z_f = \left[\begin{array}{c}z_{1,1,f}\\ \vdots\\ z_{N,1,f}\\z_{1,2,f}\\
    %   \vdots\\ z_{N,2,f}\\ \vdots\\z_{N,P,F}\end{array}\right]$$
    %
    % According to the martix formulation, the number of non-zero elements in |As|
    % is equal to |F*(P+nnz(Hs))|, where |nnz(Hs)| is equal to the number of
    % nonzero elements in |I_dc_path|.
    %
    % NOTE: it is not necessary to evaluate As_res each time when visiting it. so, we
    % define a normal function to update the property |As_res|.
    function As = getAs_res(this, flow_owner, alpha_f)
			slice = this.hs;
      Nsn = slice.NumberServiceNodes;
      Np = slice.NumberPaths;
      Nvf = slice.NumberVNFs;       % For 'single-function', NV=1.
      nnz_As = Nvf*(Np+nnz(this.I_dc_path));
      num_lcon = Np*Nvf;
      As = spalloc(num_lcon, this.NumberVariables, nnz_As);
      row_index = 1:Np;
      for f = 1:Nvf
        if nargin >= 3   % used when treat all VNFs as one function.
          for p = 1:Np
            As(p,p) = alpha_f(flow_owner(slice.path_owner(p))); %#ok<SPRIX>
          end
        else
          af = slice.Parent.VNFTable{slice.VNFList(f),{'ProcessEfficiency'}};
          As(row_index,1:Np) = af * speye(Np); %#ok<SPRIX>
        end
        row_index = row_index + Np;
      end
      col_index = 1:Nsn;
      Hst = zeros(Np, Nsn*Np);     % |Hst| is a staircase-like diagnoal.
      for p = 1:Np
        Hst(p, col_index) = -this.I_dc_path(:,p);
        col_index = col_index + Nsn;
      end
      As(:, (Np+1):end) = block_diag(Hst, Nvf);
      this.As_res = As;
    end
    
    % equation:
    %    $\sum_{p,f}{b_{np}z_{npf}} \le V_{n}$
    %
    % For the _node capacity constraint_, the coefficient matrix is filled row-by-row.
    % On row i, the non-zero elements located at (i-1)+(1:NC:((NP-1)*NC+1)) for the
    % first |NC*NP| columns, and then the first |NC*NP| columns are duplicated for
    % |NV| times, resulting in |NC*NP*NV| columns.
    %
    % The base matrix is the same as in <Hdiag>.
    function H = getHrep(this)
      Np = this.hs.NumberPaths;
      Nsn = this.hs.NumberServiceNodes;
      Nvf = this.hs.NumberVNFs;
      
      H = spalloc(Nsn, this.num_varz, nnz(this.I_dc_path)*Nvf);
      col_index = (1:Nsn:((Np-1)*Nsn+1))';      % the vnf variable index related to the first node
      col_index = repmat(col_index, 1, Nvf);   % derive other node's vnf variable index
      for v = 2:Nvf
        col_index(:,v) = col_index(:,v-1) + Nsn*Np;
      end
      col_index = col_index(:);
      for n = 1:Nsn    % all path that use node n is counted
        H(n, col_index) = repmat(this.I_dc_path(n,:), 1, Nvf);  %#ok<SPRIX>
        col_index = col_index + 1;
      end
      this.Hrep = H;
    end
    
    % equation:
    %    $\sum_{p}{b_{np}z_{npf}} \le v_{nf}$
    %
    % For the node capacity constraint for VNF |f|, the coeffiect matrix Hs is
    % filled block-by-block. In the |k|th block (NC*NC), we put the |k|th column
    % of H_np to the diagnal of the block. Then |H_np| is duplicated into a larger
    % diagonal corresponding to all VNFs.
    %
    % The capacity of VNF instance |VNFCapacity|, has been calculate after
    % allocating resource to the slice. Some component of |VNFCapacity| might be
    % zero, since VNF f may not be located at node n.
    %
    % See also <Hrep>.
    function Hs = getHdiag(this)
      global DEBUG; %#ok<NUSED>
      Nsn = this.hs.NumberServiceNodes;
      Np = this.hs.NumberPaths;
      Hs_np = spalloc(Nsn, Nsn*Np, nnz(this.I_dc_path));
      col_index = 1:Nsn;
      for p = 1:Np
        Hs_np(1:Nsn,col_index) = diag(this.I_dc_path(:,p));  %#ok<SPRIX>
        col_index = col_index + Nsn;
      end
      Hs = block_diag(Hs_np, this.hs.NumberVNFs);
      this.Hdiag = Hs;
    end
    
    % Called by <optimalFlowRate>, if 'bCompact = true', options should be speicifed with
    % 'num_orig_vars'.
    % See also <SlicingMethod>.
    function [x, fval] = optimize(this, options)
			if options.SlicingMethod.IsSingle
				% 'normal', 'single-function'
				fobj = @SimpleSliceOptimizer.fcnSocialWelfare;
			elseif options.SlicingMethod.IsPricing
				fobj = @SimpleSliceOptimizer.fcnProfit;
			else
				error('error: unrecognized SlicingMethod.');
			end
      
			[xs, fval, exitflag, output] = ...
				fmincon(@(x)fobj(x, this, options), ...
				this.problem.x0, this.problem.As, this.problem.bs, this.problem.Aeq, ...
				this.problem.beq, this.problem.lb, this.problem.ub, [], options.fmincon_opt);
      SimpleSliceOptimizer.interpretExitflag(exitflag, output);
      if isfield(options, 'bCompact') && options.bCompact
        x = zeros(options.num_orig_vars, 1);
        x(this.I_active_variables) = xs;
      else
        x = xs;
      end
      assert(this.checkFeasible(x, ...
        struct('ConstraintTolerance', options.fmincon_opt.ConstraintTolerance)), ...
        'error: infeasible solution.');
    end
    

	end
		
end