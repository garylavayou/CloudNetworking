%% priceOptimalFlowRate
% Override <NormalSliceOptimizer>.<priceOptimalFlowRate>, Considering resource
% reconfiguration cost. This function has also implemented the parent method, so no need
% to call the parent method.
% 
% NOTE: update price before call this function.
function [output, loads, fval] = priceOptimalFlowRate(this, x0, options) %#ok<INUSL>
global DEBUG INFO;
if nargin <= 2
	options = Dictionary();
else
	options = Dictionary(options);
end
setdefault(options, struct('bParallel', false, 'bInitialize', false));
pardata = this.pardata;
Nvnf = pardata.NumberVNFs;
Nl = pardata.NumberLinks;
Nsn = pardata.NumberServiceNodes;


%% List of constaints
% (see also <DynamicSlice.fastReconfigure2>):
%		(1) flow conservation law: Nvf*Nan, <Equality>
%   (2) flow processing requirement: Nsn*Nvnf*Nf (num_lcon_res);
%   (3) virtual link capacity constraint: Nl (resource reservation);
%   (4) VNF instance capacity constraint: Nvnf*Nsn (num_varv); virtual node capacity constraint;
%   (5) link reconfiguration cost constraint: 2*num_varx;
%   (6) VNF allocation reconfiguration cost constraint: 2*num_varz;
%   (7) VNF instance capacity reconfiguration cost constraint: 2*num_varv;
%   (8) lower bound of virtual link capacity (optional);
%   (9) lower bound of virtual node capacity (dimconfig2): sum capacity of VNF isntance
%       larger than the lowbound;
%   (10)lower bound of VNF instance capacity (dimconfig1);
%   (11)resource utilization lower bound for resource reservation (optional);
%   (12)node capacity constraint with resource reservation (optional)
%
%           1   2   3   4   5   6   7   8   9  10  11
% InitRSV   +   +   +   -   -   -   +   -   -   -   +
% Redim     +   +   +   +   +   +   +   -   -   -   -
% RedimRSV  +   +   +   +   +   +   +   +  (+) (+)  -
% RedimRSV+ +   +   +   +   +   +   +   +  (+) (+)  +
%% Number of Variables
%   x: path; z: VNF instances assignment; r: user flow rate; v: VNF instances capacity;
%   tx;      tz;                                             tv;
%   c: link capacity [MUST be the LAST variable]
%           x   z   r   v   tx   tz   tv    c
% InitRSV   +   +   +   +    -    -    -    +
% Redim     +   +   +   +    +    +    +    +
% RedimRSV  +   +   +   +    +    +    +    +
% RedimRSV+ +   +   +   +    +    +    +    +
if options.bInitialize
	Nan = this.NumberAugmentedNodes;
	Nal = this.NumberAugmentedLinks;
	Nvf = pardata.NumberFlowSections;   % Removed flows not involved in optimization.
	Nf = pardata.NumberFlows;
	this.num_vars = [	Nvf*Nal; 	Nvnf*Nf*Nsn; 	Nf; Nvnf*Nsn];  % MUST be reset here in the continue updating process
	theta0 = this.GLOBAL_OPTIONS.theta0;
	options = structmerge(options, ...
		getstructfields(this.options, {'Reserve'}, 'default', {theta0}), 'exclude');
	if this.b_derive_vnf || this.options.ReconfigMethod == ReconfigMethod.DimBaseline
		% Initialization
		% For the fisrt time, or for the case that no reconfigration cost is considered.
		this.invoke_method = 1;
	elseif isempty(this.lower_bounds)  % bReserve is true
		this.invoke_method = 2;  % no resource reservation (theta = 0.99 to tolerate error)
	elseif options.Reserve >= 1
		this.invoke_method = 3;  % implicit resource reservation
	elseif options.Reserve > 0
		this.invoke_method = 4;  % explicit resource reservation
	elseif options.Reserve < 0 && options.Reserve > -1
		this.invoke_method = 5;
		options.Reserve = -options.Reserve;
	else
		error('%s: invalid value ''Reserve=%.2f''for resource reservation.', calledby, options.Reserve);
	end
	num_var_xzrv = sum(this.num_vars(1:4));
	num_vars_prim = sum(this.num_vars(1:3));
	num_varxz = sum(this.num_vars(1:2));
	num_varx = this.num_vars(1);
	num_varz = this.num_vars(2);
	num_varv = this.num_vars(4);
	num_lcon_proc = size(this.As_proc,1);
	num_lcon_le = num_lcon_proc + Nvnf*Nsn + num_varv;
	if this.invoke_method ~= 1
		this.num_vars = [this.num_vars; num_varx; num_varz; num_varv];
		num_lcon_le = num_lcon_le + 2*num_varx + 2*num_varz + 2*num_varv;
	end
	this.num_vars = [this.num_vars; Nl];
	num_vars = sum(this.num_vars);		% Total number of variables including temp ones.
	num_lcon_le = num_lcon_le + Nl;   % link capacity constraint
	
	if this.invoke_method == 3 || this.invoke_method == 4
		if isfield(this.lower_bounds, 'node')
			if isfield(this.options, 'ResidualCapacity')
				node_idx = this.options.ResidualCapacity.Node-this.lower_bounds.node>1;
			else
				node_idx = true(Nsn,1);
			end
			nn_bounds = nnz(node_idx);
			num_lcon_le = num_lcon_le + nn_bounds;
		end
	end
	if (this.invoke_method == 1 || this.invoke_method == 4) && this.options.bReserve
		% overall resource utilization constraint is mandated.
		theta = options.Reserve;
		num_lcon_le = num_lcon_le + 2;
	end
	if this.invoke_method == 1 || this.invoke_method == 4
		switch this.options.bReserve
			case 1 % only resource reservation to achieve overall utilization ratio
				theta1 = theta0;
			case 2 % resource reservation for individual resource (per node)
				theta1 = min((1+theta)/2, theta0);
				num_lcon_le = num_lcon_le + Nsn;
			case 3 % resource reservation for individual resource (per VNF instance)
				% add coefficient to constraint (2), no new constraints.
				theta1 = min((1+theta)/2, theta0);
				theta0 = theta1;
			case 0 % no resource reservation
				theta1 = theta0;
		end
	else
		% If |theta0| is set to 1, no resource is reserved for individual resource;
		% Otherwise, we set |theta0| very close to 1, to handle the precision error.
		theta1 = theta0;
	end
	
	%% Initialize Problem
	switch pardata.Weight
		case 50
			unit = 0.02;
		case 10
			unit = 0.05;
		case 300
			unit = 0.005;
		otherwise
			unit = 1;
	end
	this.problem = Dictionary('bCompact', true); prbm = this.problem;
	prbm.Aeq = sparse(Nan*Nvf, num_vars);  % (1) flow conservation law
	prbm.Aeq(:, 1:num_vars_prim) = [this.As_flow, sparse(Nan*Nvf, Nsn*Nvnf*Nf), -this.Ids];
	prbm.beq = sparse(Nan*Nvf, 1);
	prbm.A = sparse(num_lcon_le, num_vars);
	prbm.b = sparse(num_lcon_le,1);
	prbm.lb = sparse(num_vars, 1);
	row_offset = 0;
	prbm.A(1:num_lcon_proc, 1:(num_varxz)) = [this.As_proc, -this.As_procz]; % (2) processing reuirement
	row_offset = row_offset + num_lcon_proc;
	prbm.A(row_offset+(1:num_varv), num_varx+(1:num_varz)) = repmat(speye(num_varv), 1, Nf);  % (3) VNF capacity constraint
	prbm.A(row_offset+(1:num_varv), num_vars_prim+(1:num_varv)) = -theta0*speye(num_varv);    % VNF capacity is variables
	row_offset = row_offset + num_varv;
	prbm.A(row_offset+(1:Nl), 1:num_varx) = repmat(speye(Nl, Nal),1, Nvf); %% (4) virtual link capacity constraint
	prbm.A(row_offset+(1:Nl), num_vars-Nl+(1:Nl)) = -theta1*speye(Nl);  % resource reservation for individual links
	row_offset = row_offset + Nl;
	if this.invoke_method~=1
		prbm.A(row_offset+(1:num_varx), 1:num_varx) = speye(num_varx); 	% (5) link reconfiguration cost constraint
		prbm.A(row_offset+(1:num_varx), num_var_xzrv+(1:num_varx)) = -speye(num_varx);
		prbm.b(row_offset+(1:num_varx)) = this.topts.old_variables.x/unit; % which have been pre-processed, so it can be compared with the current states.
		row_offset = row_offset + num_varx;
		prbm.A(row_offset+(1:num_varx), 1:num_varx) = -speye(num_varx);
		prbm.A(row_offset+(1:num_varx), num_var_xzrv+(1:num_varx)) = -speye(num_varx);
		prbm.b(row_offset+(1:num_varx)) = -this.topts.old_variables.x/unit;
		row_offset = row_offset + num_varx;   	% (6) VNF allocation reconfiguration cost constraint =>
		prbm.A(row_offset+(1:num_varz), (num_varx+1):num_varxz) = speye(num_varz);
		prbm.A(row_offset+(1:num_varz), (num_var_xzrv+num_varx)+(1:num_varz)) = -speye(num_varz);
		prbm.b(row_offset+(1:num_varz)) = this.topts.old_variables.z/unit;
		row_offset = row_offset + num_varz;
		prbm.A(row_offset+(1:num_varz), (num_varx+1):num_varxz) = -speye(num_varz);
		prbm.A(row_offset+(1:num_varz), (num_var_xzrv+num_varx)+(1:num_varz)) = -speye(num_varz);
		prbm.b(row_offset+(1:num_varz)) = -this.topts.old_variables.z/unit;
		row_offset = row_offset + num_varz;   	% (7) VNF instance capacity reconfiguration cost constraint =>
		prbm.A(row_offset+(1:num_varv), num_vars_prim+(1:num_varv)) = speye(num_varv);
		prbm.A(row_offset+(1:num_varv), (num_var_xzrv+num_varxz)+(1:num_varv)) = -speye(num_varv);
		prbm.b(row_offset+(1:num_varv)) = this.old_variables.v/unit;
		row_offset = row_offset + num_varv;
		prbm.A(row_offset+(1:num_varv), num_vars_prim+(1:num_varv)) = -speye(num_varv);
		prbm.A(row_offset+(1:num_varv), (num_var_xzrv+num_varxz)+(1:num_varv)) = -speye(num_varv);
		prbm.b(row_offset+(1:num_varv)) = -this.old_variables.v/unit;
		row_offset = row_offset + num_varv;
	end
	if this.invoke_method == 3 || this.invoke_method == 4  % (8) lower bound of node capacity
		if isfield(this.lower_bounds, 'node')
			reduced_eye = speye(Nsn);
			reduced_eye = reduced_eye(node_idx,:);
			prbm.A(row_offset+find(node_idx), num_vars_prim+(1:num_varv)) = -repmat(reduced_eye,1,Nvnf);
			prbm.b(row_offset+find(node_idx)) = -this.lower_bounds.node(node_idx)/unit;
			row_offset = row_offset + nn_bounds;
		end
	end
	%% (10) resource utilization lower bound for resource reservation
	I_active_vars = this.updateActiveVariables();  % xzr
	if this.invoke_method == 1 || this.invoke_method == 4
		if this.options.bReserve  % overall resource utilization constraint
			prbm.A(row_offset+1, [num_varx+(1:num_varz), num_vars_prim+(1:num_varv)]) = ...
				[this.I_active_variables((num_varx+1):num_varxz), -theta*ones(1,num_varv)];  % z_nvf => (node, function, flow)
			prbm.A(row_offset+2, [1:num_varx, (num_vars-Nl+1):num_vars]) = ...
				[~this.I_fake_edgevars', -theta*ones(1, Nl)];  % masked later by this.I_active_variables(1:num_varx)
			row_offset = row_offset + 2;
		end
		if this.options.bReserve == 2 % per node
			%% (11) virtual node capacity constraint for resource reservation
			% (optional) resource reservation for individual nodes: node's capacity is the sum of
			% its VNF capacity.
			prbm.A(row_offset+(1:Nsn), num_varx+(1:num_varz)) = repmat(speye(Nsn), 1, Nvnf*Nf);
			prbm.A(row_offset+(1:Nsn), num_vars_prim+(1:num_varv)) = -theta1*repmat(speye(Nsn),1,Nvnf);
			% bs(row_offset+(1:NC)) = zeros(NC,1);
			row_offset = row_offset + Nsn; %#ok<NASGU>
		elseif this.options.bReserve == 3
			warning('incomplete implementation!');
		end
		% (per VNF instance is appied to constriant (2))
	end
	if (this.invoke_method == 3 || this.invoke_method == 4)  %% (7)(9)  lower-bounds
		if isfield(this.options, 'ResidualCapacity')
			link_idx = this.options.ResidualCapacity.Link-this.lower_bounds.link>1;
		else
			link_idx = true(Nl,1);
		end
		prbm.lb(num_vars-Nl+find(link_idx)) = this.lower_bounds.link(link_idx)/unit;
		if isfield(this.lower_bounds, 'VNF')
			vnf_idx = repmat(node_idx, Nvnf, 1);
			prbm.lb(num_vars_prim+find(vnf_idx)) = this.lower_bounds.VNF(vnf_idx)/unit;
		end
	end
	
	prbm.x0 = sparse(num_vars,1);
	assert(this.checkFeasible(prbm.x0/unit), 'error: infeasible start point.');
	
	%% Compress
	% eleminate the columns (unused variables)
	if this.invoke_method == 1
		I_active_vars = [I_active_vars, true(1,num_varv+Nl)]; %xzrvc
		prbm.A(:, ~logical(I_active_vars)) = 0;		% Aeq not been updated since initialize
	else
		I_active_vars = [I_active_vars, true(1,num_varv), I_active_vars, true(1,num_varv+Nl)];
		prbm.A(:, ~I_active_vars) = 0;  % xzrvxzvc  % Aeq not been updated since initialize
		% eleminate the rows (fake variables): not related to reconfiguration cost
		row_offset = num_lcon_proc + num_varv + Nl;
		idx_unused_rows = [row_offset+find(this.topts.x_reconfig_cost==0); ...
			row_offset+num_varx+find(this.topts.x_reconfig_cost==0);...
			row_offset+num_varx*2+find(this.topts.z_reconfig_cost==0);...
			row_offset+num_varx*2+num_varz+find(this.topts.z_reconfig_cost==0)];
		prbm.A(idx_unused_rows, :) = 0;		% CLEARED
		prbm.b(idx_unused_rows, :) = 0;
	end
	I_unused_rows = sum(abs(prbm.A), 2)==0;
	% No need to assert, see also <fastReconfigure>.
	% assert(isempty(find(prbm.b(I_unused_rows)>eps,1)), 'error: infeasible inequality constraints.');
	% assert(isempty(find(prbm.beq(I_unused_rows)>eps,1)), 'error: infeasible equality constraints.');
	prbm.A(I_unused_rows, :) = [];		% REMOVED
	prbm.b(I_unused_rows, :) = [];
	I_unused_rows = sum(abs(prbm.Aeq), 2)==0;
	prbm.Aeq(I_unused_rows, :) = [];		% REMOVED
	prbm.beq(I_unused_rows, :) = [];
	% Determine the active variables in the problem, MUST use [A; Aeq] both.
	I_active_variables = sparse(sum(abs([prbm.A; prbm.Aeq]),1)~=0);
	% NOTE that the |reconfig_cost| cosefficient is not filtered like the
	% <fastReconfigure>. The objective function <fcnProfitReconfigureSlicing> will use the
	% orginal cost coefficient.
	prbm.Aeq = prbm.Aeq(:, I_active_variables);
	prbm.A = prbm.A(:, I_active_variables);
	prbm.x0 = prbm.x0(I_active_variables);
	prbm.lb = prbm.lb(I_active_variables);
	prbm.num_vars = [nnz(I_active_variables(1:num_varx));
		nnz(I_active_variables((num_varx+1):num_varxz));
		nnz(I_active_variables((num_varxz+1):num_vars_prim));
		nnz(I_active_variables(num_vars_prim+(1:num_varv)))];
	if this.invoke_method ~= 1
		prbm.num_vars = [prbm.num_vars;
			nnz(I_active_variables(num_var_xzrv+(1:num_varx)));
			nnz(I_active_variables(num_var_xzrv+num_varx+(1:num_varz)));
			nnz(I_active_variables(num_var_xzrv+num_varxz+(1:num_varv)))];
	end
	prbm.num_vars = [prbm.num_vars; nnz(I_active_variables(num_vars-Nl+(1:Nl)))];
	% options.num_orig_vars = num_vars;   % original number of variables in the problem
	if this.invoke_method ~= 1
		% The |x_reconfig_cost| also include the items for the "fake variables", see also
		% <update_reconfig_costinfo>.
		% As the |I_active_vars| masks the unused varaiables |x| and |z|, the cost vector will
		% be significantly shrinked.
		this.topts.x_reconfig_cost = this.topts.x_reconfig_cost(...
			I_active_variables(num_var_xzrv+(1:num_varx)) );
		this.topts.z_reconfig_cost = this.topts.z_reconfig_cost(...
			I_active_variables(num_var_xzrv+num_varx+(1:num_varz)) );
		% 'vnf_reconfig_cost' cannot be compressed
	end

	%% Optimization options
	% * *Algorithm* : since the problem contains linear constraints and bound
	% constraints, then |trust-region-reflective| method is not applicable. Hence,
	% we choose the |interior point| method. As a result the Hessian matrix should
	% be computed separately.
	% is directly returned from the objective function as the second derivatives.
	% * *HessianFcn* : we compute Hessian using the objective function.
	% Therefore, this option is set as |'objective'|.
	% * *SpecifyObjectiveGradient* : the gradient can be directly computed from
	% the objective function, so this option is set to |true|.
	% * *SpecifyConstraintGradient* : since this problem does not contain nonlinear
	% constraint, this option is set to |false|.
	% * *Display information* : use |'iter'| to display iteration information for
	% debug. use |'notify-detailed'| to only display exception message.
	prbm.minopts = optimoptions(@fmincon);
	prbm.minopts.Algorithm = 'interior-point';
	prbm.minopts.SpecifyObjectiveGradient = true;
	prbm.minopts.MaxIterations = 60;
	prbm.minopts.MaxFunctionEvaluations = 300;
	prbm.minopts.Display = 'iter';   %'notify-detailed'; %'iter'; 'notify'
	% minopts.OptimalityTolerance = 10^-5;
	% minopts.CheckGradients = true;
	% minopts.FiniteDifferenceType = 'central';
	% minopts.Diagnostics = 'on';
	%% Save information
	prbm.Reserve = options.Reserve;
	prbm.unit = unit;
	prbm.I_active_full_variables = I_active_variables;   % not update to the optimizer
	if options.bParallel  % update to the optimizer
		prbm.num_full_vars = this.num_vars;  
		prbm.I_active_variables = this.I_active_variables;	 
		prbm.invoke_method = this.invoke_method;
	end
else
	prbm = this.problem(1);
	num_varx = this.num_vars(1);
	num_varz = this.num_vars(2);
	num_varxz = sum(this.num_vars(1:2));
	num_varv = this.num_vars(4);
	num_vars_prim = sum(this.num_vars(1:3));
	options.Reserve = prbm.Reserve;
	I_active_variables = prbm.I_active_full_variables;
end
options = structupdate(options, prbm, {'bCompact', 'unit'});
options = setdefault(options, this.pardata, {'PricingPolicy'});
options.bFinal = false;
x0 = prbm.x0/options.unit;
if this.invoke_method == 1
	prbm.minopts.HessianFcn = ...
		@(x,lambda)NormalDynamicSliceOptimizer.hessSlicing(x, lambda, this, options);
	% 'hessSlicing' is the same as hessInitialSlicing as we carefully treat the link
	% capacity variables.
	[xs, fval, exitflag, foutput] = ...
		fmincon(@(x)NormalDynamicSliceOptimizer.fcnProfitReserveSlicing(x, this, options), ...
		x0, prbm.A, prbm.b, prbm.Aeq, prbm.beq, prbm.lb, [], [], prbm.minopts);
else
	prbm.minopts.HessianFcn = ...
		@(x,lambda)NormalDynamicSliceOptimizer.hessSlicing(x, lambda, this, options);
	[xs, fval, exitflag, foutput] = ...
		fmincon(@(x)NormalDynamicSliceOptimizer.fcnProfitReconfigureSlicing(x, this, options), ...
		x0, prbm.A, prbm.b, prbm.Aeq, prbm.beq, prbm.lb, [], [], prbm.minopts);
	if exitflag < 0 
		warning([strtok(foutput.message, newline), 'Retry.']);
		x0= [this.topts.old_variables.x;
			this.topts.old_variables.z;
			this.topts.old_variables.r;
			this.old_variables.v;
			this.topts.old_variables.x;
			this.topts.old_variables.z;
			this.old_variables.v;
			this.old_variables.c];
		x0 = x0(I_active_variables)/options.unit;
		[xs, fval, exitflag, foutput] = ...
			fmincon(@(x)NormalDynamicSliceOptimizer.fcnProfitReconfigureSlicing(x, this, options), ...
			x0, prbm.A, prbm.b, prbm.Aeq, prbm.beq, prbm.lb, [], [], prbm.minopts);
	end
end
xs = xs*options.unit;
num_vars = sum(this.num_vars);
if options.bCompact
	x = sparse(num_vars, 1);
	x(I_active_variables) = xs;
else
	x = xs;
end
this.interpretExitflag(exitflag, foutput.message);
if (~isempty(DEBUG) && DEBUG) || (~isempty(INFO) && INFO)
	fprintf('\tThe optimal net profit of the slice: %G.\n', -fval);
end

%% Output solution
options.ConstraintTolerance = prbm.minopts.ConstraintTolerance;
assert(this.checkFeasible(x, options), 'error: infeasible solution.');
output = Dictionary();
output.iterations = foutput.iterations;
output.funcCount = foutput.funcCount;
%%%
% The post processing like <optimalFlowRate> is not needed, since the
% objective function (reconfiguration cost) will force those variables to be zero.
num_varr = this.num_vars(3);
this.temp_vars = struct('x', x(1:num_varx), ...
	'z', x(num_varx+(1:num_varz)), ...
	'r', x(num_varxz+(1:num_varr)), ...
	'v', x(num_vars_prim+(1:num_varv)));
if this.invoke_method >= 2
	num_var_xzrv = sum(this.num_vars(1:4));
	this.temp_vars.tx = x(num_var_xzrv+(1:num_varx));
	this.temp_vars.tz = x(num_var_xzrv+num_varx+(1:num_varz));
	this.temp_vars.tv = x(num_var_xzrv+num_varxz+(1:num_varv));
end
this.temp_vars.c = x(num_vars-Nl+(1:Nl));

%% perform resource reservation
if this.invoke_method == 5
	if options.bInitialize
		warning('Debug Required!');
		prbm = Dictionary();
		prbm.num_vars = [Nl;			Nsn*Nvnf;			Nsn*Nvnf];
		
		theta = options.Reserve;
		loads.Node = this.getNodeLoad(false);
		loads.Link = this.getLinkLoad(false);
		prbm.x0r = [this.temp_vars.c;
			this.old_variables.v;		% v
			zeros(Nsn*Nvnf,1)];			% tv
		prbm.lbr = [max(this.lower_bounds.link, loads.Link);      % options.bDisableLowerBound
			this.getVNFLoad();
			zeros(Nsn*Nvnf,1)];
		%%
		num_vars = sum(prbm.num_vars);
		num_lcons = Nsn + 2 + 2*Nsn*Nvnf;
		row_offset = 0;
		prbm.Ar = sparse(num_lcons, num_vars);
		prbm.br = sparse(num_lcons, 1);
		prbm.br(1:Nsn) = -max(this.lower_bounds.node, loads.Node);    % options.bDisableLowerBound
		prbm.Ar(1:Nsn,Nl+(1:Nsn*Nvnf)) = -repmat(speye(Nsn),1,Nvnf);
		row_offset = row_offset + Nsn;
		prbm.br(row_offset+1) = -sum(loads.Link)/theta;
		prbm.Ar(row_offset+1,1:Nl) = -ones(1,Nl);
		row_offset = row_offset + 1;
		prbm.br(row_offset+1) = -sum(loads.Node)/theta;
		prbm.Ar(row_offset+1,Nl+(1:Nsn*Nvnf)) = -ones(1,Nsn*Nvnf);
		row_offset = row_offset + 1;
		prbm.Ar(row_offset+(1:Nsn*Nvnf), Nl+(1:2*Nsn*Nvnf)) = [speye(Nsn*Nvnf), -speye(Nsn*Nvnf)];
		prbm.br(row_offset+(1:Nsn*Nvnf)) = this.old_variables.v;
		row_offset = row_offset + Nsn*Nvnf;
		prbm.Ar(row_offset+(1:Nsn*Nvnf), Nl+(1:2*Nsn*Nvnf)) = [-speye(Nsn*Nvnf), -speye(Nsn*Nvnf)];
		prbm.br(row_offset+(1:Nsn*Nvnf)) = -this.old_variables.v;
		%     row_offset = row_offset + Nsn*Nvnf;
		%%
		% we might add a resource utilization ratio for individual resources.
		prbm.minopts = optimoptions(@fmincon);
		prbm.minopts.Algorithm = 'interior-point';
		prbm.minopts.SpecifyObjectiveGradient = true;
		prbm.minopts.Display = 'iter';   %'notify-detailed'; %'iter'; 'notify'
		%     fr_opt.CheckGradients = true;
		%     fr_opt.FiniteDifferenceType = 'central';
		this.problem(2) = prbm;
	else
		prbm = this.problem(2);
	end
	%%
	% we might add a resource utilization ratio for individual resources.
	prbm.minopts.HessianFcn = @(x,lambda)hessReconfigCost(x, lambda, this, options);
	[x, ~, exitflag, foutput] = ...
		fmincon(@(x)fcnReconfigCost(x, this, options), ...
		prbm.x0r, prbm.Ar, prbm.br, [], [], prbm.lbr, [], [], prbm.minopts);
	this.interpretExitflag(exitflag, foutput.message);
	this.temp_vars.c = x(1:Nl);
	this.temp_vars.v = x(Nl+(1:Nsn*Nvnf));
	this.temp_vars.tv = x(Nl+Nsn*Nvnf+(1:Nsn*Nvnf));
	output.iterations = output.iterations + foutput.iterations;
	output.funcCount = output.funcCount + foutput.funcCount;
end

this.x0 = x;
output.net_profit = -fval;
if options.bParallel
	output.x0 = x;
	output.temp_vars = this.temp_vars;
	if options.bInitialize
		output.problem = prbm;
	end
end
end

function [f,g] = fcnReconfigCost(vars, this, options)
Nsn = this.pardata.NumberServiceNodes;
Nl = this.pardata.NumberLinks;
Nvnf = this.pardata.NumberVNFs;
link_cap = vars(1:Nl);
vnf_cap = vars(Nl+(1:Nsn*Nvnf));
vnf_cap_diff = vars(Nl+Nsn*Nvnf+(1:Nsn*Nvnf));
load.Node = sum(reshape(vnf_cap, Nsn, Nvnf),2);
switch options.PricingPolicy
	case {'quadratic-price', 'quadratic'}
		[link_payment,link_price_grad] = this.fcnLinkPricing(this.prices.Link, link_cap);
		[node_payment,node_price_grad] = this.fcnNodePricing(this.prices.Node, load.Node);
		f = link_payment + node_payment;
	otherwise
		error('%s: invalid pricing policy', calledby);
end
f = f + dot(vnf_cap_diff, this.topts.vnf_reconfig_cost);

g = [link_price_grad; repmat(node_price_grad, Nvnf, 1); this.topts.vnf_reconfig_cost];

end

function h = hessReconfigCost(vars, lambda, this, options) %#ok<INUSL>
Nsn = this.pardata.NumberServiceNodes;
Nl = this.pardata.NumberLinks;
Nvnf = this.pardata.NumberVNFs;
link_cap = vars(1:Nl);
vnf_cap = vars(Nl+(1:Nsn*Nvnf));
load.Node = sum(reshape(vnf_cap, Nsn, Nvnf),2);
num_vars = length(vars);
h = spalloc(num_vars,num_vars, num_vars);
switch options.PricingPolicy
	case {'quadratic-price', 'quadratic'}
		[~,~,lph] = this.fcnLinkPricing(this.prices.Link, link_cap);
		[~,~,nph] = this.fcnNodePricing(this.prices.Node, load.Node);
		h(1:Nl,1:Nl) = diag(lph);
		h(Nl+(1:Nsn*Nvnf),Nl+(1:Nsn*Nvnf)) = block_diag(diag(nph), Nvnf);
	otherwise
		error('%s: invalid pricing policy', calledby);
end
end