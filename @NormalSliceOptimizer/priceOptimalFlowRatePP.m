function [gamma, fval, x_prox, output] = ...
	priceOptimalFlowRatePP(this, x0, lambda, q, xp, options)
global DEBUG INFO; %#ok<NUSED>
if nargin <= 2
	options = Dictionary();
else
	options = Dictionary(options);
end
pardata = this.pardata;
if options.InterSlicePenalty == 0
	options.InterSlicePenalty = pardata.Weight;
end
setdefault(options, struct('bParallel', false, 'bInitialize', false));
Nsn = pardata.NumberServiceNodes;
Nf = pardata.NumberFlows;
%%
idx = [pardata.PhysicalLink; length(pardata.PhysicalLink)+pardata.PhysicalDataCenter];
lambda = lambda(idx);
q = q(idx);

if options.bInitialize
	Nal = this.NumberAugmentedLinks;
	Nvf = pardata.NumberFlowSections;
	Nvnf = pardata.NumberVNFs;
	num_lcon_res = size(this.As_proc, 1);
	num_lcon_flow = size(this.As_flow,1);
	
	this.num_vars = [Nvf*Nal; Nvnf*Nf*Nsn; Nf; Nvnf*Nsn]; % MUST be reset here in the continue updating process
	num_varx = this.num_vars(1);
	num_varz = this.num_vars(2);
	num_varxz = sum(this.num_vars(1:2));
	num_vars = sum(this.num_vars);
	%% precondigtioning
	% see also <priceOptimalFlowRate>.
	unit = 1;
	%% Flow reservation constraints and processing resource constraints
	% see also <priceOptimalFlowRate>.
	this.problem = Dictionary('bCompact', true); prbm = this.problem;
	prbm.Aeq = [this.As_flow, sparse(num_lcon_flow, num_varz), -this.Ids];
	prbm.A = [this.As_proc, -this.As_procz, sparse(num_lcon_res, Nf)];
	prbm.beq = sparse(num_lcon_flow, 1);
	prbm.b = sparse(num_lcon_res, 1);
	this.updateActiveVariables();
	prbm.Aeq = prbm.Aeq(:, this.I_active_variables);
	prbm.A = prbm.A(:, this.I_active_variables);
	I_effect_rows = sum(abs(prbm.Aeq), 2)~=0;
	prbm.Aeq = prbm.Aeq(I_effect_rows,:);
	prbm.beq = prbm.beq(I_effect_rows,:);
	I_effect_rows = sum(abs(prbm.A), 2)~=0;
	prbm.A = prbm.A(I_effect_rows,:);
	prbm.b = prbm.b(I_effect_rows,:);
	
	prbm.num_vars = [...
		nnz(this.I_active_variables(1:num_varx));...
		nnz(this.I_active_variables(num_varx+1 : num_varxz));...
		Nf];  % this.num_vars is the original number of variables
	prbm.lb = sparse(sum(prbm.num_vars),1);
	if strcmpi(this.options.OptimizationTool, 'matlab')
		%% Set the optimization options		
		prbm.minopts = optimoptions('fmincon');
		prbm.minopts.Algorithm = 'interior-point';
		%prbm.minopts.Algorithm = 'sqp';
		prbm.minopts.SpecifyObjectiveGradient = true;
		prbm.minopts.OptimalityTolerance = 1e-5;
		prbm.minopts.ConstraintTolerance = 1e-3;
		prbm.minopts.MaxIterations = 100;
		% prbm.minopts.SubproblemAlgorithm = 'cg'; % only take cg steps.
		%prbm.minopts.CheckGradients = true;
		%prbm.minopts.FiniteDifferenceType = 'central';
	end
	%% Save information
	prbm.unit = unit;
	if options.bParallel  % update to the optimizer
		prbm.num_full_vars = this.num_vars;
		prbm.I_active_variables = this.I_active_variables;
	end
else
	prbm = this.problem;
	num_vars = sum(this.num_vars);
	num_varx = this.num_vars(1);
	num_varz = this.num_vars(2);
end
if isfield(options, 'Display')
	prbm.minopts.Display = options.Display;
else
	prbm.minopts.Display = 'off';   %{'notify-detailed'|'notify'|'iter'}
end

%%
lambda = lambda(pardata.DualIndex);
q = q(pardata.DualIndex);
xp = xp(this.I_active_variables);
if strcmpi(this.options.OptimizationTool, 'matlab')	
	if nargin >= 2 && ~isempty(x0)
		prbm.x0 = x0;
	else
		%var0 = this.x0;
		% var0 = rand(this.num_flow_vars,1);
		%var0 = zeros(this.num_flow_vars,1);
		prbm.x0 = xp;
	end
	options.TrueValue = false;
	if strcmpi(prbm.minopts.Algorithm, 'interior-point')
		prbm.minopts.HessianFcn = ...
			@(x,lbd)NormalSliceOptimizer.fcnHessianPP(x, lbd, this, lambda, q, options);
	end
	warning('off')
	[xs, fval, exitflag, foutput] = ...
		fmincon(@(x)NormalSliceOptimizer.fcnProfitPrimalPP(x, this, lambda, q, xp, options), ...
		prbm.x0, prbm.A, prbm.b, prbm.Aeq, prbm.beq, prbm.lb, [], [], prbm.minopts);
	warning('on')
	this.interpretExitflag(exitflag, foutput);
elseif strcmpi(this.options.OptimizationTool, 'cvx')
	[xs, fval, exitflag, foutput] = cvx_method(lambda, q, xp);
	if exitflag<1
		error('error: CVX failed to solve the problem.');
	end
else
	error('error: Un-recognized optimization tool.');
end

%% CVX
	function [xs, fval, exitflag, output] = cvx_method(lambda, q, xp)
		w = pardata.Weight;
		p = [this.prices.Link; this.prices.Node];
		r = options.InterSlicePenalty;
		Ns = options.NumberSlices;
		c = [this.capacities.Link; this.capacities.Node];
		Al = this.As_load(:,this.problem.I_active_xz);
		cvx_problem = []; cvx_optbnd = []; %#ok<NASGU>
		cvx_optval = []; cvx_status = [];
		cvx_slvitr = []; cvx_slvtol = []; cvx_cputime = []; %#ok<NASGU>
		if options.bCompact
			c_nx = this.num_vars(1);
			c_nz = this.num_vars(2);
			c_nr = this.num_vars(3);
		else
			c_nx = this.num_vars(1);
			c_nz = this.num_vars(2);
			c_nr = this.num_vars(3);
		end
		c_nxz = c_nx + c_nz;
		c_x = []; c_z = []; c_r = [];
		cvx_begin
		variable c_x(c_nx,1);
		variable c_z(c_nz,1);
		variable c_r(c_nr,1);
		minimize( -w*sum(log(1+c_r))+p'*Al*[c_x;c_z]+ 1/(2*r*Ns^2)*sum(pow_abs([c_x;c_z;c_r]-xp,2)) + r/2*sum(pow_pos(lambda+1/r*(Al*[c_x;c_z]-c+q),2)) )
		subject to
		prbm.A(:,1:c_nxz) * [c_x;c_z] <= 0;  %#ok<VUNUS>
		prbm.Aeq(:,[1:c_nx, c_nxz+(1:c_nr)]) * [c_x;c_r] == 0;  %#ok<EQEFF>
		[c_x;c_z;c_r] >= 0; %#ok<VUNUS>
		cvx_end
		fval = cvx_optval;
		xs = [c_x;c_z;c_r];
		switch cvx_status
			case 'Solved'
				exitflag = 1;
			case 'Inaccurate/Solved'
				exitflag = 2;
			case 'Failed'
				exitflag = -1;
			otherwise
				exitflag = 0;
		end
		output.iterations = cvx_slvitr;
	end

%% output solution
% assert(this.checkFeasible())
xs(xs<10^-4*norm(xs)) = 0;
x_prox = zeros(num_vars, 1);
x_prox(this.I_active_variables) = xs;
this.x0 = x_prox;
this.temp_vars.x = x_prox(1:num_varx);
this.temp_vars.z = x_prox(num_varx+(1:num_varz));
num_varr = this.num_vars(3);
this.temp_vars.r = x_prox(num_varx+num_varz+(1:num_varr));
output = Dictionary();
if options.bParallel
	output.x0 = x_prox;
	output.temp_vars = this.temp_vars;
	if options.bInitialize
		output.problem = prbm;
	end
end
options.TrueValue = true;
fval = fval + NormalSliceOptimizer.fcnProfitPrimalPP(xs, this, lambda, q, xp, options);
output.net_profit = -fval;
output.iterations = foutput.iterations;
if strcmpi(this.options.OptimizationTool, 'cvx')
	output.funcCount = foutput.iterations;
else
	output.funcCount = foutput.funcCount;
end

gamma = zeros(pardata.NumberDuals,1);
num_varxz = sum(prbm.num_vars(1:2));
gamma(pardata.DualIndex) = max(0, NormalSliceOptimizer.fcnPenalty(xs(1:num_varxz), this, lambda, q, ...
	options.InterSlicePenalty));
end