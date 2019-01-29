%%
% We introduce the quadratic penalty item of the dual-ADMM method.
% The variables are the flow-edge variables  and the flow-node
% variables.
%
function [gamma, fval, output] = priceOptimalFlowRateDA(this, x0, lambda, q, options)
global DEBUG INFO; %#ok<NUSED>
if nargin <= 2
	options = Dictionary();
else
	options = Dictionary(options);
end
setdefault(options, struct('bParallel', false, 'bInitialize', false));
pardata = this.pardata;
Nsn = pardata.NumberServiceNodes;
Nf = pardata.NumberFlows;

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
		prbm.minopts.ConstraintTolerance = 1e-4;
		prbm.minopts.MaxIterations = 60;
		% prbm.minopts.SubproblemAlgorithm = 'cg'; % only take cg steps.
		% 		prbm.minopts.CheckGradients = true;
		% 		prbm.minopts.FiniteDifferenceType = 'central';
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
	prbm.minopts.Display = 'off';   %{'notify-detailed'|'notify'|'iter'|'off'}
end

%%
lambda = lambda(pardata.DualIndex);
q = q(pardata.DualIndex);
if strcmpi(this.options.OptimizationTool, 'matlab')
	if nargin >= 2 && ~isempty(x0)
		prbm.x0 = x0(this.I_active_variables);
	else
		%var0 = this.x0;
		% var0 = rand(this.num_flow_vars,1);
		prbm.x0 = sparse(sum(prbm.num_vars),1);
	end
	if strcmpi(prbm.minopts.Algorithm, 'interior-point')
		% The method could be overload by subclasses, so invoking it via class name
		% should be safe.
		prbm.minopts.HessianFcn = ...
			@(x,lbd)NormalSliceOptimizer.fcnHessian(x, lbd, this, lambda, q, options);
	end
	warning('off')
	[xs, fval, exitflag, foutput] = ...
		fmincon(@(x)NormalSliceOptimizer.fcnProfitPrimal(x, this, lambda, q, options), ...
		prbm.x0, prbm.A, prbm.b, prbm.Aeq, prbm.beq, prbm.lb, [], [], prbm.minopts);
	warning('on')
	this.interpretExitflag(exitflag, foutput);
elseif strcmpi(this.options.OptimizationTool, 'cvx')
	[xs, fval, exitflag, foutput] = cvx_method(lambda, q);
	if exitflag<1
		error('error: CVX failed to solve the problem.');
	end
else
	error('error: Un-recognized optimization tool.');
end
%% CVX
	function [xs, fval, exitflag, output] = cvx_method(lambda, q)
		w = pardata.Weight;
		p = [this.prices.Link; this.prices.Node];
		plt = options.InterSlicePenalty;
		c = [this.capacities.Link; this.capacities.Node];
		Al = this.As_load(:, prbm.I_active_xz);
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
		%{
					cvx_begin
					variable xs(n,1);
					%maximize( w*sum(log(1+xs((nxz+1):n)))- p'*Al*xs(1:nxz)-plt/2*norm(max(0,1/plt*(Al*xs(1:nxz)-c+q)))^2 )
					maximize( w*sum(log(1+xs((nxz+1):n)))- p'*Al*xs(1:nxz)-plt/2*sum(pow_pos(lambda+1/plt*(Al*xs(1:nxz)-c+q),2)) )
					subject to
					this.problem.A * xs <= this.problem.b;  %#ok<VUNUS>
					this.problem.Aeq * xs == this.problem.beq;  %#ok<EQEFF>
					xs >= lbs; %#ok<VUNUS>
					cvx_end
		%}
		c_x = []; c_z = []; c_r = [];
		cvx_begin
		variable c_x(c_nx,1);
		variable c_z(c_nz,1);
		variable c_r(c_nr,1);
		minimize( -w*sum(log(1+c_r))+p'*Al*[c_x;c_z]+plt/2*sum(pow_pos(lambda+1/plt*(Al*[c_x;c_z]-c+q),2)) )
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
x = sparse(num_vars, 1);
x(this.I_active_variables) = xs;
this.x0 = x;
this.temp_vars.x = x(1:num_varx);
this.temp_vars.z = x(num_varx+(1:num_varz));
num_varr = this.num_vars(3);
this.temp_vars.r = x(num_varx+num_varz+(1:num_varr));
output = Dictionary();
if options.bParallel
	output.x0 = x;
	output.temp_vars = this.temp_vars;
	if options.bInitialize
		output.problem = prbm;
	end
end
output.net_profit = -fval;
output.iterations = foutput.iterations;
if strcmpi(this.options.OptimizationTool, 'cvx')
	output.funcCount = foutput.iterations;
else
	output.funcCount = foutput.funcCount;
end
gamma = zeros(pardata.NumberDuals,1);
num_varxz = sum(prbm.num_vars(1:2));
gamma(pardata.DualIndex) = max(0, ...
	NormalSliceOptimizer.fcnPenalty(xs(1:num_varxz), this, lambda, q, options.InterSlicePenalty));
end