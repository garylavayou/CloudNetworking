function output = priceOptimalFlowRate(this, x0, options)
global DEBUG INFO; %#ok<NUSED>
if nargin <= 2
	options = Dictionary();
else
	options = Dictionary(options);
end
setdefault(options, struct('bParallel', false, 'bInitialize', false));
pardata = this.pardata;
Nl = pardata.NumberLinks;
Nsn = pardata.NumberServiceNodes;
Nf = pardata.NumberFlows;
% Number of variables in the problem, which has been initialized. See also
% <inititalState>, <onAddlingFlow>, <onRemovingFlow>.  

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
	% Slice     54   44
	%	Weight    50   10
	% Unit    0.02  0.01
	if ~isfield(options, 'unit')
		switch pardata.Weight
			case 10
				unit = 0.2;
			case 25
				unit = 0.03;
			case 50
				unit = 0.02;
			case 300
				unit = 0.01;
			otherwise
				unit = 1;
		end			%%
	else
		unit = options.unit;
	end
	
	%% Flow reservation constraints and processing resource constraints
	%		[As -Ir] * [F; r]  = 0   ==> As*F  = r
	%		[Al -Iz] * [F; z] <= 0   ==> Al*F <= z
	% ===>
	%								 ┌F┐										┌F┐
	%		[As 0 -Ir] * |z| = 0,  [Al -Iz 0] * |z| <= 0
	%								 └r┘										└r┘
	%
	% Several situation should be considered:
	%	1. intermediate forwarding nodes;
	%	2. source forwarding nodes of the first segment;
	% 3. destiniation forwarding nodes of the last segment;
	%	4. data center nodes that potentially process flows.
	this.problem = Dictionary('bCompact', true); prbm = this.problem;
	prbm.Aeq = [this.As_flow, sparse(num_lcon_flow, num_varz), -this.Ids];
	prbm.A = [this.As_proc, -this.As_procz, sparse(num_lcon_res, Nf)];
	prbm.beq = sparse(num_lcon_flow, 1);
	prbm.b = sparse(num_lcon_res, 1);
	
	%% Compact
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
	if nargin >= 2 && ~isempty(x0)
		prbm.x0 = x0;
	else
		prbm.x0 = sparse(sum(prbm.num_vars),1);
	end
	prbm.lb = sparse(sum(prbm.num_vars),1);
	
	if isfield(options, 'CapacityConstrained') && options.CapacityConstrained
		num_varxz = sum(prbm.num_vars(1:2));
		prbm.A(num_lcon_res+(1:Nl+Nsn), 1:num_varxz) = this.As_load(:,prbm.I_active_xz);
		prbm.b(num_lcon_res+(1:Nl+Nsn)) = [this.capacities.Link; this.capacities.Node]/unit;
		idx_off = find(sum(this.As_load(:,prbm.I_active_xz),2)==0);
		prbm.A(num_lcon_res+idx_off, :) = [];
		prbm.b(num_lcon_res+idx_off) = [];
	end
	
	%% Set the optimization options
	prbm.minopts = optimoptions(@fmincon);
	prbm.minopts.Algorithm = 'interior-point';
	prbm.minopts.SpecifyObjectiveGradient = true;
	% prbm.minopts.OptimalityTolerance = 1e-5;
	prbm.minopts.ConstraintTolerance = 1e-4;
	prbm.minopts.MaxIterations = 50;
% 	prbm.minopts.CheckGradients = true;
% 	prbm.minopts.FiniteDifferenceType = 'central';
	% prbm.minopts.InitTrustRegionRadius = sqrt(sum(prbm.num_vars)/100); [little effect]
	% prbm.minopts.ScaleProblem = 'obj-and-constr'; [even slower]
	% prbm.minopts.SubproblemAlgorithm = 'cg'; % only take cg steps. [even slower]
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
	prbm.minopts.Display = 'iter';   %{'notify-detailed'|'notify'|'iter'}
end
options = structupdate(options, prbm, {'bCompact', 'unit'});
setdefault(options, pardata, {'PricingPolicy'}, 'warning');
options.bFinal = false;
if strcmpi(this.options.OptimizationTool, 'matlab')
	if strcmpi(prbm.minopts.Algorithm, 'interior-point')
		prbm.minopts.HessianFcn = ...
			@(x,lbd)NormalSliceOptimizer.fcnHessianCC(x, lbd, this, options);
	end
	if isfield(options, 'Warning') && strcmpi(options.Warning, 'off')
		warning('off');
	end
	x0 = prbm.x0/prbm.unit;
	[xs, fval, exitflag, foutput] = ...
		fmincon(@(x)NormalSliceOptimizer.fcnProfitPrimalCC(x, this, options), ...
		x0, prbm.A, prbm.b, prbm.Aeq, prbm.beq, prbm.lb, [], [], prbm.minopts);
	if isfield(options, 'Warning') && strcmpi(options.Warning, 'off')
		warning('on');
	end
	this.interpretExitflag(exitflag, foutput);
	xs = xs*options.unit;
else
	error('error: Un-recognized optimization tool.');
end

%% output solution
% assert(this.checkFeasible())
% 			tol_zero = this.options.NonzeroTolerance;
% 			xs(xs<tol_zero*mean(xs)) = 0;
x = sparse(num_vars, 1);
x(this.I_active_variables) = xs;
this.x0 = x;
this.temp_vars.x = x(1:num_varx);
this.temp_vars.z = x(num_varx+(1:num_varz));
this.temp_vars.r = x(num_varx+num_varz+(1:Nf));
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
output.funcCount = foutput.funcCount;
end