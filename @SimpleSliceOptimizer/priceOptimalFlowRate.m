%% priceOptimalFlowRate
% priceOptimalFlowRate return the optimal flow rate for each flow in the slice, and the
% net profit of the slice. 
% TODO: by setting the capacity as |inf|, this method is equivalent to _optimalFlowRate_.
% NOTE: prices should be updated before calling this function.

%% Function Prototype
%   [output, loads, fval] = priceOptimalFlowRate(this, x0, options)
% |options|: if price is not provided in |options|, then this method will use the price
% stored in the slice.
function output = priceOptimalFlowRate(this, x0, options)
global DEBUG INFO;
if nargin <= 2
	options = Dictionary();
else
	options = Dictionary(options);
end
setdefault(options, struct('bParallel', false, 'bInitialize', false));
pardata = this.pardata;
Nsn = pardata.NumberServiceNodes;
Np = pardata.NumberPaths;
Nvnf = pardata.NumberVNFs;

if options.bInitialize
	this.num_vars = [Np; Nvnf*Nsn*Np];
	num_vars = sum(this.num_vars(1:2));
	num_varx = this.num_vars(1);
	num_varz = this.num_vars(2);
	Nl = pardata.NumberLinks;
	num_lcon_res = size(this.As_res,1);
	
	this.problem = Dictionary(); prbm = this.problem;
	prbm.x0 = zeros(num_vars,1);
	prbm.x0(1:num_varx) = 1;
	alpha_max = max(pardata.ProcessEfficiency);
	prbm.x0((num_varx+1):end) = alpha_max;
	assert(this.checkFeasible(prbm.x0), 'error: infeasible start point.');	
	prbm.A = this.As_res;
	prbm.b = sparse(num_lcon_res,1);
	if isfield(options, 'CapacityConstrained') && options.CapacityConstrained
		% See <SimpleSliceOptimizer.optimalFlowRate> for the capacity constraints.
		prbm.A = [prbm.A; this.I_edge_path, sparse(Nl, num_varz); ...
			sparse(Nsn, Np); this.Hrep];
		prbm.b = [prbm.b; this.capacities.Link; this.capacities.Node];
	end
	
	%% Compress
	% options.Form = 'normal';
	if strcmpi(this.options.Form, 'compact')
		%     isequal(this.I_active_variables', sum(this.As_res,1)~=0)
		z_filter = sparse(repmat(...
			reshape(logical(this.I_dc_path), numel(this.I_dc_path),1), Nvnf, 1));
		this.I_active_variables = [true(num_varx,1) ;  z_filter];
		prbm.A = prbm.A(:, this.I_active_variables);
		prbm.x0 = prbm.x0(this.I_active_variables);
		prbm.lb = sparse(length(prbm.x0),1);
		prbm.bCompact = true;
	else
		prbm.lb = sparse(num_vars,1);
		prbm.bCompact = false;
	end
	%% Set the optimization options
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
	% diagnostics
	prbm.minopts.Display = 'iter';   % {'notify-detailed'|'iter'|'notify'};
	% minopts.Diagnostics = 'on';
	% minopts.CheckGradients = true;
	% minopts.FiniteDifferenceType = 'central';
	% minopts.FiniteDifferenceStepSize = 1e-10;
	%% Save information
	if options.bParallel
		prbm.num_full_vars = this.num_vars;
		if prbm.bCompact
			prbm.I_active_variables = this.I_active_variables; % update to optimizer
		end
	end
else
	prbm = this.problem;
	num_vars = sum(this.num_vars(1:2));
	num_varx = this.num_vars(1);
	num_varz = this.num_vars(2);
end
options = structupdate(options, prbm, {'bCompact'}, 'ignore');
if options.bCompact
	options.num_orig_vars = sum(this.num_vars);
end
options = setdefault(options, this.pardata, {'PricingPolicy'});
options.bFinal = false;

if nargin >= 2 && ~isempty(x0)
	assert(this.checkFeasible(x0, options), 'error: infeasible start point.');
	if options.bCompact
		prbm.x0 = x0(this.I_active_variables);
	else
		prbm.x0 = x0;
	end
end

prbm.minopts.HessianFcn = ...
    @(x,lambda)SimpleSliceOptimizer.fcnHessian(x, lambda, this, options);
[xs, fval, exitflag, foutput] = fmincon(@(x)SimpleSliceOptimizer.fcnProfit(x, this, options), ...
    prbm.x0, prbm.A, prbm.b, [], [], prbm.lb, [], [], prbm.minopts);
if options.bCompact
	x = sparse(num_vars, 1);
	x(this.I_active_variables) = xs;
else
	x = xs;
end
SimpleSliceOptimizer.interpretExitflag(exitflag, foutput.message);
if (~isempty(DEBUG) && DEBUG) || (~isempty(INFO) && INFO)
    fprintf('\tThe optimal net profit of the slice: %G.\n', -fval);
end
% tol_zero = this.Parent.options.NonzeroTolerance;
% this.temp_vars.x(this.temp_vars.x<tol_zero*max(this.temp_vars.x)) = 0;
% this.temp_vars.z(this.temp_vars.z<tol_zero*max(this.temp_vars.z)) = 0;
% if ~this.checkFeasible([this.temp_vars.x; this.temp_vars.z])
%         warning('priceOptimalFlowRate: the rounding of variables %s', ...
%             'with small quantity will make the solution infeasible.');
% end

%% Output solution
options.ConstraintTolerance = prbm.minopts.ConstraintTolerance;
assert(this.checkFeasible(x, options), 'error: infeasible solution.');
this.x0 = x;
this.temp_vars.x = x(1:num_varx);
this.temp_vars.z = x((num_varx+1):(num_varx+num_varz));
nz = Nsn*Np;
z_index = 1:nz;
for f = 1:Nvnf
	% When compute node load, z_Npf corresponding to h_Np = 0 has been set as zero.
	this.temp_vars.z(z_index) = this.I_dc_path(:).*this.temp_vars.z(z_index);
	z_index = z_index + nz;
end
output = Dictionary();
output.net_profit = -fval;
output.iterations = foutput.iterations;
output.funcCount = foutput.funcCount;
if options.bParallel
	output.x0 = this.x0;
	output.temp_vars = this.temp_vars;
	if options.bInitialize
		output.problem = prbm;
	end
end
end
