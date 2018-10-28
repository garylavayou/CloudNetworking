%% priceOptimalFlowRate
% priceOptimalFlowRate return the optimal flow rate for each flow in the slice, and the
% net profit of the slice. 
% TODO: by setting the capacity as |inf|, this method is equivalent to _optimalFlowRate_.
% NOTE: prices should be updated before calling this function.

%% Function Prototype
%   [output, loads, fval] = priceOptimalFlowRate(this, x0, options)
% |options|: if price is not provided in |options|, then this method will use the price
% stored in the slice.
function [output, loads, fval] = priceOptimalFlowRate(this, x0, options)
global DEBUG INFO;
options = structmerge(...
    getstructfields(options, 'PricingPolicy', 'default', {'quadratic'}),...
    getstructfields(this.Parent.options, 'Form', 'default', {'normal'}));

Np = this.hs.NumberPaths;
Nvf = this.hs.NumberVNFs;
Ne = this.hs.NumberLinks;
Nsn = this.hs.NumberServiceNodes;

if nargin >= 2 && ~isempty(x0)
    this.x0 = x0;
else%if isempty(this.x0)
    this.x0 = zeros(this.NumberVariables,1);
    this.x0(1:Np) = 1;
    alpha_max = max(this.Parent.VNFTable.ProcessEfficiency(this.VNFList));
    this.x0((Np+1):end) = alpha_max;
end
num_vars = length(this.x0);
assert(this.checkFeasible(this.x0), 'error: infeasible start point.');

bs = sparse(this.NumberLinearConstraints,1);

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
fmincon_opt = optimoptions(@fmincon);
fmincon_opt.Algorithm = 'interior-point';
fmincon_opt.SpecifyObjectiveGradient = true;
%% diagnostics
fmincon_opt.Display = 'notify';   %'notify-detailed'; %'iter';
% fmincon_opt.CheckGradients = true;
% fmincon_opt.FiniteDifferenceType = 'central';
% fmincon_opt.FiniteDifferenceStepSize = 1e-10;
% fmincon_opt.Diagnostics = 'on';
%%
% options.Form = 'normal';
if strcmpi(options.Form, 'compact')
    %     isequal(this.I_active_variables', sum(this.As_res,1)~=0)
    z_filter = sparse(repmat(...
        reshape(logical(this.I_dc_path), numel(this.I_dc_path),1), Nvf, 1));
    this.I_active_variables = [true(Np,1) ;  z_filter];
    As = this.As_res(:, this.I_active_variables);
    var0 = this.x0(this.I_active_variables);
    lbs = sparse(length(var0),1);
    options.num_orig_vars = this.NumberVariables;
    options.bCompact = true;
else
    lbs = sparse(this.NumberVariables,1);
    As = this.As_res;
    var0 = this.x0;
end
if isfield(options, 'CapacityConstrained') && options.CapacityConstrained
    % See <SimpleSliceOptimizer.optimalFlowRate> for the capacity constraints.
    As = [As; this.I_edge_path, sparse(Ne, this.num_varz); ...
    sparse(Nsn, Np); this.Hrep];
    bs = [bs; this.capacities.Link; this.capacities.Node];
end
fmincon_opt.HessianFcn = ...
    @(x,lambda)SimpleSlice.fcnHessian(x, lambda, this, options);
[xs, fval, exitflag, foutput] = fmincon(@(x)SimpleSlice.fcNprofit(x, this, options), ...
    var0, As, bs, [], [], lbs, [], [], fmincon_opt);
if strcmpi(options.Form, 'compact')
    x = zeros(num_vars, 1);
    x(this.I_active_variables) = xs;
else
    x = xs;
end
this.interpretExitflag(exitflag, foutput.message);
if (~isempty(DEBUG) && DEBUG) || (~isempty(INFO) && INFO)
    fprintf('\tThe optimal net profit of the slice: %G.\n', -fval);
end

%% Output solution
assert(this.checkFeasible(x, struct('ConstraintTolerance', fmincon_opt.ConstraintTolerance)),...
    'error: infeasible solution.');
output.temp_vars.x = x(1:this.NumberPaths);
output.temp_vars.z = x((this.NumberPaths+1):end);
%%%
% When compute node load, z_Npf corresponding to h_Np = 0 has been set as zero.
nz = Nsn*Np;
z_index = 1:nz;
for f = 1:Vnf
    this.temp_vars.z(z_index) = this.I_dc_path(:).*this.temp_vars.z(z_index);
    z_index = z_index + nz;
end
% tol_zero = this.Parent.options.NonzeroTolerance;
% this.temp_vars.x(this.temp_vars.x<tol_zero*max(this.temp_vars.x)) = 0;
% this.temp_vars.z(this.temp_vars.z<tol_zero*max(this.temp_vars.z)) = 0;
% if ~this.checkFeasible([this.temp_vars.x; this.temp_vars.z])
%         warning('priceOptimalFlowRate: the rounding of variables %s', ...
%             'with small quantity will make the solution infeasible.');
% end
output.x0 = x;

if nargout >= 2
    loads = zeros(this.hs.Parent.NumberLinks+this.hs.Parent.NumberDataCenters, 1);
    idx = [this.hs.Links.PhysicalLink; this.Parent.NumberLinks+this.getDCPI];
    loads(idx) = [this.hs.getLinkLoad(this.temp_vars.x); this.hs.getNodeLoad(this.temp_vars.z)];
end
output.net_profit = -fval;
output.iterations = foutput.iterations;
output.funcCount = foutput.funcCount;
output.flow_rate = this.getFlowRate(this.temp_vars.x);
%%%
% FOR DEBUG
% this.setPathBandwidth(this.temp_vars.x);
end
