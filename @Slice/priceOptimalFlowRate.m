%% priceOptimalFlowRate
% priceOptimalFlowRate return the optimal flow rate for each flow in the slice, and the
% net profit of the slice. 
% TODO: by setting the capacity as |inf|, this method is equivalent to _optimalFlowRate_.
% NOTE: prices should be updated before calling this function.

%% Function Prototype
%   [net_profit, node_load, link_load] = priceOptimalFlowRate(this, x0, options)
% |options|: if price is not provided in |options|, then this method will use the price
% stored in the slice.
function [net_profit, node_load, link_load] = priceOptimalFlowRate(this, x0, options)
global InfoLevel;
options = structmerge(...
    getstructfields(options, 'PricingPolicy', 'empty-ignore'),...
    getstructfields(this.Parent.options, 'Form', 'empty-ignore'));
if isempty(options.Form)
    options.Form = 'normal';
    warning('Slice.calculateOutput: <Form> set to default (%s).', options.Form);
end
if isempty(options.PricingPolicy)
    options.PricingPolicy = 'linear';
    warning('Slice.calculateOutput: <PricingPolicy> set to default (%s).', options.PricingPolicy);
end

NP = this.NumberPaths;
NV = this.NumberVNFs;

if nargin == 2 && ~isempty(x0)
    this.x0 = x0;
else%if isempty(this.x0)
    this.x0 = zeros(this.num_vars,1);
    this.x0(1:NP) = 1;
    alpha_max = max(this.Parent.VNFTable.ProcessEfficiency(this.VNFList));
    this.x0((NP+1):end) = alpha_max;
end
assert(this.checkFeasible(this.x0), 'error: infeasible start point.');

bs = sparse(this.num_lcon_res,1);

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
fmincon_opt.Display = InfoLevel.InnerModel.char;   %'notify-detailed'; %'iter';
if isfield(options, 'Form') && strcmpi(options.Form, 'compact')
    %     isequal(this.I_active_variables', sum(this.As_res,1)~=0)
    z_filter = sparse(repmat(...
        reshape(logical(this.I_node_path), numel(this.I_node_path),1), NV, 1));
    this.I_active_variables = [true(NP,1) ;  z_filter];
    As_compact = this.As_res(:, this.I_active_variables);
    var0_compact = this.x0(this.I_active_variables);
    lbs = sparse(length(var0_compact),1);
    options.num_orig_vars = this.num_vars;
    fmincon_opt.HessianFcn = ...
        @(x,lambda)Slice.fcnHessianCompact(x, lambda, this, options);
    [x_compact, fval, exitflag] = ...
        fmincon(@(x)Slice.fcnProfitCompact(x,this, options), ...
        var0_compact, As_compact, bs, [], [], lbs, [], [], fmincon_opt);
    x = zeros(this.num_vars, 1);
    x(this.I_active_variables) = x_compact;
else
    lbs = sparse(this.num_vars,1);
    fmincon_opt.HessianFcn = @(x,lambda)Slice.fcnHessian(x, lambda, this, options);
    [x, fval, exitflag] = fmincon(@(x)Slice.fcnProfit(x, this, options), ...
        this.x0, this.As_res, bs, [], [], lbs, [], [], fmincon_opt);
end
if InfoLevel.UserModelDebug == DisplayLevel.Iteration
    fprintf('\tThe optimal net profit of the slice: %G.\n', -fval);
end
if exitflag == 0
    if InfoLevel.UserModel == DisplayLevel.Notify
        warning('reaching maximum number of iterations.');
    end
elseif exitflag < 0
    error('abnormal exit with flag %d.',exitflag);
elseif exitflag ~= 1
    if InfoLevel.UserModel == DisplayLevel.Notify
        warning('local optimal solution found.');
    end
end

%% Output solution
assert(this.checkFeasible(x, struct('ConstraintTolerance', fmincon_opt.ConstraintTolerance)),...
    'error: infeasible solution.');
this.temp_vars.x = x(1:this.NumberPaths);
this.temp_vars.z = x((this.NumberPaths+1):end);
%%%
% When compute node load, z_npf corresponding to h_np = 0 has been set as zero.
nz = this.NumberDataCenters*this.NumberPaths;
z_index = 1:nz;
for f = 1:this.NumberVNFs
    this.temp_vars.z(z_index) = this.I_node_path(:).*this.temp_vars.z(z_index);
    z_index = z_index + nz;
end
% tol_zero = this.Parent.options.NonzeroTolerance;
% this.temp_vars.x(this.temp_vars.x<tol_zero*max(this.temp_vars.x)) = 0;
% this.temp_vars.z(this.temp_vars.z<tol_zero*max(this.temp_vars.z)) = 0;
% if ~this.checkFeasible([this.temp_vars.x; this.temp_vars.z])
%     if InfoLevel.UserModel >= DisplayLevel.Notify
%         warning('priceOptimalFlowRate: the rounding of variables %s', ...
%             'with small quantity will make the solution infeasible.');
%     end
% end
this.x0 = x;

if nargout >= 2
    node_load = zeros(this.Parent.NumberDataCenters,1);
    data_center_id = this.getDCPI;
    node_load(data_center_id) = this.getNodeLoad(this.temp_vars.z);
end
if nargout >= 3
    link_load = zeros(this.Parent.NumberLinks,1);
    link_load(this.VirtualLinks.PhysicalLink) = this.getLinkLoad(this.temp_vars.x);    
end
this.flow_rate = this.getFlowRate(this.temp_vars.x);
net_profit = -fval;
%%%
% FOR DEBUG
% this.setPathBandwidth(this.temp_vars.x);
end