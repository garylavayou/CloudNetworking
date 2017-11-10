%% priceOptimalFlowRate
% Extending <Slice.priceOptimalFlowRate>, Considering resource reconfiguration cost.
% NOTE: update price before call this function.
function [net_profit, node_load, link_load] = priceOptimalFlowRate(this, x0, options)
if isempty(fieldnames(this.net_changes)) && (false == this.b_dim)
    if nargout <= 1
        net_profit = priceOptimalFlowRate@Slice(this, x0, options);
    else
        [net_profit, node_load, link_load] = priceOptimalFlowRate@Slice(this, x0, options);
    end
    return;
end
options = structmerge(...
    getstructfields(options, 'PricingPolicy', 'empty-ignore'),...
    getstructfields(this.Parent.options, 'Form', 'empty-ignore'));
if isempty(options.Form)
    options.Form = 'normal';
    warning('Slice.calculateOutput: <Form> set to default (%s).', options.Form);
end
if isempty(options.PricingPolicy)
    options.PricingPolicy = 'linear';
    warning('Slice.calculateOutput: <PricingPolicy> set to default (%s).', ...
        options.PricingPolicy);
end
options.Method = this.options.Method;

global InfoLevel;
NP = this.NumberPaths;
NV = this.NumberVNFs;

Hd = this.Hdiag;
As_res = this.As_res;        % update As_res

%%
% List of constaints (see also <DynamicSlice.fastReconfigure2>):
%   (1) flow processing requirement: NP*NV (this.num_lcon_res);
%   (2) VNF instance capacity constraint: NV*NN (this.num_varv);
%   (3) link reconfiguration cost constraint: 2*NP;
%   (4) node reconfiguration cost constraint: 2*NN*NP*NV (2*this.num_varz);
%   (5) VNF instance reconfiguration cost constraint: 2*NV*NN;
% Number of Variables
%   x: path; z: VNF instances assignment; v: VNF instances capacity;
%   tx;      tx;                          tv;     
num_lcon = this.num_lcon_res + this.num_varv + ...
    2*NP + 2*this.num_varz + 2*this.num_varv;
nnz_As = nnz(As_res) + (nnz(Hd)+this.num_varv) + ...
    + 4*NP + 4*this.num_varz+ 4*this.num_varv;
num_vars = 2*this.num_vars + 2*this.num_varv;
As = spalloc(num_lcon, num_vars, nnz_As);
As(1:this.num_lcon_res,1:this.num_vars) = As_res;
row_offset = this.num_lcon_res;
As(row_offset+(1:this.num_varv), NP+(1:this.num_varz)) = Hd;
As(row_offset+(1:this.num_varv), this.num_vars+(1:this.num_varv)) = -eye(this.num_varv);
row_offset = row_offset + this.num_varv;
As(row_offset+(1:NP), 1:NP) = eye(NP);
As(row_offset+(1:NP), num_vars/2+(1:NP)) = -eye(NP);
row_offset = row_offset + NP;
As(row_offset+(1:NP), 1:NP) = -eye(NP);
As(row_offset+(1:NP), num_vars/2+(1:NP)) = -eye(NP);
row_offset = row_offset + NP;
As(row_offset+(1:this.num_varz), (NP+1):this.num_vars) = eye(this.num_varz);
As(row_offset+(1:this.num_varz), (num_vars/2+NP)+(1:this.num_varz)) = -eye(this.num_varz);
row_offset = row_offset + this.num_varz;
As(row_offset+(1:this.num_varz), (NP+1):this.num_vars) = -eye(this.num_varz);
As(row_offset+(1:this.num_varz), (num_vars/2+NP)+(1:this.num_varz)) = -eye(this.num_varz);
row_offset = row_offset + this.num_varz;
As(row_offset+(1:this.num_varv), (this.num_vars+1):num_vars/2) = eye(this.num_varv);
As(row_offset+(1:this.num_varv), (num_vars/2+this.num_vars+1):end) = -eye(this.num_varv);
row_offset = row_offset + this.num_varv;
As(row_offset+(1:this.num_varv), (this.num_vars+1):num_vars/2) = -eye(this.num_varv);
As(row_offset+(1:this.num_varv), (num_vars/2+this.num_vars+1):end) = -eye(this.num_varv);

bs = [sparse(this.num_lcon_res+this.num_varv,1);
    this.topts.old_variables_x;       % which have been pre-processed, so it can be
    -this.topts.old_variables_x;      % compared with the current states.
    this.topts.old_variables_z;
    -this.topts.old_variables_z;
    this.old_variables.v;
    -this.old_variables.v];

var0 = [this.topts.old_variables_x;
    this.topts.old_variables_z;
    this.old_variables.v;
    this.topts.old_variables_x;
    this.topts.old_variables_z;
    this.old_variables.v];
assert(this.checkFeasible(var0), 'error: infeasible start point.');

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
fmincon_opt.Display = InfoLevel.InnerModelDebug.char;   %'notify-detailed'; %'iter';
% options.Form = 'normal';
options.num_varx = this.NumberPaths;
options.num_varz = this.num_varz;
options.num_varv = this.num_varv;
if isfield(options, 'Form') && strcmpi(options.Form, 'compact')
    % column/variables: isequal(this.I_active_variables', sum(this.As_res,1)~=0)
    % row/constraints:  isequal(active_rows, find(sum(this.As_res,2)~=0))  
    z_filter = sparse(repmat(...
        reshape(logical(this.I_node_path), numel(this.I_node_path),1), NV, 1));
    this.I_active_variables = [true(NP,1);  z_filter;  true(this.num_varv,1);...
        true(NP,1); z_filter; true(this.num_varv,1)];
    row_offset = this.num_lcon_res + this.num_varv + 2*NP;
    active_rows = [(1:row_offset)'; row_offset+find(z_filter); ...
        row_offset+this.num_varz+find(z_filter); ...
        ((num_lcon-2*this.num_varv+1):num_lcon)'];
    As_compact = As(active_rows, this.I_active_variables);
    var0_compact = var0(this.I_active_variables);
    bs = bs(active_rows);
    lbs = sparse(length(var0_compact),1);
%     old_z_reconfig_cost = this.z_reconfig_cost;
%     this.z_reconfig_cost = this.z_reconfig_cost(z_filter);
    options.num_orig_vars = this.num_vars+this.num_varv;
    fmincon_opt.HessianFcn = ...
        @(x,lambda)DynamicSlice.fcnHessianCompact(x, lambda, this, options);
    [x_compact, fval, exitflag, output] = ...
        fmincon(@(x)DynamicSlice.fcnProfitCompact(x,this, options), ...
        var0_compact, As_compact, bs, [], [], lbs, [], [], fmincon_opt);
    x = zeros(num_vars, 1);
    x(this.I_active_variables) = x_compact;
%     this.z_reconfig_cost = old_z_reconfig_cost;     % recover.
else
    lbs = sparse(num_vars,1);
    fmincon_opt.HessianFcn = @(x,lambda)DynamicSlice.fcnHessian(x, lambda, this, options);
    [x, fval, exitflag, output] = fmincon(@(x)DynamicSlice.fcnProfit(x, this, options), ...
        var0, As, bs, [], [], lbs, [], [], fmincon_opt);
end
this.interpretExitflag(exitflag, output.message);
if InfoLevel.UserModelDebug >= DisplayLevel.Iteration
    fprintf('\tThe optimal net profit of the slice: %G.\n', -fval);
end

%% Output solution
assert(this.checkFeasible(x, struct('ConstraintTolerance', fmincon_opt.ConstraintTolerance)), ...
    'error: infeasible solution.');
this.x0 = x;
this.convert(x);
%%%
% The post processing like <optimalFlowRate> is not needed, since the
% objective function (reconfiguration cost) will force those variables to be zero.

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