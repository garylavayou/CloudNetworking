%% priceOptimalFlowRate
% Extending <Slice.priceOptimalFlowRate>, Considering resource reconfiguration cost.
% NOTE: update price before call this function.
function [net_profit, node_load, link_load] = priceOptimalFlowRate(this, x0, new_opts)
global DEBUG INFO;
if nargin <= 2
    new_opts = struct;
end
if strcmpi(this.options.ReconfigMethod, 'dimconfig0') ||... 
        (isempty(fieldnames(this.net_changes)) && (false == this.b_dim))
    if nargout <= 1
        net_profit = priceOptimalFlowRate@Slice(this, x0, new_opts);
    else
        [net_profit, node_load, link_load] = priceOptimalFlowRate@Slice(this, x0, new_opts);
    end
    return;
end

options = structmerge(new_opts, ...
    getstructfields(this.Parent.options, 'Form'), ...
    getstructfields(this.options, 'PricingPolicy', 'default', 'linear'),...
    'exclude');     

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
%   (6) virtual link capacity constraint: NL (optional);
%   (7) lower bound of virtual link capacity (optional);
%   (8) lower bound of virtual node capacity (dimconfig2): sum capacity of VNF isntance
%       larger than the lowbound;
%   (9) lower bound of VNF instance capacity (dimconfig1).
%
% Number of Variables
%   x: path; z: VNF instances assignment; v: VNF instances capacity;
%   tx;      tx;                          tv;     
%   c: link capacity (optional)
num_lcon = this.num_lcon_res + this.num_varv + ...
    2*NP + 2*this.num_varz + 2*this.num_varv;
nnz_As = nnz(As_res) + (nnz(Hd)+this.num_varv) + ...
    + 4*NP + 4*this.num_varz+ 4*this.num_varv;
num_base_vars = this.num_vars + this.num_varv;
num_vars = 2*num_base_vars;
if ~isempty(this.lower_bounds)
    NL = this.NumberVirtualLinks;
    num_lcon0 = num_lcon;
    num_lcon = num_lcon + NL;
    num_vars = num_vars + NL;
    if isfield(this.lower_bounds, 'node')
        NC = this.NumberDataCenters;
        num_lcon = num_lcon + NC;
    end
end
As = spalloc(num_lcon, num_vars, nnz_As);
%% (1)
As(1:this.num_lcon_res,1:this.num_vars) = As_res;
row_offset = this.num_lcon_res;
%% (2)
As(row_offset+(1:this.num_varv), NP+(1:this.num_varz)) = Hd;
As(row_offset+(1:this.num_varv), this.num_vars+(1:this.num_varv)) = -eye(this.num_varv);
row_offset = row_offset + this.num_varv;
%% (3)
As(row_offset+(1:NP), 1:NP) = eye(NP);
As(row_offset+(1:NP), num_base_vars+(1:NP)) = -eye(NP);
row_offset = row_offset + NP;
As(row_offset+(1:NP), 1:NP) = -eye(NP);
As(row_offset+(1:NP), num_base_vars+(1:NP)) = -eye(NP);
row_offset = row_offset + NP;
%% (4)
As(row_offset+(1:this.num_varz), (NP+1):this.num_vars) = eye(this.num_varz);
As(row_offset+(1:this.num_varz), (num_base_vars+NP)+(1:this.num_varz)) = -eye(this.num_varz);
row_offset = row_offset + this.num_varz;
As(row_offset+(1:this.num_varz), (NP+1):this.num_vars) = -eye(this.num_varz);
As(row_offset+(1:this.num_varz), (num_base_vars+NP)+(1:this.num_varz)) = -eye(this.num_varz);
row_offset = row_offset + this.num_varz;
%% (5)
As(row_offset+(1:this.num_varv), (this.num_vars+1):num_base_vars) = eye(this.num_varv);
As(row_offset+(1:this.num_varv), (num_base_vars+this.num_vars+1):num_base_vars*2) = -eye(this.num_varv);
row_offset = row_offset + this.num_varv;
As(row_offset+(1:this.num_varv), (this.num_vars+1):num_base_vars) = -eye(this.num_varv);
As(row_offset+(1:this.num_varv), (num_base_vars+this.num_vars+1):num_base_vars*2) = -eye(this.num_varv);
%% (6)
row_offset = row_offset + this.num_varv;
if ~isempty(this.lower_bounds)
    As(row_offset+(1:NL), 1:NP) = this.I_edge_path;
    As(row_offset+(1:NL), 2*num_base_vars+(1:NL)) = -eye(NL);
    row_offset = row_offset + NL;
    if isfield(this.lower_bounds, 'node')
        As(row_offset+(1:NC), this.num_vars+(1:this.num_varv)) = -repmat(eye(NC),1,NV);
    end    
end

bs = [sparse(this.num_lcon_res+this.num_varv,1);
    this.topts.old_variables_x;       % which have been pre-processed, so it can be
    -this.topts.old_variables_x;      % compared with the current states.
    this.topts.old_variables_z;
    -this.topts.old_variables_z;
    this.old_variables.v;
    -this.old_variables.v];
if ~isempty(this.lower_bounds)
    bs = [bs; zeros(NL,1)];
    if isfield(this.lower_bounds, 'node')
        bs = [bs; -this.lower_bounds.node];
    end
end

lbs = sparse(num_vars, 1);
if ~isempty(this.lower_bounds)
    lbs(2*num_base_vars+(1:NL)) = this.lower_bounds.link;
    if isfield(this.lower_bounds, 'VNF')
        lbs(this.num_vars+(1:this.num_varv)) = this.lower_bounds.VNF;
    end
end

var0 = [this.topts.old_variables_x;
    this.topts.old_variables_z;
    this.old_variables.v;
    this.topts.old_variables_x;
    this.topts.old_variables_z;
    this.old_variables.v];
if ~isempty(this.lower_bounds)
    var0 = [var0; this.old_state.link_capacity];
end
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
fmincon_opt.MaxIterations = 60;
fmincon_opt.MaxFunctionEvaluations = 180;
fmincon_opt.Display = 'notify';   %'notify-detailed'; %'iter'; 'notify'
% fmincon_opt.CheckGradients = true;
% fmincon_opt.FiniteDifferenceType = 'central';
% fmincon_opt.Diagnostics = 'on';
options.num_varx = this.NumberPaths;
options.num_varz = this.num_varz;
options.num_varv = this.num_varv;
% options.Form = 'normal';
if isfield(options, 'Form') && strcmpi(options.Form, 'compact')
    % column/variables: isequal(this.I_active_variables', sum(this.As_res,1)~=0)
    % row/constraints:  isequal(active_rows, find(sum(this.As_res,2)~=0))  
    z_filter = sparse(repmat(...
        reshape(logical(this.I_node_path), numel(this.I_node_path),1), NV, 1));
    this.I_active_variables = [true(NP,1);  z_filter;  true(this.num_varv,1);...
        true(NP,1); z_filter; true(this.num_varv,1)];
    row_offset = this.num_lcon_res + this.num_varv + 2*NP;
    active_rows = [(1:row_offset)'; row_offset+find(z_filter); ...
        row_offset+this.num_varz+find(z_filter)];
    if ~isempty(this.lower_bounds)
        this.I_active_variables = [this.I_active_variables; true(NL,1)];
        active_rows = [active_rows; ((num_lcon0-2*this.num_varv+1):num_lcon0)';...
            num_lcon0+(1:NL)'];
        if isfield(this.lower_bounds, 'node')
            active_rows = [active_rows; num_lcon0+NL+(1:NC)'];
        end
    else
        active_rows = [active_rows; ((num_lcon-2*this.num_varv+1):num_lcon)'];
    end
    As_compact = As(active_rows, this.I_active_variables);
    var0_compact = var0(this.I_active_variables);
    lbs = lbs(this.I_active_variables);
    bs = bs(active_rows);
    options.num_orig_vars = this.num_vars+this.num_varv;
    fmincon_opt.HessianFcn = ...
        @(x,lambda)DynamicSlice.fcnHessianCompact(x, lambda, this, options);
    [x_compact, fval, exitflag, output] = ...
        fmincon(@(x)DynamicSlice.fcnProfitCompact(x,this, options), ...
        var0_compact, As_compact, bs, [], [], lbs, [], [], fmincon_opt);
    x = zeros(num_vars, 1);
    x(this.I_active_variables) = x_compact;
else
    fmincon_opt.HessianFcn = @(x,lambda)DynamicSlice.fcnHessian(x, lambda, this, options);
    [x, fval, exitflag, output] = fmincon(@(x)DynamicSlice.fcnProfit(x, this, options), ...
        var0, As, bs, [], [], lbs, [], [], fmincon_opt);
end
this.interpretExitflag(exitflag, output.message);
if (~isempty(DEBUG) && DEBUG) || (~isempty(INFO) && INFO)
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

% Passing out the load information is not necessary, since the network only reuqire the
% information of occupied capacity.
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