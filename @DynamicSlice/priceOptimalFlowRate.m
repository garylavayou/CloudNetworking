%% priceOptimalFlowRate
% Extending <Slice.priceOptimalFlowRate>, Considering resource reconfiguration cost.
% NOTE: update price before call this function.
function [net_profit, node_load, link_load] = priceOptimalFlowRate(this, x0, new_opts)
global DEBUG INFO;
if nargin <= 2
    new_opts = struct;
end
if strcmpi(this.options.ReconfigMethod, 'dimconfig0') ||... 
        (isempty(fieldnames(this.net_changes)) && (false == this.b_dim)) &&...
        (~isfield(this.options, 'Reserve') || isempty(this.options.Reserve))
    if nargout <= 1
        net_profit = priceOptimalFlowRate@Slice(this, x0, new_opts);
    else
        [net_profit, node_load, link_load] = priceOptimalFlowRate@Slice(this, x0, new_opts);
    end
    return;
end

options = structmerge(new_opts, ...
    getstructfields(this.Parent.options, 'Form', 'default', 'normal'), ...
    getstructfields(this.options, {'PricingPolicy','Reserve'}, 'default', {'quadratic', 1}),...
    'exclude');     

NP = this.NumberPaths;
NV = this.NumberVNFs;
NL = this.NumberVirtualLinks;
NC = this.NumberDataCenters;
As_res = this.As_res;        % update As_res

%% List of constaints
% (see also <DynamicSlice.fastReconfigure2>):
%   (1) flow processing requirement: NP*NV (this.num_lcon_res);
%   (2) VNF instance capacity constraint: NV*NN (this.num_varv); virtual node capacity constraint;
%   (3) link reconfiguration cost constraint: 2*NP;
%   (4) VNF allocation reconfiguration cost constraint: 2*NN*NP*NV (2*this.num_varz);
%   (5) VNF instance capacity reconfiguration cost constraint: 2*NV*NN;
%   (6) virtual link capacity constraint: NL (optional, resource reservation);
%   (7) lower bound of virtual link capacity (optional);
%   (8) lower bound of virtual node capacity (dimconfig2): sum capacity of VNF isntance
%       larger than the lowbound; 
%   (9) lower bound of VNF instance capacity (dimconfig1);
%   (10)resource utilization lower bound for resource reservation (optional);
%   (11)node capacity constraint with resource reservation (optional)
%
%           1	2	3	4	5	6	7	8	9	10 
% InitRSV   +   +   -   -   -   +   -   -   -   +   
% Redim     +   +   +   +   +   +   -   -   -   -   
% RedimRSV  +   +   +   +   +   +   +  (+) (+)  -  
% RedimRSV+ +   +   +   +   +   +   +  (+) (+)  +  
%% Number of Variables
%   x: path; z: VNF instances assignment; v: VNF instances capacity;
%   tx;      tz;                          tv;     
%   c: link capacity (optional)
%           x    z    v 	tx    tz     tv     c
% InitRSV   +    +    + 	-     -      -      +
% Redim     +    +    + 	+     +      +      +
% RedimRSV  +    +    +     +     +      +      +
% RedimRSV+ +    +    +     +     +      +      +
if this.b_derive_vnf    % && strcmpi(options.Reserve, 'err')
    this.invoke_method = 1;  % initialization
elseif isempty(this.lower_bounds)
    this.invoke_method = 2;  % no resource reservation
elseif options.Reserve >= 1  
    this.invoke_method = 3;  % implicit resource reservation
elseif options.Reserve > 0
    this.invoke_method = 4;  % explicit resource reservation
else 
    error('%s: invalid value ''Reserve=%.2f''for resource reservation.', calledby, options.Reserve);
end
num_base_vars = this.num_vars + this.num_varv;
num_vars = num_base_vars;
num_lcon = this.num_lcon_res + this.num_varv;
nnz_As = nnz(As_res) + (nnz(this.Hdiag)+this.num_varv);
if this.invoke_method ~= 1
    num_vars = num_vars + num_base_vars;
    num_lcon = num_lcon + ...
        2*NP + 2*this.num_varz + 2*this.num_varv;
    nnz_As = nnz_As + 4*NP + 4*this.num_varz+ 4*this.num_varv;
end
num_vars = num_vars + NL;
num_lcon = num_lcon + NL;
nnz_As = nnz_As + nnz(this.I_edge_path) + NL;

if this.invoke_method == 3 || this.invoke_method == 4
    if isfield(this.lower_bounds, 'node')
        num_lcon = num_lcon + NC;
        nnz_As = nnz_As + NC*NV;
    end
end
if this.invoke_method == 1 || this.invoke_method == 4
    num_lcon = num_lcon + 2;
    %     nnz_As = nnz_As + this.num_varz + this.num_varv + NP + NL;
    nnz_As = nnz_As + nnz(this.I_node_path)*NV + this.num_varv + NP + NL;
    num_lcon = num_lcon + NC;
    nnz_As = nnz_As + nnz(this.Hrep) + this.num_varv;
end

As = spalloc(num_lcon, num_vars, nnz_As);
%% (1) flow processing requirement
As(1:this.num_lcon_res,1:this.num_vars) = As_res;
row_offset = this.num_lcon_res;
%% (2) VNF instance capacity constraint
As(row_offset+(1:this.num_varv), NP+(1:this.num_varz)) = this.Hdiag;
As(row_offset+(1:this.num_varv), this.num_vars+(1:this.num_varv)) = -eye(this.num_varv);
row_offset = row_offset + this.num_varv;
if this.invoke_method~=1
    %% (3) link reconfiguration cost constraint
    As(row_offset+(1:NP), 1:NP) = eye(NP);
    As(row_offset+(1:NP), num_base_vars+(1:NP)) = -eye(NP);
    row_offset = row_offset + NP;
    As(row_offset+(1:NP), 1:NP) = -eye(NP);
    As(row_offset+(1:NP), num_base_vars+(1:NP)) = -eye(NP);
    row_offset = row_offset + NP;
    %% (4) VNF allocation reconfiguration cost constraint
    As(row_offset+(1:this.num_varz), (NP+1):this.num_vars) = eye(this.num_varz);
    As(row_offset+(1:this.num_varz), (num_base_vars+NP)+(1:this.num_varz)) = -eye(this.num_varz);
    row_offset = row_offset + this.num_varz;
    As(row_offset+(1:this.num_varz), (NP+1):this.num_vars) = -eye(this.num_varz);
    As(row_offset+(1:this.num_varz), (num_base_vars+NP)+(1:this.num_varz)) = -eye(this.num_varz);
    row_offset = row_offset + this.num_varz;
    %% (5) VNF instance capacity reconfiguration cost constraint
    As(row_offset+(1:this.num_varv), (this.num_vars+1):num_base_vars) = eye(this.num_varv);
    As(row_offset+(1:this.num_varv), (num_base_vars+this.num_vars+1):num_base_vars*2) = -eye(this.num_varv);
    row_offset = row_offset + this.num_varv;
    As(row_offset+(1:this.num_varv), (this.num_vars+1):num_base_vars) = -eye(this.num_varv);
    As(row_offset+(1:this.num_varv), (num_base_vars+this.num_vars+1):num_base_vars*2) = -eye(this.num_varv);
    row_offset = row_offset + this.num_varv;
end
%% (6) virtual link capacity constraint
% (optional) resource reservation for individual links
As(row_offset+(1:NL), 1:NP) = this.I_edge_path;
if this.invoke_method == 1 || this.invoke_method == 4
    theta = options.Reserve;
    theta1 = (1+theta)/2;
else
    theta1 = 1;
end
if this.invoke_method == 1
    As(row_offset+(1:NL), num_base_vars+(1:NL)) = -theta1*eye(NL);
else
    As(row_offset+(1:NL), 2*num_base_vars+(1:NL)) = -theta1*eye(NL);
end
row_offset = row_offset + NL;
%% (8) lower bound of node capacity
if this.invoke_method == 3 || this.invoke_method == 4
    if isfield(this.lower_bounds, 'node')
        As(row_offset+(1:NC), this.num_vars+(1:this.num_varv)) = -repmat(eye(NC),1,NV);
        row_offset = row_offset + NC;
    end
end

%% (10) resource utilization lower bound for resource reservation
if this.invoke_method == 1 || this.invoke_method == 4
    %     As(row_offset+1, [NP+(1:this.num_varz), this.num_vars+(1:this.num_varv)]) = ...
    %         [ones(1,this.num_varz), -theta*ones(1,this.num_varv)];
    As(row_offset+1, [NP+(1:this.num_varz), this.num_vars+(1:this.num_varv)]) = ...
        [repmat((this.I_node_path(:))', 1, NV), -theta*ones(1,this.num_varv)];
    As(row_offset+2, [1:NP, (num_vars-NL+1):num_vars]) = ...
        [sum(this.I_edge_path,1), -theta*ones(1, NL)];
    row_offset = row_offset + 2;
    %% (11) virtual node capacity constraint for resource reservation
    % (optional) resource reservation for individual nodes
    As(row_offset+(1:NC), NP+(1:this.num_varz)) = this.Hrep;
    As(row_offset+(1:NC), this.num_vars+(1:this.num_varv)) = -theta1*repmat(eye(NC),1,NV);
end

%% (1)(2) right side
bs = sparse(this.num_lcon_res+this.num_varv,1);
%% (3)(4)(5) right side
if this.invoke_method ~= 1
    bs = [bs;   
        this.topts.old_variables_x;       % which have been pre-processed, so it can be
        -this.topts.old_variables_x;      % compared with the current states.
        this.topts.old_variables_z;
        -this.topts.old_variables_z;
        this.old_variables.v;
        -this.old_variables.v];
end
%% (6) right side
bs = [bs; zeros(NL,1)];
%% (8) right side
if this.invoke_method == 3 || this.invoke_method == 4 
    if isfield(this.lower_bounds, 'node')
        bs = [bs; -this.lower_bounds.node];
    end
end
%% (10)(11) right side
if this.invoke_method == 1 || this.invoke_method == 4
    bs = [bs; 0; 0; zeros(NC,1)];
end

%% (7)(9)
lbs = sparse(num_vars, 1);
if this.invoke_method == 3 || this.invoke_method == 4
    lbs(2*num_base_vars+(1:NL)) = this.lower_bounds.link;
    if isfield(this.lower_bounds, 'VNF')
        lbs(this.num_vars+(1:this.num_varv)) = this.lower_bounds.VNF;
    end
end

if this.invoke_method == 1
    var0 = zeros(num_vars,1);
else
    var0 = [this.topts.old_variables_x;
        this.topts.old_variables_z;
        this.old_variables.v;
        this.topts.old_variables_x;
        this.topts.old_variables_z;
        this.old_variables.v];
end
if this.invoke_method == 3 || this.invoke_method == 4
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
fmincon_opt.Display = 'iter';   %'notify-detailed'; %'iter'; 'notify'
% fmincon_opt.CheckGradients = true;
% fmincon_opt.FiniteDifferenceType = 'central';
% fmincon_opt.Diagnostics = 'on';
% options.Form = 'normal';
if strcmpi(options.Form, 'compact')
    % column/variables: isequal(this.I_active_variables', sum(this.As_res,1)~=0)
    % row/constraints:  isequal(active_rows, find(sum(this.As_res,2)~=0))  
    z_filter = sparse(repmat(...
        reshape(logical(this.I_node_path), numel(this.I_node_path),1), NV, 1));
    this.I_active_variables = [true(NP,1);  z_filter;  true(this.num_varv,1)];
    if this.invoke_method ~= 1
        this.I_active_variables = [this.I_active_variables;...
        true(NP,1); z_filter; true(this.num_varv,1)];    
    end
    this.I_active_variables = [this.I_active_variables; true(NL,1)];
    if this.invoke_method ~= 1
        row_offset = this.num_lcon_res + this.num_varv + 2*NP;
        active_rows = [(1:row_offset)'; row_offset+find(z_filter); ...
            row_offset+this.num_varz+find(z_filter)];
        row_offset = row_offset + 2*this.num_varz;
        active_rows = [active_rows; row_offset+(1:2*this.num_varv)';...
            row_offset+2*this.num_varv+(1:NL)'];
        row_offset = row_offset + 2*this.num_varv + NL;
        if isfield(this.lower_bounds, 'node')
            active_rows = [active_rows; row_offset+(1:NC)'];
            row_offset = row_offset + NC;
        end
        if this.invoke_method == 4
            active_rows = [active_rows; row_offset+(1:2)'; row_offset+2+(1:NC)'];
        end
    else
        active_rows = (1:num_lcon)';
    end
    As = As(active_rows, this.I_active_variables);
    var0 = var0(this.I_active_variables);
    lbs = lbs(this.I_active_variables);
    bs = bs(active_rows);
    options.num_orig_vars = num_vars;   % original number of variables in the problem
    options.bCompact = true;
end
options.num_varx = this.NumberPaths;
options.num_varz = this.num_varz;   % ONLY need the original number of variables
options.num_varv = this.num_varv;
if this.invoke_method == 1
    fmincon_opt.HessianFcn = ...
        @(x,lambda)DynamicSlice.hessInitialSlicing(x, lambda, this, options);
    [xs, fval, exitflag, output] = ...
        fmincon(@(x)DynamicSlice.fcnProfitReserveSlicing(x, this, options), ...
        var0, As, bs, [], [], lbs, [], [], fmincon_opt);
else
    fmincon_opt.HessianFcn = ...
        @(x,lambda)DynamicSlice.hessSlicing(x, lambda, this, options);
    [xs, fval, exitflag, output] = ...
        fmincon(@(x)DynamicSlice.fcnProfitReconfigureSlicing(x, this, options), ...
        var0, As, bs, [], [], lbs, [], [], fmincon_opt);    
end
if strcmpi(options.Form, 'compact')
    x = zeros(num_vars, 1);
    x(this.I_active_variables) = xs;
else
    x = xs;
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
%% ISSUE
% if we use the normal form, since we have not filter those zero variables beforehand, 
% in the objective function we must calculate load counting all variables (instead of 
% using <getNodeLoad> and <getLinkLoad>).

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