%% Fast Reconfiguration with VNF Instance Re-scaling
% Comparing with <_fastReconfigure_>, this method also considers VNF instance re-scaling,
% i.e., the capacity of VNF instances could change during reconfiguration.
% See also <fastReconfigure>.
function profit = fastReconfigure2(this, action, options)
global InfoLevel;
NL = this.NumberVirtualLinks;
NN = this.NumberDataCenters;
NP = this.NumberPaths;
NV = this.NumberVNFs;

%%% Formulate input for convex optimization (fmincon).
% The problem has multiple inequalities, and the lowerbounds for variables.
global DEBUG; %#ok<NUSED>
Hd = this.Hdiag;
Hr = this.Hrep;
As_res = this.As_res;        % update As_res
%%
% List of constaints:
%   (1) flow processing requirement: NP*NV (this.num_lcon_res);
%   (2) VNF instance capacity constraint: NV*NN (this.num_varv);
%   (3) Node capacity constraint: NN;
%   (4) Link capacity constraint: NL;
%   (5) link reconfiguration cost constraint: 2*NP;
%   (6) node reconfiguration cost constraint: 2*NN*NP*NV (2*this.num_varz);
%   (5) VNF instance reconfiguration cost constraint: 2*NV*NN;
num_lcon = this.num_lcon_res + this.num_varv + NN + NL + ...
    2*NP + 2*this.num_varz + 2*this.num_varv;
nnz_As = nnz(As_res) + (nnz(Hd)+this.num_varv) + nnz(Hr) + nnz(this.I_edge_path) + ...
    + 4*NP + 4*this.num_varz+ 4*this.num_varv;
num_vars = 2*this.num_vars + 2*this.num_varv;
As = spalloc(num_lcon, num_vars, nnz_As);
As(1:this.num_lcon_res,1:this.num_vars) = As_res;
row_offset = this.num_lcon_res;
As(row_offset+(1:this.num_varv), NP+(1:this.num_varz)) = Hd;
As(row_offset+(1:this.num_varv), this.num_vars+(1:this.num_varv)) = -eye(this.num_varv);
row_offset = row_offset + this.num_varv;
As(row_offset+(1:NN), NP+(1:this.num_varz)) = Hr;
row_offset = row_offset + NN;
As(row_offset+(1:NL), 1:NP) = this.I_edge_path;
row_offset = row_offset + NL;
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
    this.VirtualDataCenters.Capacity; % The field will only change in slice dimensionging.
    this.VirtualLinks.Capacity;
    this.topts.old_variables_x;       % which have been pre-processed, so it can be
    -this.topts.old_variables_x;      % compared with the current states.
    this.topts.old_variables_z;
    -this.topts.old_variables_z;
    this.old_variables.v;      % Equal to last stage's VNF capacity, size not change
    -this.old_variables.v];

var0 = [this.topts.old_variables_x;
    this.topts.old_variables_z;
    this.old_variables.v;
    this.topts.old_variables_x;
    this.topts.old_variables_z;
    this.old_variables.v];
assert(this.checkFeasible(var0), 'error: infeasible start point.');

%% Perform optimization
fmincon_opt = optimoptions(@fmincon);
fmincon_opt.Algorithm = 'interior-point';
fmincon_opt.SpecifyObjectiveGradient = true;
fmincon_opt.Display = InfoLevel.InnerModel.char;
if nargin >= 2 && strcmpi(options.Form, 'compact')
    %% get the compact formulation
    % There are lots of zeros in $z_npf$, which could be determined by $h_np$.
    % If a row is all zero, this row is in-active, and could be removed.
    z_filter = sparse(repmat(...
        reshape(logical(this.I_node_path), numel(this.I_node_path),1),...
        NV,1));
    this.I_active_variables = [true(NP,1);  z_filter;  true(this.num_varv,1);...
        true(NP,1); z_filter; true(this.num_varv,1)];
    row_offset = this.num_lcon_res + this.num_varv + NN + NL + 2*NP;
    active_rows = [(1:row_offset)'; row_offset+find(z_filter); ...
        row_offset+this.num_varz+find(z_filter); ...
        ((num_lcon-2*this.num_varv+1):num_lcon)'];
    As_compact = As(active_rows, this.I_active_variables);
    var0_compact = var0(this.I_active_variables);
    bs = bs(active_rows);
    fcn_opts.num_varx = this.NumberPaths;
    fcn_opts.num_varz = nnz(z_filter);
    fcn_opts.num_varv = this.num_varv;
    this.topts.z_reconfig_cost = this.topts.z_reconfig_cost(z_filter);
    lbs = sparse(length(var0_compact),1);
    fmincon_opt.HessianFcn = ...
        @(x,lambda)DynamicSlice.fcnHessianCompact(x,lambda,this,fcn_opts);
    [x_compact, fval, exitflag] = ...
        fmincon(@(x)DynamicSlice.fcnFastConfigProfitCompact(x,this,fcn_opts), ...
        var0_compact, As_compact, bs, [], [], lbs, [], [], fmincon_opt);
    x = zeros(num_vars, 1);
    x(this.I_active_variables) = x_compact;
else
    lbs = sparse(2*this.num_vars,1);
    fmincon_opt.HessianFcn = @(x,lambda)Slice.fcnHessian(x,lambda,this);
    [x, fval, exitflag] = ...
        fmincon(@(x)DynamicSlice.fcnFastConfigProfit(x,this), ...
        var0, As, bs, [], [], lbs, [], [], fmincon_opt);
end
%% Reconfiguration Cost in Problem Formulation
% comparing the old variables with new variables to decide the reconfiguration
% cost.
% Since the two vector has different number of elements, we should comparing
% it accordingly, and set aside the variables of new arriving/departing flow.

% x is a local solution to the problem when exitflag is positive.
if exitflag == 0
    warning('reaching maximum number of iterations.');
elseif exitflag < 0
    error('abnormal exit with flag %d.',exitflag);
elseif exitflag ~= 1
    warning('(exitflag = %d) local optimal solution found.', exitflag);
end
options.Action = action;    % This might be used when check feasible solution.
options.ConstraintTolerance = fmincon_opt.ConstraintTolerance;
assert(this.checkFeasible(x,options), 'error: infeasible solution.');

%%%
% The post processing like <optimalFlowRate> is not needed, since the
% objective function (reconfiguration cost) will force those variables to be zero.
this.convert(x, 0);
this.flow_rate = this.getFlowRate(this.temp_vars.x);
this.postProcessing(struct('VNFCapacity',true));
this.setPathBandwidth;
this.FlowTable.Rate = this.getFlowRate;
this.VirtualLinks.Load = this.getLinkLoad;
this.VirtualDataCenters.Load = this.getNodeLoad;
profit = -fval - this.getSliceCost(options.PricingPolicy);
end