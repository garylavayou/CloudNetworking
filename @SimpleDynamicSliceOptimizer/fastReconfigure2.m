%% Fast Reconfiguration with VNF Instance Re-scaling
% Comparing with <_fastReconfigure_>, this method also considers VNF instance re-scaling,
% i.e., the capacity of VNF instances could change during reconfiguration.
% See also <fastReconfigure>.
function [profit,cost] = fastReconfigure2(this, action, options)
global DEBUG; %#ok<NUSED>
if nargin <= 2
	options = Dictionary;
else
	options = Dictionary(options);
end
slice = this.hs;
options = setdefault(options, slice.options, {'PricingPolicy'});   

if this.NumberFlows == 0
    [profit, cost] = this.handle_zero_flow(options);
    return;
end

Nl = slice.NumberLinks;
Nsn = slice.NumberServiceNodes;
Np = slice.NumberPaths;
Nvnf = slice.NumberVNFs;

%%% Formulate input for convex optimization (fmincon).
% The problem has multiple inequalities, and the lowerbounds for variables.
As_res = this.As_res;        % update As_res
%% Save the VNF capacity to the previous state, berfore optimization
% 'FastReconfig2' reconfigure VNF instance. After the slice is created, |VNFCapacity| is recorded.
% After each optimization, the |VNFCapacity| is updated.
this.old_variables.v = this.Variables.v;  
this.update_reconfig_costvinfo();              % update reconfigure cost with scaler.
%%
% List of constraints:
%   (1) flow processing requirement: NP*NV (num_linconstr);
%   (2) VNF instance capacity constraint: NV*NN (this.num_varv); VNF load
%       (y_nf) is no more than VNF instance capacity (v_nf);
%   (3) Node capacity constraint: NN; VNF instance capacity is no more than
%       node capacity. Since there is reconfiguration cost, we cannot use
%       VNF load to express the constraint (Hrep*z). Instead we directly
%       use VNF instance capacity as variables to express it
%       (sum(v_nf)<=v_n). See also <Slice.getHrep>. 
%   (4) Link capacity constraint: NL;
%   (5) link reconfiguration cost constraint: 2*NP;
%   (6) node reconfiguration cost constraint: 2*NN*NP*NV (2*this.num_varz);
%   (5) VNF instance reconfiguration cost constraint: 2*NV*NN;
this.num_vars = [Np; Nvnf*Nsn*Np; Nsn*Nvnf];		% reconfigure x_p and z_npf.
num_vars = sum(this.num_vars)*2;
num_varz = this.num_vars(2);
num_varv = this.num_vars(3);
num_varxz = sum(this.num_vars(1:2));
num_varxzv = sum(this.num_vars(1:3));
num_lcon_res = size(this.As_res,1);
num_lcon = num_lcon_res + num_varv + Nsn + Nl + 2*Np + 2*num_varz + 2*num_varv;
nnz_As = nnz(As_res) + (nnz(this.Hdiag)+num_varv) + nnz(this.Hrep) + nnz(this.I_edge_path) + ...
	+ 4*Np + 4*num_varz+ 4*num_varv;
As = spalloc(num_lcon, num_vars, nnz_As);
As(1:num_lcon_res,1:num_varxz) = As_res;
row_offset = num_lcon_res;
As(row_offset+(1:num_varv), Np+(1:num_varz)) = this.Hdiag;
As(row_offset+(1:num_varv), num_varxz+(1:num_varv)) = -eye(num_varv);
row_offset = row_offset + num_varv;
As(row_offset+(1:Nsn), num_varxz+(1:this.num_varv)) = repmat(eye(Nsn),1,Nvnf);
row_offset = row_offset + Nsn;
As(row_offset+(1:Nl), 1:Np) = this.I_edge_path;
row_offset = row_offset + Nl;
As(row_offset+(1:Np), 1:Np) = eye(Np);
As(row_offset+(1:Np), num_varxzv+(1:Np)) = -eye(Np);
row_offset = row_offset + Np;
As(row_offset+(1:Np), 1:Np) = -eye(Np);
As(row_offset+(1:Np), num_varxzv+(1:Np)) = -eye(Np);
row_offset = row_offset + Np;
As(row_offset+(1:num_varz), (Np+1):num_varxz) = eye(num_varz);
As(row_offset+(1:num_varz), (num_varxzv+Np)+(1:num_varz)) = -eye(num_varz);
row_offset = row_offset + num_varz;
As(row_offset+(1:num_varz), (Np+1):num_varxz) = -eye(num_varz);
As(row_offset+(1:num_varz), (num_varxzv+Np)+(1:num_varz)) = -eye(num_varz);
row_offset = row_offset + num_varz;
As(row_offset+(1:num_varv), (num_varxz+1):num_varxzv) = eye(num_varv);
As(row_offset+(1:num_varv), (num_varxzv+num_varxz+1):num_vars) = -eye(num_varv);
row_offset = row_offset + num_varv;
As(row_offset+(1:num_varv), (num_varxz+1):num_varxzv) = -eye(this.num_varv);
As(row_offset+(1:this.num_varv), (num_varxzv+num_varxz+1):num_vars) = -eye(num_varv);

bs = [sparse(num_lcon_res+num_varv,1);
    this.ServiceNodes.Capacity; % The field will only change in slice dimensionging.
    this.Links.Capacity;
    this.topts.old_variables.x;       % which have been pre-processed, so it can be
    -this.topts.old_variables.x;      % compared with the current states.
    this.topts.old_variables.z;
    -this.topts.old_variables.z;
    this.old_variables.v;      % Equal to last stage's VNF capacity, size not change
    -this.old_variables.v];

var0 = [this.topts.old_variables.x;
    this.topts.old_variables.z;
    this.old_variables.v;
    %     this.topts.old_variables.x;
    %     this.topts.old_variables.z;
    %     this.old_variables.v
    zeros(num_varxzv,1)
    ];
assert(this.checkFeasible(var0), 'error: infeasible start point.');

%% Perform optimization
minopts = optimoptions(@fmincon);
minopts.Algorithm = 'interior-point';
minopts.SpecifyObjectiveGradient = true;
minopts.Display = 'notify';
this.problem.num_vars = this.num_vars;
if strcmpi(this.options.Form, 'compact')
    %% get the compact formulation
    % There are lots of zeros in $z_npf$, which could be determined by $h_np$.
    % If a row is all zero, this row is in-active, and could be removed.
    z_filter = sparse(repmat(...
        reshape(logical(this.I_node_path), numel(this.I_node_path),1),...
        Nvnf,1));
    this.I_active_variables = [true(Np,1);  z_filter;  true(num_varv,1);...
        true(Np,1); z_filter; true(num_varv,1)];
    row_offset = num_lcon_res + num_varv + Nsn + Nl + 2*Np;
    active_rows = [(1:row_offset)'; row_offset+find(z_filter); ...
        row_offset+num_varz+find(z_filter); ...
        ((num_lcon-2*num_varv+1):num_lcon)'];
    As = As(active_rows, this.I_active_variables);
    var0 = var0(this.I_active_variables);
    bs = bs(active_rows);
    this.problem.num_vars(2) = nnz(z_filter);
    this.topts.z_reconfig_cost = this.topts.z_reconfig_cost(z_filter);
    options.bCompact = true;
else
		options.bCompact = false; 
end
lbs = sparse(length(var0),1);
minopts.HessianFcn = ...
    @(x,lambda)SimpleDynamicSliceOptimizer.hessReconfigure(x, lambda, this, options);
[xs, fval, exitflag, output] = ...
    fmincon(@(x)SimpleDynamicSliceOptimizer.fcnFastConfigProfit(x, this, options), ...
    var0, As, bs, [], [], lbs, [], [], minopts);
if options.bCompact
    x = zeros(num_vars, 1);
    x(this.I_active_variables) = xs;
else
    x = xs;
end

%% Reconfiguration Cost in Problem Formulation
% comparing the old variables with new variables to decide the reconfiguration
% cost.
% Since the two vector has different number of elements, we should comparing
% it accordingly, and set aside the variables of new arriving/departing flow.

% x is a local solution to the problem when exitflag is positive.
this.interpretExitflag(exitflag, output.message);
options.Action = action;    % This might be used when check feasible solution.
options.ConstraintTolerance = minopts.ConstraintTolerance;
assert(this.checkFeasible(x,options), 'error: infeasible solution.');

%%%
% The post processing like <optimalFlowRate> is not needed, since the
% objective function (reconfiguration cost) will force those variables to be zero.
this.convert(x, 0);
this.temp_vars.r = this.getFlowRate(this.temp_vars.x);
this.postProcessing();
slice.FlowTable.Rate = this.getFlowRate;
slice.Links.Load = this.getLinkLoad;
slice.ServiceNodes.Load = this.getNodeLoad;
if nargout >= 1
    cost = slice.getCost('const');
    rc_linear = this.get_reconfig_cost('linear', true);
    profit = -fval - cost + rc_linear;
end
end