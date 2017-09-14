%% Optimize flow rate in the slice
% This optimization procedure remove the unnecessary components from the independent
% variable x, so that the problem scale is cut down.
% if h(n,p) = 0, then z(n,p,f) = 0 for all f. thus the corresponding column in the
% coefficient matrix and the variables can be removed.
function x = priceOptimalFlowRateCompact(this, x0)
if nargin == 2
    this.x0 = x0;
else%if isempty(this.x0)
    % give a feasible start point
    this.x0 = zeros(this.num_vars,1);
    this.x0(1:this.NumberPaths) = 1;
    max_alpha_f = max(this.Parent.VNFTable.ProcessEfficiency(this.VNFList));
    this.x0((this.NumberPaths+1):end) = 1*this.NumberPaths*max_alpha_f;
end
% the number of linear constraint without considering the bound constraints.
nnz_As = this.NumberVNFs*(this.NumberPaths+nnz(this.I_node_path));
As = spalloc(this.num_lcon_res, this.num_vars, nnz_As);
row_index = 1:this.NumberPaths;
col_index = this.NumberPaths+(1:this.NumberVirtualNodes);
for f = 1:this.NumberVNFs
    alpha_f = this.Parent.VNFTable.ProcessEfficiency(this.VNFList(f));
    As(row_index,1:this.NumberPaths) = alpha_f * eyesparse(this.NumberPaths); %#ok<SPRIX>
    for p = 1:this.NumberPaths
        As(row_index(p), col_index) = -this.I_node_path(:,p); %#ok<SPRIX>
        col_index = col_index + this.NumberVirtualNodes;
    end
    row_index = row_index + this.NumberPaths;
end
% remove the all-zero column from As
this.I_active_variable = sum(As,1)~=0;
As(:, ~this.I_active_variable) = [];
bs = sparse(this.num_lcon_res,1);
num_act_vars = nnz(this.I_active_variable);
lbs = sparse(num_act_vars,1);
x0 = this.x0(this.I_active_variable);
%% Set the optimization options
% * *Algorithm* -- since the problem contains linear constraints and bound
% constraints, then |trust-region-reflective| method is not applicable. Hence,
% we choose the |interior point| method. As a result the Hessian matrix should
% be computed separately.
% is directly returned from the objective function as the second derivatives.
% * *HessianFcn* -- we compute Hessian using the objective function.
% Therefore, this option is set as |'objective'|.
% * *SpecifyObjectiveGradient* -- the gradient can be directly computed from
% the objective function, so this option is set to |true|.
% * *SpecifyConstraintGradient*: since this problem does not contain nonlinear
% constraint, this option is set to |false|.
% * *Display information*: use |'iter'| to display iteration information for
% debug.
fmincon_opt = optimoptions(@fmincon);
fmincon_opt.Algorithm = 'interior-point';
fmincon_opt.HessianFcn = @(x,lambda)Slice.hessianfcn(x,lambda,this);
fmincon_opt.SpecifyObjectiveGradient = true;
fmincon_opt.SpecifyConstraintGradient = false;
fmincon_opt.Display = 'iter';
[x, fval, exitflag] = fmincon(@(x)Slice.funUtilityCompact(x,this), ...
    x0, As, bs, [], [], lbs, [], [], fmincon_opt);
disp(fval);
disp(exitflag);
if exitflag == 1
    %% Output solution
    % the output
    this.x0 = zeros(this.num_vars,1);
    this.x0(I_active_variable) = x;
    this.x_path = this.x0(1:this.NumberPaths);
    this.z_npf = this.x0((this.NumberPaths+1):end);
    this.x_path(this.x_path<(10^-4)*max(this.x_path)) = 0;
    this.z_npf(this.z_npf<(10^-4)*max(this.z_npf)) = 0;
    this.z_node_path_vnf = reshape(this.z_npf, ...
        this.NumberVirtualNodes, this.NumberPaths, this.NumberVNFs);
    this.link_load = this.getLinkLoad(this.x0(1:this.NumberPaths));
    this.node_load = this.getNodeLoad(this.x0((this.NumberPaths+1):end));
else
    warning('%d',exitflag);
end
end