function profit = fastReconfigure(this, action, options)
NL = this.NumberVirtualLinks;
NN = this.NumberDataCenters;
NV = this.NumberVNFs;
NP = this.NumberPaths;

%%% Formulate input for convex optimization (fmincon).
% The problem has multiple inequalities, and the lowerbounds for variables.
global DEBUG InfoLevel; %#ok<NUSED>
H = this.Hdiag;
As_res = this.As_res;        % update As_res
nnz_As = nnz(As_res) + nnz(H) + nnz(this.I_edge_path) + ...
    + 4*NP + 4*this.num_varz; % 2 I_x, 2 I_tx, 2 I_z, 2 I_tz
num_lcon = this.num_lcon_res + NN*NV + NL + 2*NP + 2*this.num_varz;
num_vars = 2*this.num_vars;     % including axuilary variables.
As = spalloc(num_lcon, num_vars, nnz_As);
As(1:this.num_lcon_res,1:this.num_vars) = As_res;
row_offset = this.num_lcon_res;
As(row_offset+(1:NN*NV), (1:this.num_varz)+NP) = H;
row_offset = row_offset + NN*NV;
As(row_offset+(1:NL), 1:NP) = this.I_edge_path;
row_offset = row_offset + NL;
As(row_offset+(1:NP), 1:NP) = eye(NP);
As(row_offset+(1:NP), this.num_vars+(1:NP)) = -eye(NP);
row_offset = row_offset + NP;
As(row_offset+(1:NP), 1:NP) = -eye(NP);
As(row_offset+(1:NP), this.num_vars+(1:NP)) = -eye(NP);
row_offset = row_offset + NP;
As(row_offset+(1:this.num_varz), (NP+1):this.num_vars) = eye(this.num_varz);
As(row_offset+(1:this.num_varz), (this.num_vars+NP+1):end) = -eye(this.num_varz);
row_offset = row_offset + this.num_varz;
As(row_offset+(1:this.num_varz), (NP+1):this.num_vars) = -eye(this.num_varz);
As(row_offset+(1:this.num_varz), (this.num_vars+NP+1):end) = -eye(this.num_varz);

bs = [sparse(this.num_lcon_res,1);
    this.VNFCapacity(:);
    this.VirtualLinks.Capacity;
    this.topts.old_variables_x;       % which have been pre-processed, so it can be
    -this.topts.old_variables_x;      % compared with the current states.
    this.topts.old_variables_z;
    -this.topts.old_variables_z];

var0 = [this.topts.old_variables_x;
    this.topts.old_variables_z;
    this.topts.old_variables_x;
    this.topts.old_variables_z];
assert(this.checkFeasible(var0), 'error: infeasible solution.');

%% Perform optimization
fmincon_opt = optimoptions(@fmincon);
fmincon_opt.Algorithm = 'interior-point';
fmincon_opt.SpecifyObjectiveGradient = true;
fmincon_opt.Display = InfoLevel.InnerModel.char;
if nargin >= 2 && strcmpi(options.Form, 'compact')
    %% get the compact formulation
    % There are lots of zeros in z_npf, which could be determined by h_np.
    z_filter = sparse(repmat(...
        reshape(logical(this.I_node_path), numel(this.I_node_path),1),...
        this.NumberVNFs,1));
    this.I_active_variables = [true(this.NumberPaths,1)  ;  z_filter;  ...
        true(this.NumberPaths,1); z_filter];
    row_offset = this.num_lcon_res + NN*NV + NL + 2*NP;
    active_rows = [(1:row_offset)'; row_offset+find(z_filter); ...
        row_offset+this.num_varz+find(z_filter)];
    As_compact = As(active_rows, this.I_active_variables);
    var0_compact = var0(this.I_active_variables);
    bs = bs(active_rows);
    fcn_opts.num_varx = this.NumberPaths;
    fcn_opts.num_varz = nnz(z_filter);
    this.topts.z_reconfig_cost = this.topts.z_reconfig_cost(z_filter);
    lbs = sparse(2*(fcn_opts.num_varx+fcn_opts.num_varz),1);
    fmincon_opt.HessianFcn = ...
        @(x,lambda)DynamicSlice.fcnHessianCompact(x,lambda,this, fcn_opts);
    [x_compact, fval, exitflag] = ...
        fmincon(@(x)DynamicSlice.fcnFastConfigProfitCompact(x,this,fcn_opts), ...
        var0_compact, As_compact, bs, [], [], lbs, [], [], fmincon_opt);
    x = zeros(num_vars, 1);
    x(this.I_active_variables) = x_compact;
else
    lbs = sparse(num_vars,1);
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
    if InfoLevel.UserModelDebug >= DisplayLevel.Notify
        warning('reaching maximum number of iterations.');
    end
elseif exitflag < 0
    error('abnormal exit with flag %d.',exitflag);
elseif exitflag ~= 1
    if InfoLevel.UserModelDebug >= DisplayLevel.Notify
        warning('(exitflag = %d) local optimal solution found.', exitflag);
    end
end
options.Action = action;    % This might be used when check feasible solution.
options.ConstraintTolerance = fmincon_opt.ConstraintTolerance;
assert(this.checkFeasible(x,options), 'error: infeasible solution.');
%%%
% The post processing like <optimalFlowRate> is not needed, since the
% objective function will force those variables to be zero.
this.temp_vars.x = x(1:NP);
this.temp_vars.z = x((NP+1):this.num_vars);
this.flow_rate = this.getFlowRate(this.temp_vars.x);
this.Variables.x = this.temp_vars.x;
this.Variables.z = this.temp_vars.z;
this.postProcessing();
% if ~b
%     if InfoLevel.UserModelDebug >= DisplayLevel.Iteration
%         warning('FastReconfigure: the rounding of variables %s\n%s', ...
%             'with small quantity will make the solution infeasible.',...
%             'successfully recoverd.');
%     end
% end
this.setPathBandwidth;
this.FlowTable.Rate = this.getFlowRate;
this.VirtualLinks.Load = this.getLinkLoad;
this.VirtualDataCenters.Load = this.getNodeLoad;
profit = -fval - this.getSliceCost(options.PricingPolicy);
end
