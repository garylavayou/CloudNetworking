%% priceOptimalFlowRate
% Extending <Slice.priceOptimalFlowRate>, Considering resource reconfiguration cost.
% NOTE: update price before call this function.
function [net_profit, node_load, link_load] = priceOptimalFlowRate(this, x0, new_opts) %#ok<INUSL>
global DEBUG INFO;
if nargin <= 2
    new_opts = struct;
end
% if this.options.ReconfigMethod == ReconfigMethod.Baseline ||...
%         this.options.ReconfigMethod == ReconfigMethod.DimBaseline || ...
%         ~this.options.bReserve && this.b_derive_vnf && isempty(fieldnames(this.net_changes))
%     if nargout <= 1
%         net_profit = priceOptimalFlowRate@Slice(this, x0, new_opts);
%     else
%         [net_profit, node_load, link_load] = priceOptimalFlowRate@Slice(this, x0, new_opts);
%     end
%     return;
% end

theta0 = 0.999;
options = structmerge(new_opts, ...
    getstructfields(this.Parent.options, 'Form', 'default', 'normal'), ...
    getstructfields(this.options, {'PricingPolicy', 'Reserve'}, 'default', {'quadratic',theta0}),...
    'exclude');     

Np = this.NumberPaths;
Nvnf = this.NumberVNFs;
Nl = this.NumberLinks;
Ndc = this.NumberServiceNodes;
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
if this.b_derive_vnf || this.options.ReconfigMethod == ReconfigMethod.DimBaseline
    this.invoke_method = 1;  % initialization
elseif isempty(this.lower_bounds)  % bReserve is true
    this.invoke_method = 2;  % no resource reservation (theta = 0.99 to tolerate error)
elseif options.Reserve >= 1  
    this.invoke_method = 3;  % implicit resource reservation
elseif options.Reserve > 0
    this.invoke_method = 4;  % explicit resource reservation
elseif options.Reserve < 0 && options.Reserve > -1
    this.invoke_method = 5;
    options.Reserve = -options.Reserve;
else 
    error('%s: invalid value ''Reserve=%.2f''for resource reservation.', calledby, options.Reserve);
end
if ~this.b_derive_vnf % this.invoke_method ~= 1
    this.update_reconfig_cost(this.sh_options.action, true);
end
num_base_vars = this.num_vars + this.num_varv;
num_vars = num_base_vars;
num_lcon = this.num_lcon_res + this.num_varv;
nnz_As = nnz(As_res) + (nnz(this.Hdiag)+this.num_varv);
if this.invoke_method ~= 1
    num_vars = num_vars + num_base_vars;
    num_lcon = num_lcon + ...
        2*Np + 2*this.num_varz + 2*this.num_varv;
    nnz_As = nnz_As + 4*Np + 4*this.num_varz+ 4*this.num_varv;
end 
num_vars = num_vars + Nl;
num_lcon = num_lcon + Nl;
nnz_As = nnz_As + nnz(this.I_edge_path) + Nl;

if this.invoke_method == 3 || this.invoke_method == 4
    if isfield(this.lower_bounds, 'node')         
        if isfield(options, 'ResidualCapacity')
            node_idx = options.ResidualCapacity.Node-this.lower_bounds.node>1;
        else
            node_idx = true(Ndc,1);
        end
        nn_bounds = nnz(node_idx);
        num_lcon = num_lcon + nn_bounds;
        nnz_As = nnz_As + nn_bounds*Nvnf;
    end
end
if (this.invoke_method == 1 || this.invoke_method == 4) && this.options.bReserve
    % overall resource utilization constraint is mandated.
    theta = options.Reserve;
    num_lcon = num_lcon + 2;
    %     nnz_As = nnz_As + this.num_varz + this.num_varv + NP + NL;
    nnz_As = nnz_As + nnz(this.I_dc_path)*Nvnf + this.num_varv + Np + Nl;
end
if this.invoke_method == 1 || this.invoke_method == 4
    switch this.options.bReserve 
        case 1 % only resource reservation to achieve overall utilization ratio
            theta1 = theta0;
        case 2 % resource reservation for individual resource (per node)
            theta1 = min((1+theta)/2, theta0);
            num_lcon = num_lcon + Ndc;
            nnz_As = nnz_As + nnz(this.Hrep) + this.num_varv;
        case 3 % resource reservation for individual resource (per VNF instance)
            % add coefficient to constraint (2), no new constraints.
            theta1 = min((1+theta)/2, theta0);
            theta0 = theta1;
        case 0 % no resource reservation
            theta1 = theta0;
    end
else
    % If |theta0| is set to 1, no resource is reserved for individual resource;
    % Otherwise, we set |theta0| very close to 1, to handle the precision error.
    theta1 = theta0;
end

% if this.invoke_method == 4
%     idx_up_node = find(this.upper_bounds.node>0);
%     num_lcon = num_lcon + numel(idx_up_node);
%     nnz_As = nnz_As + numel(idx_up_node)*NV;
%     idx_up_link = find(this.upper_bounds.link>0);
%     num_lcon = num_lcon + numel(idx_up_link);
%     nnz_As = nnz_As + nnz(this.I_edge_path(idx_up_link,:));
% end

As = spalloc(num_lcon, num_vars, nnz_As);
bs = sparse(num_lcon,1);
lbs = sparse(num_vars, 1);
%% (1) flow processing requirement
As(1:this.num_lcon_res,1:this.num_vars) = As_res;
% bs(1:this.num_lcon_res+(1:NL)) = zeros(this.num_lcon_res,1);
row_offset = this.num_lcon_res;
%% (2) VNF instance capacity constraint
As(row_offset+(1:this.num_varv), Np+(1:this.num_varz)) = this.Hdiag;
As(row_offset+(1:this.num_varv), this.num_vars+(1:this.num_varv)) = -theta0*eye(this.num_varv);
% bs(row_offset+(1:this.num_varv)) = zeros(this.num_varv,1);
row_offset = row_offset + this.num_varv;
if this.invoke_method~=1
    %% (3) link reconfiguration cost constraint
    As(row_offset+(1:Np), 1:Np) = eye(Np);
    As(row_offset+(1:Np), num_base_vars+(1:Np)) = -eye(Np);
    bs(row_offset+(1:Np)) = this.topts.old_variables_x; % which have been pre-processed, so it can be compared with the current states.
    row_offset = row_offset + Np;
    As(row_offset+(1:Np), 1:Np) = -eye(Np);
    As(row_offset+(1:Np), num_base_vars+(1:Np)) = -eye(Np);
    bs(row_offset+(1:Np)) = -this.topts.old_variables_x;
    row_offset = row_offset + Np;
    %% (4) VNF allocation reconfiguration cost constraint
    As(row_offset+(1:this.num_varz), (Np+1):this.num_vars) = eye(this.num_varz);
    As(row_offset+(1:this.num_varz), (num_base_vars+Np)+(1:this.num_varz)) = -eye(this.num_varz);
    bs(row_offset+(1:this.num_varz)) = this.topts.old_variables_z;
    row_offset = row_offset + this.num_varz;
    As(row_offset+(1:this.num_varz), (Np+1):this.num_vars) = -eye(this.num_varz);
    As(row_offset+(1:this.num_varz), (num_base_vars+Np)+(1:this.num_varz)) = -eye(this.num_varz);
    bs(row_offset+(1:this.num_varz)) = -this.topts.old_variables_z;
    row_offset = row_offset + this.num_varz;
    %% (5) VNF instance capacity reconfiguration cost constraint
    As(row_offset+(1:this.num_varv), (this.num_vars+1):num_base_vars) = eye(this.num_varv);
    As(row_offset+(1:this.num_varv), (num_base_vars+this.num_vars+1):num_base_vars*2) = -eye(this.num_varv);
    bs(row_offset+(1:this.num_varv)) = this.old_variables.v;
    row_offset = row_offset + this.num_varv;
    As(row_offset+(1:this.num_varv), (this.num_vars+1):num_base_vars) = -eye(this.num_varv);
    As(row_offset+(1:this.num_varv), (num_base_vars+this.num_vars+1):num_base_vars*2) = -eye(this.num_varv);
    bs(row_offset+(1:this.num_varv)) = -this.old_variables.v;
    row_offset = row_offset + this.num_varv;
end
%% (6) virtual link capacity constraint
% (optional) resource reservation for individual links
As(row_offset+(1:Nl), 1:Np) = this.I_edge_path;
if this.invoke_method == 1
    As(row_offset+(1:Nl), num_base_vars+(1:Nl)) = -theta1*eye(Nl);
else
    As(row_offset+(1:Nl), 2*num_base_vars+(1:Nl)) = -theta1*eye(Nl);
end
% bs(row_offset+(1:NL)) = zeros(NL,1);
row_offset = row_offset + Nl;
%% (8) lower bound of node capacity
if this.invoke_method == 3 || this.invoke_method == 4
    if isfield(this.lower_bounds, 'node')
        reduced_eye = eye(Ndc);
        reduced_eye = reduced_eye(node_idx,:);
        As(row_offset+find(node_idx), this.num_vars+(1:this.num_varv)) = -repmat(reduced_eye,1,Nvnf);
        bs(row_offset+find(node_idx)) = -this.lower_bounds.node(node_idx);
        row_offset = row_offset + nn_bounds;
    end
end

%% (10) resource utilization lower bound for resource reservation
if this.invoke_method == 1 || this.invoke_method == 4 
    if this.options.bReserve  % overall resource utilization constraint
        % As(row_offset+1, [NP+(1:this.num_varz), this.num_vars+(1:this.num_varv)]) = ...
        % 	[ones(1,this.num_varz), -theta*ones(1,this.num_varv)];
        As(row_offset+1, [Np+(1:this.num_varz), this.num_vars+(1:this.num_varv)]) = ...
            [repmat((this.I_dc_path(:))', 1, Nvnf), -theta*ones(1,this.num_varv)];
        As(row_offset+2, [1:Np, (num_vars-Nl+1):num_vars]) = ...
            [sum(this.I_edge_path,1), -theta*ones(1, Nl)];
        % bs(row_offset+(1:2)) = [0; 0];
        row_offset = row_offset + 2;
    end
    if this.options.bReserve == 2 % per node 
        %% (11) virtual node capacity constraint for resource reservation
        % (optional) resource reservation for individual nodes
        As(row_offset+(1:Ndc), Np+(1:this.num_varz)) = this.Hrep;
        As(row_offset+(1:Ndc), this.num_vars+(1:this.num_varv)) = -theta1*repmat(eye(Ndc),1,Nvnf);
        % bs(row_offset+(1:NC)) = zeros(NC,1);
        row_offset = row_offset + Ndc;
    elseif this.options.bReserve == 3
        
    end
    % (per VNF instance is appied to constriant (2))
end
%% (12)
% row_offset = row_offset + NC;
% if this.invoke_method == 4
%     if ~isempty(idx_up_node)
%         d = repmat(eye(NC),1,NV);
%         As(row_offset+(1:numel(idx_up_node)), this.num_vars+(1:this.num_varv)) = d(idx_up_node,:);
%         row_offset = row_offset + numel(idx_up_node);
%         bs = [bs; this.upper_bounds.node(idx_up_node)];
%     end
%     if ~isempty(idx_up_link)
%         As(row_offset+(1:numel(idx_up_link)), this.num_vars+(1:NP)) = this.I_edge_path(idx_up_link,:);
%         bs = [bs; this.upper_bounds.link(idx_up_link)];
%     end
% end
%% (7)(9)
if (this.invoke_method == 3 || this.invoke_method == 4)
    if isfield(options, 'ResidualCapacity')
        link_idx = options.ResidualCapacity.Link-this.lower_bounds.link>1;
    else
        link_idx = true(Nl,1);
    end
    lbs(2*num_base_vars+find(link_idx)) = this.lower_bounds.link(link_idx);
    if isfield(this.lower_bounds, 'VNF')
        vnf_idx = repmat(node_idx, Nvnf, 1);
        lbs(this.num_vars+find(vnf_idx)) = this.lower_bounds.VNF(vnf_idx);
    end
end

if this.invoke_method == 1 
    var0 = zeros(num_vars,1);
else
    var0 = [this.topts.old_variables_x/4;
        this.topts.old_variables_z/2;
        this.old_variables.v;
        zeros(num_base_vars,1); 
        %         this.topts.old_variables_x;
        %         this.topts.old_variables_z;
        %         this.old_variables.v;
        this.old_state.link_capacity];
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
% options.Form = 'normal';
if strcmpi(options.Form, 'compact')
    % column/variables: isequal(this.I_active_variables', sum(this.As_res,1)~=0)
    % row/constraints:  isequal(active_rows, find(sum(this.As_res,2)~=0))  
    z_filter = sparse(repmat(logical(this.I_dc_path(:)), Nvnf, 1));
    this.I_active_variables = [true(Np,1);  z_filter;  true(this.num_varv,1)];
    if this.invoke_method ~= 1
        %         tx_filter = this.topts.x_reconfig_cost~=0;
        %         tz_filter = z_filter&(this.topts.z_reconfig_cost~=0);
        tx_filter = true(numel(this.topts.x_reconfig_cost),1);
        tz_filter = z_filter;
        this.I_active_variables = [this.I_active_variables;...
            tx_filter; tz_filter; true(this.num_varv,1)];
    end
    this.I_active_variables = [this.I_active_variables; true(Nl,1)];
    if this.invoke_method ~= 1
        row_offset = this.num_lcon_res + this.num_varv;
        active_rows = [true(row_offset,1); tx_filter; tx_filter; tz_filter; tz_filter;
            true(2*this.num_varv,1); true(Nl,1)];
        if (this.invoke_method == 4 || this.invoke_method == 3) && ...
                isfield(this.lower_bounds, 'node')
            active_rows = [active_rows; true(nn_bounds,1)];
        end
        if this.invoke_method == 4
            active_rows = [active_rows; true(2,1)]; 
            if this.options.bReserve == 2
                active_rows = [active_rows; true(Ndc,1)];                
            end
        end
    else
        active_rows = true(num_lcon,1);
    end
    % NOTE that the |reconfig_cost| cosefficient is not filtered like the
    % <fastReconfigure>. The objective function <fcnProfitReconfigureSlicing> will use the
    % orginal cost coefficient.
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

%% perform resource reservation
if this.invoke_method == 5
    theta = options.Reserve;
    node_load = this.getNodeLoad(this.temp_vars.z);
    link_load = this.getLinkLoad(this.temp_vars.x);
    num_vars = Nl+Ndc*Nvnf+Ndc*Nvnf;
    var0r = [this.temp_vars.c;
        this.old_variables.v;
        zeros(Ndc*Nvnf,1)];
    lbr = [max(this.lower_bounds.link, link_load);      % options.bDisableLowerBound
        this.Hdiag*this.temp_vars.z;
        zeros(Ndc*Nvnf,1)];
    %%
    num_lcons = Ndc + 2 + 2*Ndc*Nvnf;
    Ar = sparse(num_lcons, num_vars);
    br = sparse(num_lcons, 1);
    Ar(1:Ndc,Nl+(1:Ndc*Nvnf)) = -repmat(eye(Ndc),1,Nvnf);
    br(1:Ndc) = -max(this.lower_bounds.node, node_load);    % options.bDisableLowerBound
    row_offset = row_offset + Ndc;
    Ar(row_offset+1,1:Nl) = -ones(1,Nl);
    br(row_offset+1) = -sum(link_load)/theta;
    row_offset = row_offset + 1;
    Ar(row_offset+1,Nl+(1:Ndc*Nvnf)) = -ones(1,Ndc*Nvnf);
    br(row_offset+1) = -sum(node_load)/theta;
    row_offset = row_offset + 1;
    Ar(row_offset+(1:Ndc*Nvnf), Nl+(1:2*Ndc*Nvnf)) = [eye(Ndc*Nvnf), -eye(Ndc*Nvnf)];
    br(row_offset+(1:Ndc*Nvnf)) = this.old_variables.v;
    row_offset = row_offset + Ndc*Nvnf;
    Ar(row_offset+(1:Ndc*Nvnf), Nl+(1:2*Ndc*Nvnf)) = [-eye(Ndc*Nvnf), -eye(Ndc*Nvnf)];
    br(row_offset+(1:Ndc*Nvnf)) = -this.old_variables.v;
    %%
    % we might add a resource utilization ratio for individual resources.
    fr_opt = optimoptions(@fmincon);
    fr_opt.Algorithm = 'interior-point';
    fr_opt.SpecifyObjectiveGradient = true;
    fr_opt.Display = 'iter';   %'notify-detailed'; %'iter'; 'notify'
    %     fr_opt.CheckGradients = true;
    %     fr_opt.FiniteDifferenceType = 'central';
    fr_opt.HessianFcn = @(x,lambda)hessReconfigCost(x, lambda, this, options);
    [x, ~, exitflag, output] = ...
        fmincon(@(x)fcnReconfigCost(x, this, options), ...
        var0r, Ar, br, [], [], lbr, [], [], fr_opt);
    this.interpretExitflag(exitflag, output.message);
    this.temp_vars.c = x(1:Nl);
    this.temp_vars.v = x(Nl+(1:Ndc*Nvnf));
    this.temp_vars.tv = x(Nl+Ndc*Nvnf+(1:Ndc*Nvnf));
end

% Passing out the load information is not necessary, since the network only reuqire the
% information of occupied capacity.
if nargout >= 2 
    node_load = zeros(this.Parent.NumberServiceNodes,1);
    data_center_id = this.getDCPI;
    node_load(data_center_id) = this.getNodeCapacity(false);
end
if nargout >= 3 
    link_load = zeros(this.Parent.NumberLinks,1);
    link_load(this.Links.PhysicalLink) = this.getLinkCapacity(false);
end
this.flow_rate = this.getFlowRate(this.temp_vars.x);
net_profit = this.getProfit(options);
%%%
% FOR DEBUG
% this.setPathBandwidth(this.temp_vars.x);
end

%% TEST
%{
tx_index = (1:NP) + num_base_vars;
tx = x(tx_index);
node_load = this.getNodeLoad(this.temp_vars.z);
node_capacity = this.getNodeCapacity(false);
disp([node_load, node_capacity])
disp(node_load./node_capacity)
link_load = this.getLinkLoad(this.temp_vars.x);
link_capacity = this.getLinkCapacity(false);
disp('link load | link capacity |  old link capacity');
disp([link_load, link_capacity, this.old_state.link_capacity, ...
    this.lower_bounds.link, this.upper_bounds.link]);
sum(node_load)./sum(node_capacity)
sum(link_load)./sum(link_capacity)
%}

function [f,g] = fcnReconfigCost(vars, this, options)
ND = this.NumberServiceNodes;
NL = this.NumberLinks;
NV = this.NumberVNFs;
link_cap = vars(1:NL);
vnf_cap = vars(NL+(1:ND*NV));
vnf_cap_diff = vars(NL+ND*NV+(1:ND*NV));
node_load = sum(reshape(vnf_cap, ND, NV),2);
switch options.PricingPolicy
    case {'quadratic-price', 'quadratic'}
        [link_payment,link_price_grad] = this.fcnLinkPricing(this.prices.Link, link_cap);
        [node_payment,node_price_grad] = this.fcnNodePricing(this.prices.Node, node_load);
        f = link_payment + node_payment;
    otherwise
        error('%s: invalid pricing policy', calledby);
end
f = f + dot(vnf_cap_diff, this.topts.vnf_reconfig_cost);

g = [link_price_grad; repmat(node_price_grad, NV, 1); this.topts.vnf_reconfig_cost];

end

function h = hessReconfigCost(vars, lambda, this, options) %#ok<INUSL>
ND = this.NumberServiceNodes;
NL = this.NumberLinks;
NV = this.NumberVNFs;
link_cap = vars(1:NL);
vnf_cap = vars(NL+(1:ND*NV));
node_load = sum(reshape(vnf_cap, ND, NV),2);
num_vars = length(vars);
h = spalloc(num_vars,num_vars, num_vars);
switch options.PricingPolicy
    case {'quadratic-price', 'quadratic'}
        [~,~,lph] = this.fcnLinkPricing(this.prices.Link, link_cap);
        [~,~,nph] = this.fcnNodePricing(this.prices.Node, node_load);
        h(1:NL,1:NL) = diag(lph);
        h(NL+(1:ND*NV),NL+(1:ND*NV)) = block_diag(diag(nph), NV);
    otherwise
        error('%s: invalid pricing policy', calledby);
end
end