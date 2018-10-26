%% Reconfiguration Cost in Problem Formulation
% comparing the old variables with new variables to decide the reconfiguration
% cost.
% Since the two vector has different number of elements, we should comparing
% it accordingly, and set aside the variables of new arriving/departing flow.
function [profit, cost] = fastReconfigure(this, action, new_opts)
global computime ITER_LIMIT event_num NUMBER_ITERS;
if nargin <= 2
    new_opts = struct;
end
theta0 = 0.999;
options = structmerge(new_opts, ...
    getstructfields(this.Parent.options, 'Form'), ...
    getstructfields(this.options, 'PricingPolicy', 'default', {'quadratic'}), ...
    getstructfields(this.options, 'ReconfigMethod'), ...
    'exclude');
if ~isfield(options, 'bEnforceReserve')
    options.bEnforceReserve = false;
end
if isfield(this.options, 'penalty') && ~isempty(this.options.penalty)
    options.bDistributed = true;
elseif ~isfield(options, 'bDistributed') || isempty(options.bDistributed)
    options.bDistributed = false;
end

Nf = this.NumberFlows;
if Nf == 0
    [profit, cost] = this.handle_zero_flow(options);
    return;
end
Nl = this.NumberLinks;
Ndc = this.NumberServiceNodes;
Nvnf = this.NumberVNFs;
Np = this.NumberPaths;
this.max_flow_rate = this.FlowTable.Rate;
this.max_flow_rate(this.max_flow_rate==0) = inf;
%%% Formulate input for convex optimization (fmincon).
% The problem has multiple inequalities, and the lowerbounds for variables.
%%
% List of constaints:num_lcon_res
%   (1) flow processing requirement: Np*Nvnf (this.NumberLinearConstraints);
%   (2) VNF instance capacity constraint: NV*Ndc (this.num_varv);
%   (3) Link capacity constraint: NL;
%   (4) link reconfiguration cost constraint: 2*NP;
%   (5) node reconfiguration cost constraint: 2*Ndc*NP*NV (2*this.num_varz);
As_res = this.As_res;        % update As_res
nnz_As = nnz(As_res) + nnz(this.Hdiag) + nnz(this.I_edge_path) + ...
    + 4*Np + 4*this.num_varz; % 2 I_x, 2 I_tx, 2 I_z, 2 I_tz
num_lcon = this.NumberLinearConstraints + Ndc*Nvnf + Nl + 2*Np + 2*this.num_varz;
if this.options.bReserve
    if ~options.bEnforceReserve
        if strcmpi(action, 'add')
            nnz_As = nnz_As + nnz(this.I_flow_path(1:(Nf-1)));
        else
            nnz_As = nnz_As + nnz(this.I_flow_path(1:end));
        end
        num_lcon = num_lcon + Nf;
    end
    if options.bEnforceReserve
        % Fast reconfigure after slice dimensioning.
        theta = this.options.Reserve;
        num_lcon = num_lcon + 2;
        nnz_As = nnz_As + nnz(this.I_dc_path)*Nvnf + Np;
        if this.options.bReserve == 2
            theta1 = min((1+theta)/2, theta0);
            num_lcon = num_lcon + Ndc;
            nnz_As = nnz_As + nnz(this.Hrep);
        elseif this.options.bReserve == 3
            theta1 = min((1+theta)/2, theta0);
            theta0 = theta1;
        else
            theta1 = theta0;
        end
        % elseif
        % Independent fast reconfigure
        % if add
        %% Overall resource utilization constraint and Individual resource utilization;
        %         theta = (this.utilizationRatio()+1)/2;
        %         theta1 = min(theta0, (1+theta)/2);
        % if remove
        %% Overall resource utilization constraint and Individual resource utilization;
        % Issue: individual resource utilization might violate the reosurce granularity.
        %         if this.options.ReconfigMethod == ReconfigMethod.FastconfigReserve
        %             this.options.Reserve = this.utilizationRatio();
        %         end
        %         theta = this.options.Reserve;
        %         theta1 = min(theta0, (1+theta)/2);
    end
end

t1 = tic;
num_vars = 2*this.NumberVariables;     % including axuilary variables.
As = spalloc(num_lcon, num_vars, nnz_As);
As(1:this.NumberLinearConstraints,1:this.NumberVariables) = As_res;
row_offset = this.NumberLinearConstraints;
As(row_offset+(1:Ndc*Nvnf), (1:this.num_varz)+Np) = this.Hdiag;
row_offset = row_offset + Ndc*Nvnf;
As(row_offset+(1:Nl), 1:Np) = this.I_edge_path;
row_offset = row_offset + Nl;
As(row_offset+(1:Np), 1:Np) = speye(Np);
As(row_offset+(1:Np), this.NumberVariables+(1:Np)) = -speye(Np);
row_offset = row_offset + Np;
As(row_offset+(1:Np), 1:Np) = -speye(Np);
As(row_offset+(1:Np), this.NumberVariables+(1:Np)) = -speye(Np);
row_offset = row_offset + Np;
As(row_offset+(1:this.num_varz), (Np+1):this.NumberVariables) = speye(this.num_varz);
As(row_offset+(1:this.num_varz), (this.NumberVariables+Np+1):end) = -speye(this.num_varz);
row_offset = row_offset + this.num_varz;
As(row_offset+(1:this.num_varz), (Np+1):this.NumberVariables) = -speye(this.num_varz);
As(row_offset+(1:this.num_varz), (this.NumberVariables+Np+1):end) = -speye(this.num_varz);
row_offset = row_offset + this.num_varz;
idx_empty = [];
if this.options.bReserve
    if ~options.bEnforceReserve
        if strcmpi(action, 'add')
            %% limit the existing flows data rate
            % We assume the last flow in the flow table is the new flow
            [~, limit_idx] = sort(this.FlowTable.Rate, 'descend');
            num_limit1 = floor(0.95*Nf);
            limit_idx1 = limit_idx(1:num_limit1);
            %         mid_rate = median(this.FlowTable.Rate);
            %         limit_idx = find(this.FlowTable.Rate>=mid_rate);
            As(row_offset+(1:num_limit1), 1:Np) = this.I_flow_path(limit_idx1,:);
            %% Limit new flow/low rate flow's data rate
            row_offset = row_offset + num_limit1;
            limit_idx2 = limit_idx((num_limit1+1):(Nf-1));
            As(row_offset+(1:numel(limit_idx2)), 1:Np) = this.I_flow_path(limit_idx2,:);
            row_offset = row_offset + length(limit_idx2) + 1;
            % new arriving flow has no limit. (The corresponding row in As is 0)
            idx_empty = [idx_empty; row_offset];
        else
            %% limit the existing flows data rate
            [~, limit_idx] = sort(this.FlowTable.Rate, 'descend');
            num_limit1 = floor(0.95*this.NumberFlows);
            limit_idx1 = limit_idx(1:num_limit1);
            As(row_offset+(1:num_limit1), 1:Np) = this.I_flow_path(limit_idx1,:);
            row_offset = row_offset + num_limit1;
            limit_idx2 = limit_idx((num_limit1+1):end);
            As(row_offset+(1:numel(limit_idx2)), 1:Np) = this.I_flow_path(limit_idx2,:);
            row_offset = row_offset + length(limit_idx2);
        end
    end
    %% limit the resource utilization
    if options.bEnforceReserve
        As(row_offset+1, Np+(1:this.num_varz)) = repmat((this.I_dc_path(:))', 1, Nvnf);
        As(row_offset+2, 1:Np) = sum(this.I_edge_path,1);
        row_offset = row_offset + 2;
        if this.options.bReserve == 2
            As(row_offset+(1:Ndc), Np+(1:this.num_varz)) = this.Hrep;
            row_offset = row_offset + Ndc; %#ok<NASGU>
        elseif this.options.bReserve == 3
            
        end
    end
end

cv = this.VNFCapacity(:);              % VNFCapacity not change under 'fastconfig'
cl = this.Links.Capacity;
bs0 = [sparse(this.NumberLinearConstraints,1);
    theta0*cv;
    theta0*cl;
    this.topts.old_variables_x;       % which have been pre-processed, so it can be
    -this.topts.old_variables_x;      % compared with the current states.
    this.topts.old_variables_z;
    -this.topts.old_variables_z];
if this.options.bReserve
    if ~options.bEnforceReserve
        % The settings for 'add' and 'remove' might be different
        if strcmpi(action, 'add')
            %% Limit Flow Rate
            bs = [bs0; this.max_flow_rate(limit_idx1); ...
                max(this.max_flow_rate(limit_idx2),1.5*mean(this.FlowTable.Rate)*ones(length(limit_idx2),1));...
                0];
        else
            bs = [bs0; this.max_flow_rate(limit_idx1); ...
                max(this.max_flow_rate(limit_idx2),1.5*mean(this.FlowTable.Rate)*ones(length(limit_idx2),1))];
        end
    else
        bs = bs0;
    end
    if options.bEnforceReserve
        % The resource reservation constraint is consistent with those in <priceOptimalFlowRate>.
        % If sole FSR need resource reservation on individual resources, the constraints might need
        % to consider the "old link/VNF load", which might be higher than the "theta" items.
        bs = [bs; theta*sum(cv); theta*sum(cl)];
        Nvn = length(cv);
        if this.options.bReserve == 2  % individual link/node
            bs(this.NumberLinearConstraints+(1:Nvn)) = theta0*cv;     % VNF instance capacity
            bs(this.NumberLinearConstraints+Nvn+(1:Nl)) = theta1*cl;  % link capacity
            bs = [bs; theta1*sum(reshape(cv, Ndc, Nvnf),2)];
        elseif this.options.bReserve == 3 % % individual link/VNF instance
            % update constraints, no new ones.
            bs(this.NumberLinearConstraints+(1:(Nvn+Nl))) = theta1*[cv;cl];
        end
    end
else
    bs = bs0;
end
var0 = [this.topts.old_variables_x;
    this.topts.old_variables_z;
    this.topts.old_variables_x;
    this.topts.old_variables_z
    %     sparse(this.NumberVariables,1)
    ];
if ~this.checkFeasible(var0)
    warning('infeasible initial point.');
end
t2 = toc(t1);
fprintf('%s: initilizing arguements ... Elapsed time is %f seconds.\n', calledby(0), t2);

%% Perform optimization
fmincon_opt = optimoptions(@fmincon);
fmincon_opt.Algorithm = 'interior-point';
fmincon_opt.SpecifyObjectiveGradient = true;
fmincon_opt.Display = 'notify';     % iter
% fmincon_opt.CheckGradients = true;
% fmincon_opt.FiniteDifferenceType = 'central';
t1 = tic;
if options.bDistributed
    fmincon_opt.OptimalityTolerance = 1e-5;     % In-exact solving the sub-problems
    [x, fval, k] = distribute_optimization([],[],[],0);
else
    if strcmpi(options.Form, 'compact')
        %% get the compact formulation
        % There are lots of zeros in z_npf, which could be determined by h_np.
        z_filter = sparse(repmat(logical(this.I_dc_path(:)), Nvnf,1));
        %     tx_filter = this.topts.x_reconfig_cost~=0;
        %     tz_filter = z_filter&(this.topts.z_reconfig_cost~=0);
        tx_filter = true(numel(this.topts.x_reconfig_cost),1);
        tz_filter = z_filter;
        this.I_active_variables = [true(Np,1);  z_filter; tx_filter; tz_filter];
        row_offset = this.NumberLinearConstraints + Ndc*Nvnf + Nl;
        active_rows = [true(row_offset,1); tx_filter; tx_filter; tz_filter; tz_filter];
        if this.options.bReserve
            if ~options.bEnforceReserve
                if strcmpi(action, 'add')
                    active_rows = [active_rows; true(Nf-1,1); false];
                else
                    active_rows = [active_rows; true(Nf,1)];
                end
            end
            if options.bEnforceReserve
                active_rows = [active_rows; true(2,1)];
                if this.options.bReserve == 2
                    active_rows = [active_rows; true(Ndc,1)];
                end
            end
        end
        As = As(active_rows, this.I_active_variables);
        var0 = var0(this.I_active_variables);
        bs = bs(active_rows);
        fcn_opts.num_varz = nnz(z_filter);
        this.topts.z_reconfig_cost = this.topts.z_reconfig_cost(tz_filter);
        this.topts.x_reconfig_cost = this.topts.x_reconfig_cost(tx_filter);
        fcn_opts.bCompact = true;
    else
        fcn_opts.num_varz = this.num_varz;
        As(idx_empty,:) = [];
        bs(idx_empty,:) = [];
    end
    lbs = sparse(length(var0),1);
    fcn_opts.num_varx = this.NumberPaths;
    fmincon_opt.HessianFcn = ...
        @(x,lambda)DynamicSlice.hessReconfigure(x, lambda, this, fcn_opts);
    
    % ISSUE: Second evaluation without resource reservation constraint when flow arriving
    % b_duo_opt = true;
    % while true
    [xs, fval, exitflag, output] = ...
        fmincon(@(x)DynamicSlice.fcnFastConfigProfit(x, this ,fcn_opts), ...
        var0, As, bs, [], [], lbs, [], [], fmincon_opt);
    if strcmpi(options.Form, 'compact')
        x = zeros(num_vars, 1);
        x(this.I_active_variables) = xs;
    else
        x = xs;
    end
    
    % x is a local solution to the problem when exitflag is positive.
    this.interpretExitflag(exitflag, output.message);
end
t2 = toc;
if ~isempty(computime)
    computime(event_num-1) = t2;
end
if ~isempty(NUMBER_ITERS)
    NUMBER_ITERS(event_num-1) = k;
end
options.Action = action;    % This might be used when check feasible solution.
options.ConstraintTolerance = fmincon_opt.ConstraintTolerance;
assert(this.checkFeasible(x,options), 'error: infeasible solution.');
%%%
% The post processing like <optimalFlowRate> is not needed, since the
% objective function will force those variables to be zero.
this.temp_vars.x = x(1:Np);
this.temp_vars.z = x((Np+1):this.NumberVariables);
this.flow_rate = this.getFlowRate(this.temp_vars.x);
this.postProcessing();
this.FlowTable.Rate = this.getFlowRate;
this.max_flow_rate(this.max_flow_rate==inf) = 0;
this.max_flow_rate = max(this.FlowTable.Rate, this.max_flow_rate);
%     if strcmpi(action, 'add') && this.options.bReserve && b_duo_opt
%         fidx = find(this.FlowTable.Rate<=0.1*median(this.FlowTable.Rate), 1);
%         if ~isempty(fidx)
%             bs = [bs0; this.FlowTable.Rate(limit_idx1)];
%             bs = bs(active_rows(1:end-2-numel(limit_idx2)));
%             As = As(1:end-2-numel(limit_idx2),:);
%             b_duo_opt = false;
%             continue;
%         end
%     end
%     break;
% end

this.setPathBandwidth;
this.Links.Load = this.getLinkLoad;
this.ServiceNodes.Load = this.getNodeLoad;
if nargout >= 1
    cost = this.getCost(options.PricingPolicy, 'const');
    rc_linear = this.get_reconfig_cost('linear', true);
    profit = -fval - cost + rc_linear;
end


%% Dual ADMM to solve large scale fast slice reconfiguration problem
% NOTE: currently, we only implement the basic scheme. The schemes with resoruce reservation needs
% addition constraints to be added and the objective function of subproblems should be added with
% multiplier iterms with the new constraints.
% NOTE2: shared variables must also be copyed to each worker, since the worker might be
% distributed, and thus the shared information might not be available if not copyed.
% NOTE3: It is assumed that flow indices and path indices are continuous for each
% sub-problem.
    function [x, fval, k] = distribute_optimization(num_process, r, tols, opt_order)
        %% parameters
        if nargin<1 || isempty(num_process)
            num_process = floor(Nf);
        else
            num_process = min(num_process, floor(Nf));
        end
        % num_process <= NF
        if nargin <= 1 || isempty(r)
            if isfield(this.options, 'penalty')
                r = this.options.penalty;
            else
                r = 5;
            end
        end
        if nargin <= 2 || isempty(tols)
            eps_abs = 10^-6;
            eps_rel = 10^-3;
        else
            eps_abs = tols.abs;
            eps_rel = tols.rel;
        end
        if nargin <= 3 || isempty(opt_order)
            opt_order = 1;
        end
        if license('test', 'Distrib_Computing_Toolbox')
            p = gcp('nocreate');
        else
            p = [];
        end
        if isempty(p) || contains(fmincon_opt.Display, 'iter')
            M = 0;
        else
            M = p.NumWorkers;
        end
        %% Initialization
        % Number of flows in sub-problem: depending on how many thread is created. The flows are
        % evenly partitioned to each sub-problem.
        num_flows = floor(Nf/num_process)*ones(num_process,1);
        num_rem_flows = mod(Nf,num_process);
        num_flows(end-(num_rem_flows-1):end) = num_flows(1) + 1;
        assert(Nf==sum(num_flows));
        flow_idx_offset = cumsum([0; num_flows(1:end-1)]);
        % Number of dual variables (for basic scheme), corresponding to the VNF instance and link
        % capacity constraints.
        % NOTE: If implementing the scheme with resource reservation, there will be more dual
        % variables (corresponding to the extra constraints).
        num_dual = Ndc*Nvnf+Nl;
        var_indices = cell(num_process, 1);
        distr_num_vars = zeros(num_process,1);   % Number of primal variables for each sub-problem
        distr_xk = cell(num_process, 1);
        param_fields = {'A', 'b', 'cl', 'cr', 'lb', 'num_varx', ...
            'num_varz', 'pidx', 'fidx', 'x_reconfig_cost', 'z_reconfig_cost', ...
            'weight', 'I_flow_path', 'path_owner'};
        if strcmpi(options.Form, 'compact')
            param_fields = [param_fields, ...
                {'I_active_variables', 'active_rows'}];
        end
        distr_params = structarray(num_process, param_fields);
        I_flow_path = this.I_flow_path;
        I_dc_path = this.I_dc_path;
        topts = this.topts;
        path_owner = this.path_owner;
        num_lcon_res = this.NumberLinearConstraints;
        num_varz = this.num_varz;
        weight = this.weight;
        if strcmpi(options.Form, 'compact')
            bCompact = true;
        else
            bCompact = false;
        end
        t1 = tic;
        parfor (i = 1:num_process,M)
            %for i = 1:num_process   % DEBUG
            flow_indices = flow_idx_offset(i) + (1:num_flows(i));
            path_indices = (find(sum(I_flow_path(flow_indices,:),1)))'; %#ok<PFBNS>
            num_paths = length(path_indices);
            dup_path_idx = path_indices+(0:Np:(Nvnf-1)*Np);
            constr_indices = dup_path_idx(:);
            constr_idx_offset = num_lcon_res + Ndc*Nvnf + Nl;
            constr_indices = [constr_indices; constr_idx_offset+path_indices; ...
                constr_idx_offset+Np+path_indices];
            constr_idx_offset = constr_idx_offset + 2*Np;
            dup_z_idx = (path_indices(1)-1)*Ndc + (1:num_paths*Ndc);
            dup_z_idx = dup_z_idx(:) + (0:Ndc*Np:(Nvnf-1)*Ndc*Np);
            dup_z_idx = dup_z_idx(:);
            constr_indices = [constr_indices; constr_idx_offset+dup_z_idx; ...
                constr_idx_offset+num_varz+dup_z_idx];
            constr_idx_offset = constr_idx_offset + 2*num_varz; %#ok<NASGU>
            var_indices{i} = [path_indices; Np+dup_z_idx; Np+num_varz+path_indices;...
                Np*2+num_varz+dup_z_idx];
            distr_num_vars(i) = length(var_indices{i});
            distr_xk{i} = var0(var_indices{i});      %#ok<PFBNS> % Primal variables
            distr_params(i).lb = sparse(distr_num_vars(i),1);
            distr_params(i).A = As(constr_indices, var_indices{i}); %#ok<PFBNS>
            distr_params(i).b = bs(constr_indices); %#ok<PFBNS>
            distr_params(i).cl = As(num_lcon_res+(1:num_dual), var_indices{i});
            distr_params(i).cr = bs(num_lcon_res+(1:num_dual))*(num_paths/Np);
            distr_params(i).num_varx = num_paths;
            distr_params(i).num_varz = length(dup_z_idx);
            distr_params(i).pidx = path_indices;
            distr_params(i).fidx = flow_indices;
            distr_params(i).x_reconfig_cost = topts.x_reconfig_cost(path_indices); %#ok<PFBNS>
            distr_params(i).z_reconfig_cost = topts.z_reconfig_cost(dup_z_idx);
            distr_params(i).weight = weight;
            distr_params(i).I_flow_path = I_flow_path(flow_indices, path_indices);
            distr_params(i).path_owner = path_owner(path_indices); %#ok<PFBNS>
            if bCompact
                distr_node_path = I_dc_path(:,path_indices); %#ok<PFBNS>
                distr_z_filter = repmat(logical(distr_node_path(:)), Nvnf,1);
                distr_tx_filter = true(num_paths,1);
                distr_tz_filter = distr_z_filter;
                distr_params(i).I_active_variables = [true(num_paths,1);  ...
                    distr_z_filter; distr_tx_filter; distr_tz_filter];
                distr_params(i).active_rows = [true(num_paths*Nvnf,1); ...
                    distr_tx_filter; distr_tx_filter; distr_tz_filter; distr_tz_filter];
                %                 if this.options.bReserve
                %                     if ~options.bEnforceReserve
                %                     end
                %                     if options.bEnforceReserve
                %                     end
                %                 end
                distr_params(i).A = distr_params(i).A(distr_params(i).active_rows,...
                    distr_params(i).I_active_variables);
                distr_params(i).b = distr_params(i).b(distr_params(i).active_rows);
                distr_xk{i} = distr_xk{i}(distr_params(i).I_active_variables);
                distr_params(i).cl = distr_params(i).cl(:,distr_params(i).I_active_variables);
                distr_params(i).lb = distr_params(i).lb(distr_params(i).I_active_variables);
                distr_params(i).num_varz = nnz(distr_z_filter);
                distr_params(i).x_reconfig_cost = distr_params(i).x_reconfig_cost(distr_tx_filter);
                distr_params(i).z_reconfig_cost = distr_params(i).z_reconfig_cost(distr_tz_filter);
            end
        end
        t2 = toc(t1);
        fprintf('Dual-ADMM: distributing data ... Elapsed time is %f seconds.\n', t2);
        if ~isempty(this.init_gamma_k)
            gamma_k = repmat(this.init_gamma_k, 1, num_process);
            q_k = repmat(this.init_q_k, 1, num_process);
        else
            gamma_k = zeros(num_dual, num_process); % auxiliary variables to gamma.
            q_k = ones(num_dual, num_process); % Dual variables for the dual ADMM formulation.
        end
        fval_gamma_k = zeros(num_process,1);
        fval_lambda_k = 0;
        exitflag_k = zeros(num_process,1);
        output_k = cell(num_process,1);
        fcnObjective = @fcnAugmentedPrimalBasic;
        hessDual = @hessDualBasic;
        if opt_order~=1
            lambda_k = (sum(gamma_k,2) - sum(q_k,2)/r)/num_process;
        end
        
        k = 1;
        t1 = tic;
        while true
            fval_gamma_kminus = fval_gamma_k;
            fval_lambda_kminus = fval_lambda_k;
            if opt_order == 1
                gamma_kminus = gamma_k;
                lambda_k = (sum(gamma_k,2) - sum(q_k,2)/r)/num_process;
            end
            parfor (j = 1:num_process,M)
                %for j = 1:num_process
                fminopt = fmincon_opt;
                xk = distr_xk{j};
                fminopt.HessianFcn = @(x,lbd) ...
                    hessDual(x,lbd,lambda_k,q_k(:,j),r,distr_params(j));
                % dual_k is not the true dual function value, since the constant part of
                % the objective function has been omitted.
                [xk, fval_gamma_k(j), exitflag_k(j), output_k{j}] = ...
                    fmincon(@(x) fcnObjective(x,lambda_k,q_k(:,j),r,distr_params(j)), ...
                    xk, distr_params(j).A, distr_params(j).b, [], [], distr_params(j).lb, [], [], fminopt); %#ok<PFOUS>
                gamma_k(:,j) = max(0,...
                    lambda_k+1/r*(q_k(:,j)+distr_params(j).cl*xk-distr_params(j).cr));
                distr_xk{j} = xk;
            end
            if opt_order ~= 1
                lambda_kminus = lambda_k;
                lambda_k = (sum(gamma_k,2) - sum(q_k,2)/r)/num_process;
            end
            % fval_lambda_k = sum(q_k'*lambda_k);  % the L2-norm part is counted when computing lambda
						fval_lambda_k = sum(sum(q_k.*gamma_k)); 
            q_k = q_k + r*(lambda_k-gamma_k);
            re = gamma_k - lambda_k;
            re_norm = norm(re(:));
            if opt_order == 1
                se = r*sum(gamma_k-gamma_kminus,2);
                se_norm = norm(se);
                tol_dual = eps_abs*sqrt(num_dual) + eps_rel*norm(sum(q_k,2));       % eps_rel*|A'*y|_2
            else
                se = r*(lambda_k-lambda_kminus);
                se_norm = sqrt(num_process)*norm(se);
                tol_dual = eps_abs*sqrt(num_dual*num_process) + eps_rel*norm(q_k);       % eps_rel*|A'*y|_2
            end
            %% Stop Condition
            % Number of primal variables: in the dual-ADMM formulation, the primal
            %   variables include: (a) the dual variables of the original problems, (b)
            %   the auxiliary variables used to construct the equality system $Î»-Î³i=0$.
            %   The size of primal variables thus is |num_dual+num_dual*num_process|;
            % Number of dual variables: the dual variables are the multipliers
            %   corresponding to the equality constraints $Î»-Î³i=0$. Thus the size is
            %   |num_dual*num_process|;
            tol_primal = eps_abs*sqrt(num_dual*num_process) + ...
                eps_rel*max(sqrt(num_process)*norm(lambda_k), norm(gamma_k(:)));
            fval_change = (sum(fval_gamma_k)+fval_lambda_k)-(sum(fval_gamma_kminus)+fval_lambda_kminus);
            optimal_gap = abs(fval_change/(sum(fval_gamma_k)+fval_lambda_k));
            
            num_iters = 0;
            num_funccount = 0;
            for i = 1:num_process
                num_iters = num_iters + output_k{i}.iterations;
                num_funccount = num_funccount + output_k{i}.funcCount;
            end
            num_iters = round(num_iters/num_process);
            num_funccount = round(num_funccount/num_process);
            if mod(k,10) == 1
                fprintf('                                              Primal-                  Dual-                   Sub-       Sub- \n');
                fprintf('Iteration Step-length  Dual-change/Percent    optimality   Tolerance   optimality   Tolerance  Iterations Evaluations\n');
                cprintf('*text', ...
                        'â€”â?”â?”â?”â?”â?”â?”â?”â?? â€”â?”â?”â?”â?”â?”â?”â?”â?”â?”â?? â€”â?”â?”â?”â?”â?”â?”â?”â?”â?”â?”â?”â?”â?”â?”â?”â?”â?”â?”â?”â?”â?? â€”â?”â?”â?”â?”â?”â?”â?”â?”â?”â?”â?? â€”â?”â?”â?”â?”â?”â?”â?”â?”â?”â?? â€”â?”â?”â?”â?”â?”â?”â?”â?”â?”â?”â?? â€”â?”â?”â?”â?”â?”â?”â?”â?”â?”â?? â€”â?”â?”â?”â?”â?”â?”â?”â?”â?? â€”â?”â?”â?”â?”â?”â?”â?”â?”â?”â?”\n');
            end
            fprintf('%8d  %10.2f  %10G %10.4G   %10.4G  %10.4G   %10.4G  %10.4G  %9d   %9d \n',...
                k, r, fval_change, optimal_gap, re_norm, tol_primal, se_norm, tol_dual, num_iters, num_funccount);

            if mod(k,10) == 0
                fprintf('\n');
            end
            b_stop = false;
            if re_norm < tol_primal %&& abs(fval_change)< 100 %&& se_norm < tol_dual
                b_stop = true;
            end
            if b_stop || k>=ITER_LIMIT
                break;
            end
            k = k + 1;
        end
        t2 = toc(t1);
        fprintf('Dual-ADMM: elapsed time: %d\n', t2);
        
        x = zeros(num_vars, 1);
        for i = 1:num_process
            xi_temp = zeros(distr_num_vars(i), 1);
            if bCompact
                xi_temp(distr_params(i).I_active_variables) = distr_xk{i};
                x(var_indices{i}) = xi_temp;
            else
                x(var_indices{i}) = distr_xk{i};
            end
        end
        
        options.num_varx = Np;
        options.num_varz = this.num_varz;
        fval = this.fcnFastConfigProfit(x, this, options);
        this.init_gamma_k = mean(gamma_k, 2);
        this.init_q_k = mean(q_k, 2);
    end


end

%% The objective function used with the basic scheme without resource reservation.
% See also <fcnFastConfigProfit>.
function [profit, grad] = fcnAugmentedPrimalBasic(vars, lambda, q, r, params)
num_vars = length(vars);
num_basic_vars = params.num_varx + params.num_varz;
num_var_tx = length(params.x_reconfig_cost);
num_var_tz = length(params.z_reconfig_cost);
var_x = vars(1:params.num_varx);
var_tx_index = num_basic_vars+(1:num_var_tx);
var_tx = vars(var_tx_index);
var_tz_index = num_basic_vars+num_var_tx+(1:num_var_tz);
var_tz = vars(var_tz_index);
flow_rate = params.I_flow_path * var_x;      % see also <getFlowRate>

% objective value
profit = -params.weight*sum(fcnUtility(flow_rate));
profit = profit + dot(var_tx, params.x_reconfig_cost) + ...
    dot(var_tz, params.z_reconfig_cost);
profit = profit + r/2*sum((max(0, lambda+1/r*(q+params.cl*vars-params.cr))).^2);

% Gradient
if nargout >= 2
    grad = zeros(num_vars, 1);
    for p = 1:params.num_varx
        i = params.path_owner(p) - params.fidx(1) + 1;
        grad(p) = -params.weight/(1+params.I_flow_path(i,:)*var_x);
    end
    grad(var_tx_index) = params.x_reconfig_cost;
    grad(var_tz_index) = params.z_reconfig_cost;
    
    z_k = lambda+1/r*(q+params.cl*vars-params.cr);
    idx = z_k<=0;
    A = params.cl;
    A(idx,:) = 0;
    grad = grad + (A')*max(0, z_k);
end

end
%% The objective function used with the scheme with resource reservation
% (rate limitation).
function fcnAugmentedPrimalReserve() %#ok<DEFNU>
end
%% The objective function used with the scheme with resource reservation
% rate and resource utilization limitation
function fcnAugmentedPrimalEnforceReserve() %#ok<DEFNU>
end

function h = hessDualBasic(vars, lbd, lambda, q, r, params) %#ok<INUSL>
var_x = vars(1:params.num_varx);
num_vars = length(vars);
h = zeros(num_vars, num_vars);
for p = 1:params.num_varx
    i = params.path_owner(p) - params.fidx(1) + 1;
    h(p,1:params.num_varx) = params.weight *...
        params.I_flow_path(i,:)/(1+(params.I_flow_path(i,:)*var_x))^2;
end

z_k = lambda+1/r*(q+params.cl*vars-params.cr);
idx = z_k<=0;
A = params.cl;
A(idx,:) = 0;
h = h + (A')*A;
end

%% TEST
%{
fprintf('                                              Primal-                  Dual-                   Sub-       Sub- \n');
fprintf('Iteration Step-length       Dual-change       optimality   Tolerance   optimality   Tolerance  Iterations Evaluations\n');
cprintf('*text', ...
        'â€”â?”â?”â?”â?”â?”â?”â?”â?? â€”â?”â?”â?”â?”â?”â?”â?”â?”â?”â?? â€”â?”â?”â?”â?”â?”â?”â?”â?”â?”â?”â?”â?”â?”â?”â?”â?”â?”â?”â?”â?”â?? â€”â?”â?”â?”â?”â?”â?”â?”â?”â?”â?”â?? â€”â?”â?”â?”â?”â?”â?”â?”â?”â?”â?? â€”â?”â?”â?”â?”â?”â?”â?”â?”â?”â?”â?? â€”â?”â?”â?”â?”â?”â?”â?”â?”â?”â?? â€”â?”â?”â?”â?”â?”â?”â?”â?”â?? â€”â?”â?”â?”â?”â?”â?”â?”â?”â?”â?”\n');
fprintf('%8d  %10.2f  %10G %10.8G   %10.8G  %10.8G   %10.8G  %10.8G  %9d   %9d \n',...
    100, 15.67,  12.3456, 0.012, 0.0011, 0.00012, 1.0234, 0.12345, 20, 45);
%}
