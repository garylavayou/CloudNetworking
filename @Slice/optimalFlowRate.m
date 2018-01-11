%% Optimal Flow Rate
% optimize flow rate with link and node capacity constraints.
% Since the capacity constraints is considered in the problem, the solution can utilize
% multi-path to serve flows.
%% Model
% * *Accurate*
% * *Approximate*
% * *FixedCost*: used when slice's resource amount and resource prices are fixed, resource
% cost is fixed. See also <fcnSocialWelfare>, <DynamicSlice>.
%% Options: 
% _Methods_
%	*single-normal*: combine all flows in one slice, and solve it as a global
%           optimization; 
%   *single-function*:
%   *slice-price*: solve each slice's problem independently, with resource price known.
%    resource prices should be set beforehand.
% _CostModel_: 'fixed'.
% _Form_: 'compact'|'normal'.
% _PricingPolicy_: 'quadratic'|'linear'.
% NOTE: check options before computing.

%% Output
% If output argument is provided, _optimalFlowRate_ will calculate the final results.
% Otherwise, the results is stored in temporary variables.
function [profit, cost] = optimalFlowRate( this, new_opts )
if nargin < 2
    new_opts = struct;
end
options = structmerge(new_opts, ...
    getstructfields(this.Parent.options, 'Form'), ...
    getstructfields(this.options, 'PricingPolicy', 'default', 'linear'), ...
    getstructfields(new_opts, 'Method'),...
    'exclude');     

NL = this.NumberVirtualLinks;
NC = this.NumberDataCenters;
NV = this.NumberVNFs;
NP = this.NumberPaths;
% Coefficient for process-rate constraints
nnz_As = nnz(this.As_res) + NV*nnz(this.I_node_path) + nnz(this.I_edge_path);
num_lcon = this.num_lcon_res + NC + NL;
parameters.As = spalloc(num_lcon, this.num_vars, nnz_As);
parameters.As(1:this.num_lcon_res,:) = this.As_res;

%% Add node and link capacity constraint coefficient
row_index = this.num_lcon_res + (1:NC);
col_index = NP + (1:this.num_varz);
parameters.As(row_index, col_index) = this.Hrep;
row_index = (1:NL) + this.num_lcon_res + NC;
parameters.As(row_index, 1:NP) = this.I_edge_path; 

%% Boundary
% * *Resource Constraints*: processing-rate, the right side is all-zero; capacity
% constraints, the right side is the capacity of link and node;
parameters.bs = [sparse(this.num_lcon_res,1); 
    this.VirtualDataCenters.Capacity;
    this.VirtualLinks.Capacity];
%%%
% * *Upper Bound*: Not necessary, to facilitate the algorithm, we give a relaxed
% upper-bound.
% max(this.Parent.getLinkField('Capacity'))*ones(NP,1);...
%     max(this.Parent.getDataCenterField('Capacity'))*ones(this.num_vars-NP,1)
parameters.ub = [];
%%%
% * Remove the capacity constraints with infinity capacity;
idx = find(parameters.bs==inf);
parameters.As(idx,:) = [];         
parameters.bs(idx) = [];

%% Feasible test of start point
% * *Start Point*: in case that the capacity of a virtual link/node is zero, we initialize
% $z_{min}$ and $x_{min}$ as the nonzero minimum value.
parameters.x0 = zeros(this.num_vars,1);
switch options.Method
    case 'single-function'
        max_alpha_f = max(options.Alpha_f);
    otherwise
        max_alpha_f = max(this.Parent.VNFTable{this.VNFList, 'ProcessEfficiency'});
end
% z_min = min(this.VirtualNodes.Capacity(this.VirtualNodes.Capacity>1))/(NP*NV);
% x_min = min(this.VirtualLinks.Capacity(this.VirtualLinks.Capacity>1))/NP;
% if z_min == inf || x_min == inf
%     x0(1:NP) = 1;
%     x0((NP+1):end) = max_alpha_f;
% else
%     if z_min >= max_alpha_f*x_min  % [x,z] is feasible
%         x0(1:NP) = x_min;
%         x0((NP+1):end) = z_min;
%     else
%         x0(1:NP) = z_min/max_alpha_f;
%         x0((NP+1):end) = z_min;
%     end
% end
parameters.x0(1:NP) = 1;
parameters.x0((NP+1):end) = max_alpha_f;
assert(this.checkFeasible(parameters.x0), 'error: infeasible start point.');
parameters.Aeq = [];
parameters.beq = [];

options.fmincon_opt = optimoptions(@fmincon);
options.fmincon_opt.Algorithm = 'interior-point';
options.fmincon_opt.SpecifyObjectiveGradient = true;
options.fmincon_opt.Display = 'notify';
if strcmpi(options.Form, 'compact')
    z_filter = sparse(repmat(...
        reshape(logical(this.I_node_path), numel(this.I_node_path),1), NV, 1));
    this.I_active_variables = [true(NP,1) ;  z_filter];
    parameters.As = parameters.As(:, this.I_active_variables);
    parameters.x0 = parameters.x0(this.I_active_variables);
    parameters.lb = sparse(length(parameters.x0),1);
    if ~isempty(parameters.ub)
        parameters.ub = parameters.ub(this.I_active_variables);
    end
    options.num_orig_vars = this.num_vars; 
    options.fmincon_opt.HessianFcn = ...
        @(x,lambda)Slice.fcnHessianCompact(x, lambda, this, options);
else
    if contains(options.Method, 'price')
        options.fmincon_opt.HessianFcn = ...
            @(x,lambda)Slice.fcnHessian(x, lambda, this, options);
    else
        options.fmincon_opt.HessianFcn = ...
            @(x,lambda)Slice.fcnHessian(x, lambda, this);
    end
    parameters.lb = sparse(this.num_vars,1);
end
[x, fval] = this.optimize(parameters, options);

%% Additional Process
% Under the problem formulation, if node |n| is not used by path |p|($h_{np} = 0$) or path
% |p| does not use VNF |f|, the corresponding $z_{npf}$ should be zero. However, the
% optimization results of the related $z_{npf}$ can be arbitrary value, since the related 
% $z_{npf}$ is not counted in the cost items. As a result, we should manually set
% the value of those related $z_{npf}$ to zero.
%
% On the other hand, too small components should be rounded.
this.temp_vars.x = x(1:NP);
switch options.Method
    case 'single-function'
        %% TODO
        % allocate VNF instance by order.
        this.VNFList = 1:this.Parent.NumberVNFs;
        this.getAs_res;
        z_np = x(NP+1:end);
        z_np = reshape(this.I_node_path(:).*z_np, NC,NP);
        znpf = zeros(NC, NP, this.NumberVNFs);
        for v = 1:this.NumberVNFs
            af = this.Parent.VNFTable.ProcessEfficiency(v);
            idx_path = this.I_path_function(:,v)==1;
            p_slice = this.path_owner(idx_path);
            znpf(:, idx_path, v) = z_np(:, idx_path).*(af./this.alpha_f(p_slice))';
        end
        this.temp_vars.z = znpf(:);
    otherwise
        this.temp_vars.z = x((NP+1):end);
        nz = this.NumberDataCenters*this.NumberPaths;
        z_index = 1:nz;
        for f = 1:this.NumberVNFs
            this.temp_vars.z(z_index) = this.I_node_path(:).*this.temp_vars.z(z_index);
            z_index = z_index + nz;
        end
end
clear znpf;
this.flow_rate = this.getFlowRate(this.temp_vars.x);
profit = -fval;
if nargout >= 1    % final results
    this.Variables.x = this.temp_vars.x;
    this.Variables.z = sparse(this.temp_vars.z);
    this.postProcessing();
    this.setPathBandwidth;
    this.FlowTable.Rate = this.getFlowRate;
    this.VirtualLinks.Load = this.getLinkLoad;
    this.VirtualDataCenters.Load = this.getNodeLoad;
end
if nargout >= 2
    cost = this.getSliceCost(options.PricingPolicy);   % method is overrided by subclasses.
end
end