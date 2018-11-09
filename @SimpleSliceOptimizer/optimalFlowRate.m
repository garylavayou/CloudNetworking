%% Optimal Flow Rate
% optimize flow rate with link and node capacity constraints.
% Since the capacity constraints is considered in the problem, the solution can utilize
% multi-path to serve flows.
%% Model
% * *Accurate*
% * *Approximate*
%% Options:
% _Methods_
%	*single-normal*: combine all flows in one slice, and solve it as a global
%           optimization;
% *single-function*:
% *slice-price*: solve each slice's problem independently, with resource price known.
%    resource prices should be set beforehand.
% _CostModel_: 'fixed'.
% _Form_: 'compact'|'normal'.
% _PricingPolicy_: 'quadratic'|'linear'.
% NOTE: check options before computing.
%% Output
% If output argument is provided, _optimalFlowRate_ will calculate the final results.
% Otherwise, the results is stored in temporary variables.
function [profit, cost] = optimalFlowRate( this, new_opts )
slice = this.hs;
defaultopts = structmerge(...
	getstructfields(this.options, {'Form'}, 'error'), ...
	getstructfields(slice.options, {'PricingPolicy', 'SlicingMethod'}, 'error'));
if nargin < 2
	options = defaultopts;
else
	options = structmerge(defaultopts, new_opts);
end

Ne = slice.NumberLinks;
Nsn = slice.NumberServiceNodes;
Nvf = slice.NumberVNFs;
Np = slice.NumberPaths;
% Coefficient for process-rate constraints
nnz_As = nnz(this.As_res) + Nvf*nnz(this.I_dc_path) + nnz(this.I_edge_path);
num_lcon = this.NumberLinearConstraints + Nsn + Ne;
this.problem.As = spalloc(num_lcon, this.NumberVariables, nnz_As);
this.problem.As(1:this.NumberLinearConstraints,:) = this.As_res;

%% Add node and link capacity constraint coefficient
row_index = this.NumberLinearConstraints + (1:Nsn);
col_index = Np + (1:this.num_varz);
this.problem.As(row_index, col_index) = this.Hrep;
row_index = (1:Ne) + this.NumberLinearConstraints + Nsn;
this.problem.As(row_index, 1:Np) = this.I_edge_path;

%% Boundary
% * *Resource Constraints*: processing-rate, the right side is all-zero; capacity
% constraints, the right side is the capacity of link and node;
this.problem.bs = [sparse(this.NumberLinearConstraints,1);
	slice.ServiceNodes.Capacity;
	slice.Links.Capacity];
%%%
% * *Upper Bound*: Not necessary, to facilitate the algorithm, we give a relaxed
% upper-bound.
% max(this.Parent.readLink('Capacity'))*ones(NP,1);...
%     max(this.Parent.readDataCenter('Capacity'))*ones(this.NumberVariables-NP,1)
this.problem.ub = [];
%%%
% * Remove the capacity constraints with infinity capacity;
idx = find(this.problem.bs==inf);
this.problem.As(idx,:) = [];
this.problem.bs(idx) = [];

%% Feasible test of start point
% * *Start Point*: in case that the capacity of a virtual link/node is zero, we initialize
% $z_{min}$ and $x_{min}$ as the nonzero minimum value.
this.problem.x0 = zeros(this.NumberVariables,1);
if options.SlicingMethod == SlicingMethod.SingleFunction
	max_alpha_f = max(options.Alpha_f);
else
	max_alpha_f = max(slice.Parent.VNFTable{slice.VNFList, 'ProcessEfficiency'});
end
% z_min = min(this.Nodes.Capacity(this.Nodes.Capacity>1))/(NP*NV);
% x_min = min(this.Links.Capacity(this.Links.Capacity>1))/NP;
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
this.problem.x0(1:Np) = 1;
this.problem.x0((Np+1):end) = max_alpha_f;
assert(this.checkFeasible(this.problem.x0), 'error: infeasible start point.');
this.problem.Aeq = [];
this.problem.beq = [];

options.fmincon_opt = optimoptions(@fmincon);
options.fmincon_opt.Algorithm = 'interior-point';
options.fmincon_opt.SpecifyObjectiveGradient = true;
options.fmincon_opt.Display = 'iter';
% options.fmincon_opt.CheckGradients = true;
% options.fmincon_opt.FiniteDifferenceType = 'central';
if strcmpi(options.Form, 'compact')
	z_filter = sparse(repmat(...
		reshape(logical(this.I_dc_path), numel(this.I_dc_path),1), Nvf, 1));
	this.I_active_variables = [true(Np,1) ;  z_filter];
	this.problem.As = this.problem.As(:, this.I_active_variables);
	this.problem.x0 = this.problem.x0(this.I_active_variables);
	this.problem.lb = sparse(length(this.problem.x0),1);
	if ~isempty(this.problem.ub)
		this.problem.ub = this.problem.ub(this.I_active_variables);
	end
	options.num_orig_vars = this.NumberVariables;
	options.bCompact = true;
else
	this.problem.lb = sparse(this.NumberVariables,1);
	options.bCompact = false;
end
this.problem.num_vars = [Np; this.num_varz];
options.fmincon_opt.HessianFcn = ...
	@(x,lambda)SimpleSliceOptimizer.fcnHessian(x, lambda, this, options);
[x, fval] = this.optimize(options);

%% Additional Process
% Under the problem formulation, if node |n| is not used by path |p|($h_{np} = 0$) or path
% |p| does not use VNF |f|, the corresponding $z_{npf}$ should be zero. However, the
% optimization results of the related $z_{npf}$ can be arbitrary value, since the related
% $z_{npf}$ is not counted in the cost items. As a result, we should manually set
% the value of those related $z_{npf}$ to zero.
%
% On the other hand, too small components should be rounded.
this.temp_vars.x = x(1:Np);
if options.SlicingMethod == SlicingMethod.SingleFunction
	%% TODO
	% allocate VNF instance by order.
	slice.VNFList = 1:slice.Parent.NumberVNFs;
	this.getAs_res;
	z_np = x(Np+1:end);
	z_np = reshape(this.I_dc_path(:).*z_np, Nsn,Np);
	znpf = zeros(Nsn, Np, Nvf);
	for v = 1:Nvf
		af = slice.Parent.VNFTable.ProcessEfficiency(v);
		idx_path = this.I_path_function(:,v)==1;
		p_slice = slice.path_owner(idx_path);
		znpf(:, idx_path, v) = z_np(:, idx_path).*(af./this.alpha_f(p_slice))'; %% ISSUE: alpha_f is not defined.
	end
	this.temp_vars.z = znpf(:);
else
	this.temp_vars.z = x((Np+1):end);
	nz = Nsn*Np;
	z_index = 1:nz;
	for f = 1:Nvf
		this.temp_vars.z(z_index) = this.I_dc_path(:).*this.temp_vars.z(z_index);
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
	slice.FlowTable.Rate = this.getFlowRate;
	slice.Links.Load = this.getLinkLoad;
	slice.ServiceNodes.Load = this.getNodeLoad;
end
if nargout >= 2
	cost = this.getCost(options.PricingPolicy);   % method is overrided by subclasses.
end
end
