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
function [profit, cost, output] = optimalFlowRate( this, options )
if nargin <= 1
	options = Dictionary();
else
	options = Dictionary(options);
end
setdefault(options, struct('bParallel', false, 'bInitialize', false, ...
	'isFinalize', false));
if isempty(this.pardata)
	this.initializeParallel('optimalFlowRate', options);
end
pardata = this.pardata;
Nsn = pardata.NumberServiceNodes;
Np = pardata.NumberPaths;
Nvnf = pardata.NumberVNFs;
if options.bInitialize
	Nl = pardata.NumberLinks;
	nnz_As = nnz(this.As_res) + Nvnf*nnz(this.I_dc_path) + nnz(this.I_edge_path);
	num_lcon_res = size(this.As_res,1);
	num_lcon = num_lcon_res + Nsn + Nl;
	
	this.num_vars = [Np; Nvnf*Nsn*Np];
	num_vars = sum(this.num_vars(1:2));
	num_varz = this.num_vars(2);
	
	this.problem = Dictionary(); prbm = this.problem;
	prbm.As = spalloc(num_lcon, num_vars, nnz_As);
	prbm.As(1:num_lcon_res,:) = this.As_res;
	%% Add node and link capacity constraint coefficient
	row_index = num_lcon_res + (1:Nsn);
	col_index = Np + (1:num_varz);
	prbm.As(row_index, col_index) = this.Hrep;
	row_index = (1:Nl) + num_lcon_res + Nsn;
	prbm.As(row_index, 1:Np) = this.I_edge_path;
	%% Boundary
	% * *Resource Constraints*: processing-rate, the right side is all-zero; capacity
	% constraints, the right side is the capacity of link and node;
	prbm.bs = [sparse(num_lcon_res,1);
		pardata.NodeCapacity;
		pardata.LinkCapacity];
	%%%
	% * *Upper Bound*: Not necessary, to facilitate the algorithm, we give a relaxed
	% upper-bound.
	% max(this.Parent.readLink('Capacity'))*ones(NP,1);...
	%     max(this.Parent.readDataCenter('Capacity'))*ones(num_vars-NP,1)
	prbm.ub = [];
	%%%
	% * Remove the capacity constraints with infinity capacity;
	idx = find(prbm.bs==inf);
	prbm.As(idx,:) = [];
	prbm.bs(idx) = [];
	%% Feasible test of start point
	% * *Start Point*: in case that the capacity of a virtual link/node is zero, we initialize
	% $z_{min}$ and $x_{min}$ as the nonzero minimum value.
	prbm.x0 = zeros(num_vars,1);
	if pardata.SlicingMethod == SlicingMethod.SingleFunction
		max_alpha_f = max(options.Alpha_f);
	else
		max_alpha_f = max(pardata.ProcessEfficiency);
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
	prbm.x0(1:Np) = 1;
	prbm.x0((Np+1):(Np+num_varz)) = max_alpha_f;
	assert(this.checkFeasible(prbm.x0), 'error: infeasible start point.');
	prbm.Aeq = [];
	prbm.beq = [];
	
	%% Compress
	if strcmpi(this.options.Form, 'compact')
		z_filter = sparse(repmat(...
			reshape(logical(this.I_dc_path), numel(this.I_dc_path),1), Nvnf, 1));
		this.I_active_variables = [true(Np,1) ;  z_filter];
		prbm.As = prbm.As(:, this.I_active_variables);
		prbm.x0 = prbm.x0(this.I_active_variables);
		prbm.lb = sparse(length(prbm.x0),1);
		if ~isempty(prbm.ub)
			prbm.ub = prbm.ub(this.I_active_variables);
		end
		prbm.bCompact = true;
	else
		prbm.lb = sparse(num_vars,1);
		prbm.bCompact = false;
	end
	
	%% Optimization options
	prbm.minopts.Algorithm = 'interior-point';
	prbm.minopts.SpecifyObjectiveGradient = true;
	prbm.minopts.Display = 'iter';
	% options.minopts.CheckGradients = true;
	% options.minopts.FiniteDifferenceType = 'central';
	%% Save information
	if options.bParallel
		prbm.num_full_vars = this.num_vars;
		if prbm.bCompact
			prbm.I_active_variables = this.I_active_variables;
		end
	end
else
	prbm = this.problem;
	num_varz = this.num_vars(2);
end
options = structupdate(options, prbm, {'bCompact'}, 'ignore');
if isfield(options, 'bCompact') && options.bCompact
	options.num_orig_vars = sum(this.num_vars);
end
options = setdefault(options, this.pardata, {'PricingPolicy'});
options.bFinal = false;

prbm.minopts.HessianFcn = ...
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
num_varx = this.num_vars(1);
this.temp_vars.x = x(1:num_varx);
if pardata.SlicingMethod == SlicingMethod.SingleFunction
	%% TODO
	% allocate VNF instance by order.
	slice = this.hs;
	slice.VNFList = 1:pardata.ParentNumberVNFs;
	z_np = x(num_varx+1:end);
	z_np = reshape(this.I_dc_path(:).*z_np, Nsn, Np);
	znpf = zeros(Nsn, Np, Nvnf);
	for v = 1:Nvnf
		af = slice.Parent.VNFTable.ProcessEfficiency(v);
		idx_path = this.I_path_function(:,v)==1;
		p_slice = slice.path_owner(idx_path);
		znpf(:, idx_path, v) = z_np(:, idx_path).*(af./this.alpha_f(p_slice))'; %% ISSUE: alpha_f is not defined.
	end
	this.temp_vars.z = znpf(:);
else
	this.temp_vars.z = x((num_varx+1):(num_varx+num_varz));
	nz = Nsn*Np;
	z_index = 1:nz;
	for f = 1:Nvnf
		this.temp_vars.z(z_index) = this.I_dc_path(:).*this.temp_vars.z(z_index);
		z_index = z_index + nz;
	end
end
clear znpf;
profit = -fval;
if options.bParallel && (options.isFinalize || nargout >= 3)
	error('error: Slice is unavailable in parallel mode, cannot calculate the cost.');
end
if options.bParallel
	if nargout <= 2
		warning('results should be pass out.');
	end
	output = Dictionary();
	output.x0 = x;
	output.temp_vars = this.temp_vars;
	if options.bInitialize
		output.problem = prbm;
	end
else
	slice = this.hs; % when parallel is enabled, slice is invalid (set to empty).
	if options.isFinalize
		this.postProcessing();
		slice.FlowTable.Rate = this.getFlowRate();
		slice.Links.Load = this.getLinkLoad();
		slice.ServiceNodes.Load = this.getNodeLoad();
		pardata.erase();
		% For normal mode, the erasure operation ensures that the next computation will do
		% initialization of pardata again. The operation is performed only when the procedure
		% is specified with 'finalize', indicating the computation is to end.
		%
		% For parallel mode, the erasure operation will not effect the pardata in the master
		% process. Therefore, no need to erase the data in the worker.
	end
	if nargout >= 2
		cost = slice.getCost();   % get resource cost, method can be overrided by subclasses.
	end
end
end
