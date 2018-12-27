function [profit, cost, output] = optimalFlowRate(this, options)		% see also <SliceOptimizer>
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
if options.bInitialize
	Nl = pardata.NumberLinks;
	Nal = this.NumberAugmentedLinks;
	Nvf = this.NumberFlowSections;
	Nvnf = pardata.NumberVNFs;
	Nf = pardata.NumberFlows;
	Nsn = pardata.NumberServiceNodes;
	num_lcon_proc = size(this.As_proc, 1);
	num_lcon_flow = size(this.As_flow,1);
	
	this.num_vars = [Nvf*Nal; Nvnf*Nf*Nsn; Nf; Nvnf*Nsn]; % MUST be reset here in the continue updating process
	num_varx = this.num_vars(1);
	num_varz = this.num_vars(2);
	num_varxz = sum(this.num_vars(1:2));
	
	%% Initialize Problem
	switch pardata.Weight
		case 50
			unit = 0.02;
		case 10
			unit = 0.05;
		case 300
			unit = 0.005;
		otherwise
			unit = 1;
	end
	this.problem = Dictionary('bCompact', true); prbm = this.problem;
	prbm.Aeq = [this.As_flow, sparse(num_lcon_flow, num_varz), -this.Ids];
	prbm.A = [this.As_proc, -this.As_procz, sparse(num_lcon_proc, Nf)];
	prbm.beq = sparse(num_lcon_flow, 1);
	prbm.b = sparse(num_lcon_proc, 1);
	%% Constraints
	% Add node and link capacity constraint coefficient
	% Capacity Constraints: the right side is the capacity of link and node;
	prbm.A(num_lcon_proc+(1 : Nl+Nsn), 1:num_varxz) = this.As_load;
	prbm.b(num_lcon_proc+(1 : Nl+Nsn)) = [pardata.LinkCapacity; pardata.NodeCapacity]/unit;
	%% Compact
	this.updateActiveVariables(); % this is valid for later procedure to identify active variables.
	prbm.Aeq = prbm.Aeq(:, this.I_active_variables);
	prbm.A = prbm.A(:, this.I_active_variables);
	I_effect_rows = sum(abs(prbm.Aeq), 2)~=0;
	prbm.Aeq = prbm.Aeq(I_effect_rows,:);
	prbm.beq = prbm.beq(I_effect_rows,:);
	I_effect_rows = sum(abs(prbm.A), 2)~=0;
	prbm.A = prbm.A(I_effect_rows,:);
	prbm.b = prbm.b(I_effect_rows,:);
	
	prbm.num_vars = [...
		nnz(this.I_active_variables(1:num_varx));...
		nnz(this.I_active_variables(num_varx+1 : num_varxz));...
		Nf];  % this.num_vars is the original number of variables
	prbm.x0 = sparse(sum(prbm.num_vars),1);		% /unit
	prbm.lb = sparse(sum(prbm.num_vars),1);		% /unit
	%% [TODO] Feasible test of start point
	% * *Start Point*: in case that the capacity of a virtual link/node is zero, we initialize
	% $z_{min}$ and $x_{min}$ as the nonzero minimum value.
	%%%
	% * Remove the capacity constraints with infinity capacity;
	idx = find(prbm.b==inf);
	prbm.A(idx,:) = [];
	prbm.b(idx) = [];
	
	%% Boundary
	% * *Upper Bound*: Not necessary, to facilitate the algorithm, we give a relaxed
	% upper-bound.
	% max(this.Parent.readLink('Capacity'))*ones(NP,1);...
	%     max(this.Parent.readDataCenter('Capacity'))*ones(num_vars-NP,1)
	prbm.ub = [];			% /unit
	
	%% Optimization options
	prbm.minopts = optimoptions(@fmincon);
	prbm.minopts.Algorithm = 'interior-point';
	prbm.minopts.SpecifyObjectiveGradient = true;
	prbm.minopts.Display = 'iter';
	prbm.minopts.MaxIterations = 60;
% 	prbm.minopts.CheckGradients = true;
% 	prbm.minopts.FiniteDifferenceType = 'central';
	%% Save innformation
	prbm.unit = unit;
	if options.bParallel
		prbm.num_full_vars = this.num_vars;
		prbm.I_active_variables = this.I_active_variables;
	end
else
	prbm = this.problem;
	num_varx = this.num_vars(1);
	num_varz = this.num_vars(2);
end
options = structupdate(options, prbm, {'bCompact', 'unit'});
options = setdefault(options, this.pardata, {'PricingPolicy'});
options.bFinal = false;

prbm.minopts.HessianFcn = ...
	@(x,lambda)NormalSliceOptimizer.fcnHessianCC(x, lambda, this);  % No options
[x, fval] = this.optimize(options);
if isfield(options, 'unit')
	x = x*options.unit;
end
this.x0 = x;
this.temp_vars.x = x(1:num_varx);
this.temp_vars.z = x(num_varx+(1:num_varz));
num_varr = this.num_vars(3);
this.temp_vars.r = x(num_varxz+(1:num_varr));

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
	slice = this.hs; % when parallel is enabled, slice is set to empty.
	if options.isFinalize
		% Wether performing up-scaling or not does not influence the number of
		% reconfiguration. Hence to better utilize resources, we choose to up-scale.
		% See also <fastReconfigure>.
		this.postProcessing();
		slice.FlowTable.Rate = this.getFlowRate();
		slice.Links.Load = this.getLinkLoad();
		slice.ServiceNodes.Load = this.getNodeLoad();
		pardata.erase();   % see also <SimpleSliceOptimizer>.<optimalFlowRate>.
	end
	if nargout >= 2
		cost = slice.getCost();   % get resource cost, method can be overrided by subclasses.
	end
end
end