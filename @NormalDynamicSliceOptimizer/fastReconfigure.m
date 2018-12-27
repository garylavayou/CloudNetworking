function [profit, cost] = fastReconfigure(this, action, options)
% ITER_LIMIT event_num NUMBER_ITERS;
global computime 
if nargin <= 2
	options = Dictionary;
else
	options = Dictionary(options);
end
slice = this.hs;
theta0 = this.GLOBAL_OPTIONS.theta0;
options = setdefault(options, slice.options, {'PricingPolicy'});
if ~isfield(options, 'bEnforceReserve')
	options.bEnforceReserve = false;
end
if strcmpi(this.options.Form, 'compact')
	options.bCompact = true;
else
	options.bCompact = false;
end
options.bFinal = false;
Nf = slice.NumberFlows;
if Nf == 0
	[profit, cost] = this.handle_zero_flow(options);
	return;
end
Nvnf = slice.NumberVNFs;
Nvf = this.NumberFlowSections;   % Removed flows not involved in optimization.
Nal = this.NumberAugmentedLinks;
Nan = this.NumberAugmentedNodes;
Nsn = slice.NumberServiceNodes;
Nl = slice.NumberLinks;

this.num_vars = [...
	Nvf * Nal; ...
	Nvnf*Nf*Nsn; ...
	Nf; ...
	Nvf * Nal; ...		% tx
	Nvnf*Nf*Nsn	...		% tz
	];
num_vars = sum(this.num_vars(1:5));
num_varx = this.num_vars(1);
num_varz = this.num_vars(2);
num_varr = this.num_vars(3);
num_varxz = sum(this.num_vars(1:2));
num_varv = Nvnf*Nsn;				% VNF capacity is not variables in FSR.
num_vars_prim = sum(this.num_vars(1:3));
num_lcon_proc = size(this.As_proc,1);
num_lcon_le = num_lcon_proc + Nl + Nvnf*Nsn + 2*num_varx + 2*num_varz;
if this.options.bReserve
	if ~options.bEnforceReserve
		num_lcon_le = num_lcon_le + num_varr;
		% Nf-1: limit the flow rate, new flow is not limited
		if strcmpi(action, 'add')
			num_lcon_le = num_lcon_le - 1;		
		end
	end
	if options.bEnforceReserve
		% Fast reconfigure after slice dimensioning.
		theta = this.options.Reserve;
		num_lcon_le = num_lcon_le + 2;
		if this.options.bReserve == 2
			theta1 = min((1+theta)/2, theta0);
			num_lcon_le = num_lcon_le + Nsn;
		elseif this.options.bReserve == 3
			theta1 = min((1+theta)/2, theta0);
			theta0 = theta1;
		else
			theta1 = theta0;
		end
	end
end
%% List of constaints
%		(1) flow conservation law: Nvf*Nan, <Equality>
%   (2) flow processing requirement: Nsn*Nvnf*Nf (num_lcon_res);
%   (3) VNF instance capacity constraint: Nvnf*Nsn (num_varv);
%   (4) Link capacity constraint: Nl;
%   (5) link reconfiguration cost constraint: 2*num_varx;
%   (6) node reconfiguration cost constraint: 2*num_varz;
t1 = tic;
cprintf('Comments', 'Info:[%s] initilizing arguements ... ', calledby(0));
this.problem = Dictionary(); prbm = this.problem;
switch this.hs.Weight
	case 50
		unit = 0.02;
	case 10
		unit = 0.02;
	case 300
		unit = 0.005;
	otherwise
		unit = 1;
end
prbm.Aeq = sparse(Nan*Nvf, num_vars);  % flow conservation law
prbm.Aeq(:, 1:num_vars_prim) = [this.As_flow, sparse(Nan*Nvf, Nsn*Nvnf*Nf), -this.Ids]; 
prbm.beq = sparse(Nan*Nvf, 1);
prbm.A = sparse(num_lcon_le, num_vars);
row_offset = 0;
prbm.A(1:num_lcon_proc, 1:(num_varxz)) = [this.As_proc, -this.As_procz]; % processing reuirement
row_offset = row_offset + num_lcon_proc;
prbm.A(row_offset+(1:num_varv), num_varx+(1:num_varz)) = repmat(speye(num_varv), 1, Nf);  % VNF capacity constraint
row_offset = row_offset + num_varv;
prbm.A(row_offset+(1:Nl), 1:(num_varx)) = repmat(speye(Nl, Nal), 1, Nvf); % link capacity constraint
row_offset = row_offset + Nl;
prbm.A(row_offset+(1:num_varx), 1:num_varx) = speye(num_varx);
prbm.A(row_offset+(1:num_varx), num_vars_prim+(1:num_varx)) = -speye(num_varx);
row_offset = row_offset + num_varx;
prbm.A(row_offset+(1:num_varx), 1:num_varx) = -speye(num_varx);
prbm.A(row_offset+(1:num_varx), num_vars_prim+(1:num_varx)) = -speye(num_varx);
row_offset = row_offset + num_varx;
prbm.A(row_offset+(1:num_varz), (num_varx+1):num_varxz) = speye(num_varz);
prbm.A(row_offset+(1:num_varz), (num_vars_prim+num_varx+1):num_vars) = -speye(num_varz);
row_offset = row_offset + num_varz;
prbm.A(row_offset+(1:num_varz), (num_varx+1):num_varxz) = -speye(num_varz);
prbm.A(row_offset+(1:num_varz), (num_vars_prim+num_varx+1):num_vars) = -speye(num_varz);
row_offset = row_offset + num_varz;

cv = this.getVNFCapacity();    % VNFCapacity (Variables.v should be cleared) not change under 'fastconfig'
cl = slice.Links.Capacity;
bs0 = [ sparse(num_lcon_proc,1);
	theta0 * cv;
	theta0 * cl;
	this.topts.old_variables.x;       % which have been pre-processed, so it can be
	-this.topts.old_variables.x;      % compared with the current states.
	this.topts.old_variables.z;
	-this.topts.old_variables.z];
if this.options.bReserve
	if ~options.bEnforceReserve
		if strcmpi(action, 'add')
			%% limit the existing flows data rate
			[~, limit_idx] = sort(slice.FlowTable.Rate, 'descend');
			num_limit1 = floor(0.95*Nf);
			limit_idx1 = limit_idx(1:num_limit1);
			prbm.A(row_offset+(1:num_limit1), num_varxz+(1:num_varr)) = ...
				sparse(1:num_limit1, limit_idx1, ones(num_limit1,1), num_limit1, num_varr);
			row_offset = row_offset + num_limit1;
			%% Limit new flow/low rate flow's data rate
			% We assume the last flow in the flow table is the new flow (stable sort)
			limit_idx2 = limit_idx((num_limit1+1):(Nf-1));
			num_limit2 = numel(limit_idx2);
			prbm.A(row_offset+(1:num_limit2), num_varxz+(1:num_varr)) = ...
				sparse(1:num_limit2, limit_idx2, ones(num_limit2,1), num_limit2, num_varr);
			% new arriving flow has no limit. (The corresponding row in As is 0)
			row_offset = row_offset + length(limit_idx2); 
			prbm.b = [bs0; this.max_flow_rate(limit_idx1); ...
				max(this.max_flow_rate(limit_idx2),1.5*mean(slice.FlowTable.Rate)*ones(length(limit_idx2),1))]/unit;
		else  % action is remove
			%% limit the existing flows data rate
			[~, limit_idx] = sort(slice.FlowTable.Rate, 'descend');
			num_limit1 = floor(0.95*Nf);
			limit_idx1 = limit_idx(1:num_limit1);
			prbm.A(row_offset+(1:num_limit1), num_varxz+(1:num_varr)) = ...
				sparse(1:num_limit1, limit_idx1, ones(num_limit1,1), num_limit1, num_varr);
			row_offset = row_offset + num_limit1;
			limit_idx2 = limit_idx((num_limit1+1):Nf);
			num_limit2 = numel(limit_idx2);
			prbm.A(row_offset+(1:num_limit2), num_varxz+(1:num_varr)) = ...
				sparse(1:num_limit2, limit_idx2, ones(num_limit2,1), num_limit2, num_varr);
			row_offset = row_offset + length(limit_idx2); 
			prbm.b = [bs0; this.max_flow_rate(limit_idx1); ...
				max(this.max_flow_rate(limit_idx2),1.5*mean(slice.FlowTable.Rate)*ones(length(limit_idx2),1))]/unit;
		end
	else
		prbm.b = bs0/unit;
	end
	%% limit the resource utilization
	if options.bEnforceReserve
		% overall node / link resource utilization
		% NOTE: use 'I_active_variables' to filter the unused variables.
		% The 'I_active_variables' is only partial.
		prbm.A(row_offset+1, num_varx+(1:num_varz)) = ...
			this.I_active_variables((num_varx+1):num_varxz);
		prbm.A(row_offset+2, 1:num_varx) = ~this.I_fake_edgevars;
		row_offset = row_offset + 2;
		if this.options.bReserve == 2
			% individual node resource utilization
			prbm.A(row_offset+(1:Nsn), (num_varx+1):num_varxz) = ...
				this.As_load(Nl+(1:Nsn), (num_varx+1):num_varxz);
			row_offset = row_offset + Nsn; %#ok<NASGU>
		elseif this.options.bReserve == 3
			warning('incomplete implementation!');
		end
		%%
		% The resource reservation constraint is consistent with those in
		% <priceOptimalFlowRate>. If sole FSR need resource reservation on individual
		% resources, the constraints might need to consider the "old link/VNF load", which
		% might be higher than the "theta" items.
		prbm.b = [prbm.b; theta*sum(cv)/unit; theta*sum(cl)/unit];
		Nvn = length(cv);
		if this.options.bReserve == 2  % individual link/node
			% VNF instance capacity: as node capacity takes effect, 'theta0' is used for VNF
			% capacity.
			prbm.b(num_lcon_proc+(1:Nvn)) = theta0*cv/unit;						
			prbm.b(num_lcon_proc+Nvn+(1:Nl)) = theta1*cl/unit;						% link capacity
			prbm.b = [prbm.b; theta1*sum(reshape(cv, Nsn, Nvnf),2)/unit];% Node capacity
		elseif this.options.bReserve == 3 % % individual link/VNF instance
			% update constraints, no new ones.
			prbm.b(num_lcon_proc+(1:(Nvn+Nl))) = theta1*[cv;cl]/unit;
		end
	end
else
	prbm.b = bs0/unit;
end
% prbm.x0 = [this.topts.old_variables.x;
%   this.topts.old_variables.z;
%   this.topts.old_variables.r;
%   this.topts.old_variables.x;
%   this.topts.old_variables.z
%   ];
prbm.x0 = sparse(num_vars,1);

if ~this.checkFeasible(prbm.x0)
	warning('infeasible initial point.');
end
t2 = toc(t1);
cprintf('Comments', 'Elapsed time is %f seconds.\n', t2);

%% Perform Optimization
minopts = optimoptions(@fmincon);
minopts.Algorithm = 'interior-point';
minopts.SpecifyObjectiveGradient = true;
minopts.Display = 'iter';     % {'iter'|'notify'}
minopts.MaxIterations = 60;
% minopts.OptimalityTolerance = 1e-5;
% minopts.ScaleProblem = true;
% minopts.DiffMaxChange = 10^3;
% minopts.CheckGradients = true;
% minopts.FiniteDifferenceType = 'central';
t1 = tic;
if this.options.bDistributed
	error('error: not support.');
else
	prbm.num_vars = this.num_vars;
	if options.bCompact
		%% Get the compact formulation
		% We have determined some variables that is not used in the problem. The corresponding
		% coefficients has been set to 0. 
		% Some new constraints has been added to the problem, we need to clear the
		% coefficients of those unused (temp) variables, as to determine the active variables
		% for this problem, resulting in 'this.I_active_variables'. 
		% In addition, we need to filter out some rows of the newly added constraints, since
		% 'this.I_fake_vars' indicate that some edge variables should not be counted in
		% reconfiguration.
		
		% eleminate the columns (unused variables)
		I_active_vars = this.updateActiveVariables();  % xzr
		I_active_vars = [I_active_vars, I_active_vars(1:num_varxz)];
		prbm.A(:, ~I_active_vars) = 0; 
		% eleminate the rows (fake variables)
		% NOTE: some "active variables" are also "fake variables" that are used in the
		% flow conservation constriants and the processing requirement constraints. Therefore,
		% "fake variables" cannot be eleminated from the problem. Instead, only the
		% reconfiguration related rows are removed.
		row_offset = num_lcon_proc + num_varv + Nl;
		
		idx_unused_rows = [row_offset+find(this.topts.x_reconfig_cost==0); ...
			row_offset+num_varx+find(this.topts.x_reconfig_cost==0);...
			row_offset+num_varx*2+find(this.topts.z_reconfig_cost==0);...
			row_offset+num_varx*2+num_varz+find(this.topts.z_reconfig_cost==0)];
		prbm.A(idx_unused_rows, :) = 0;		% CLEARED
		prbm.b(idx_unused_rows, :) = 0;
		I_unused_rows = sum(abs(prbm.A), 2)==0;
		% NOTE: no need to assert. If the constraint coefficients is cleared in a row, the
		% constraint is not needed. (for VNF/Link capacity constraint, temp variables
		% constraint.) 
		% 		assert(isempty(find(prbm.b(I_unused_rows)>eps,1)),
		% 			'error: infeasible inequality constraints.');
		prbm.A(I_unused_rows, :) = [];		% REMOVED
		prbm.b(I_unused_rows, :) = [];
		I_unused_rows = sum(abs(prbm.Aeq), 2)==0;
		% NOTE: no need to assert. |beq| is all-zero.
		% 		assert(isempty(find(prbm.beq(I_unused_rows)>eps,1)), ...
		% 			'error: infeasible equality constraints.');
		prbm.Aeq(I_unused_rows, :) = [];		% REMOVED
		prbm.beq(I_unused_rows, :) = [];
		% num_lcon_le = num_lcon_le - length(idx_fake_vars);
		% determine the active variables in the problem.
		% sum(abs(prbm.A),1)~=0 should be equal to  sum(abs([prbm.Aeq; prbm.A]),1)~=0
		I_active_variables = sparse(sum(abs([prbm.A; prbm.Aeq]),1)~=0);
		prbm.A = prbm.A(:, I_active_variables);
		prbm.Aeq = prbm.Aeq(:, I_active_variables);
		prbm.x0 = prbm.x0(I_active_variables);
		prbm.num_vars = [nnz(I_active_variables(1:num_varx));
			nnz(I_active_variables((num_varx+1):num_varxz));
			nnz(I_active_variables((num_varxz+1):num_vars_prim));  % Nf
			nnz(I_active_variables(num_vars_prim+(1:num_varx)));
			nnz(I_active_variables(num_vars_prim+num_varx+(1:num_varz)))];
		assert(prbm.num_vars(3) == Nf, 'error: dimensioning of VNF vars not match.')
		% The |x_reconfig_cost| also include the items for the "fake variables", see also
		% <update_reconfig_costinfo>.
    old_opts = this.topts;
		this.topts.x_reconfig_cost = this.topts.x_reconfig_cost(...
			I_active_variables(num_vars_prim+(1:num_varx)));  
		this.topts.z_reconfig_cost = this.topts.z_reconfig_cost(...
			I_active_variables(num_vars_prim+num_varx+(1:num_varz)));
	else
		error('error: un-support.');
	end
	prbm.lb = sparse(length(prbm.x0),1);
	options.unit = unit;
	prbm.x0 = prbm.x0/unit;
	minopts.HessianFcn = ...
		@(x,lambda)NormalDynamicSliceOptimizer.hessReconfigure(x, lambda, this, options);
	[xs, fval, exitflag, output] = ...
		fmincon(@(x)NormalDynamicSliceOptimizer.fcnFastConfigProfit(x, this, options), ...
		prbm.x0, prbm.A, prbm.b, prbm.Aeq, prbm.beq, prbm.lb, [], [], minopts);
	t = 10;
  while exitflag<0 && t
		warning([strtok(output.message, newline), 'Retry.']);
		prbm.x0 = xs/unit;
		[xs, fval, exitflag, output] = ...
			fmincon(@(x)NormalDynamicSliceOptimizer.fcnFastConfigProfit(x, this, options), ...
			prbm.x0, prbm.A, prbm.b, prbm.Aeq, prbm.beq, prbm.lb, [], [], minopts);		
    t = t - 1;
	end
	xs = xs*unit;
	if options.bCompact
		x = sparse(num_vars, 1);
		x(I_active_variables) = xs;
	else
		x = xs;
	end
	% x is a local solution to the problem when exitflag is positive.
	this.interpretExitflag(exitflag, output.message);
  if options.bCompact
    this.topts = old_opts;
  end
end
t2 = toc(t1);
if ~isempty(computime)
	computime(event_num-1) = t2;
end
options.Action = action;    % This might be used when check feasible solution.
options.ConstraintTolerance = minopts.ConstraintTolerance;
assert(this.checkFeasible(x,options), 'error: infeasible solution.');

%%%
% The post processing like <optimalFlowRate> is not needed, since the
% objective function will force those variables to be zero.
this.temp_vars.x = x(1:num_varx);
this.temp_vars.z = x((num_varx+1):num_varxz);
this.temp_vars.r = x((num_varxz+1):num_vars_prim);
this.flow_rate = this.temp_vars.r;
%% Up-scaling and Reconfiguration 
% The raw solution is most close to the last stage's solution. The rounding procedure
% might lead to drop some paths, resulting in more reconfigurations. Therefore, the
% upscaling can remedy the rounded-solution to make it more close to the raw solution.
this.postProcessing();
slice.FlowTable.Rate = this.getFlowRate;
slice.Links.Load = this.getLinkLoad;
slice.ServiceNodes.Load = this.getNodeLoad;
if nargout >= 1
	cost = slice.getCost('const');
	rc_linear = this.get_reconfig_cost('linear', true);
	profit = -fval - cost + rc_linear;
end

end  % end of function