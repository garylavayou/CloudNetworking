%% priceOptimalFlowRate
% Override <SimpleSliceOptimizer.priceOptimalFlowRate>, considering resource
% reconfiguration cost. 
% NOTE: update price before call this function.
function [output, loads, fval] = priceOptimalFlowRate(this, x0, options) 
global DEBUG INFO;
if nargin <= 2
	options = Dictionary();
else
	options = Dictionary(options);
end
setdefault(options, struct('bParallel', false, 'bInitialize', false));
pardata = this.pardata;
Nsn = pardata.NumberServiceNodes;
Np = pardata.NumberPaths;
Nvnf = pardata.NumberVNFs;
Nl = pardata.NumberLinks;
%% No need to call parent method
% As this method can also handle the situtation, where no reconfiguration cost is
% considered.
%
% if this.options.ReconfigMethod == ReconfigMethod.Baseline ||...
%         this.options.ReconfigMethod == ReconfigMethod.DimBaseline || ...
%         ~this.options.bReserve && this.b_derive_vnf && isempty(fieldnames(this.net_changes))
%     if nargout <= 1
%         net_profit = priceOptimalFlowRate@Slice(this, x0, new_opts);
%     else
%         [net_profit, load] = priceOptimalFlowRate@Slice(this, x0, new_opts);
%     end
%     return;
% end

%% List of constaints
% (see also <DynamicSlice.fastReconfigure2>):
%   (1) flow processing requirement: NP*NV (num_lcon_res);
%   (2) VNF instance capacity constraint: NV*NN (num_varv); virtual node capacity constraint;
%   (3) link reconfiguration cost constraint: 2*NP;
%   (4) VNF allocation reconfiguration cost constraint: 2*NN*NP*NV (2*num_varz);
%   (5) VNF instance capacity reconfiguration cost constraint: 2*NV*NN;
%   (6) virtual link capacity constraint: NL (optional, resource reservation);
%   (7) lower bound of virtual link capacity (optional);
%   (8) lower bound of virtual node capacity (dimconfig2): sum capacity of VNF isntance
%       larger than the lowbound;
%   (9) lower bound of VNF instance capacity (dimconfig1);
%   (10)resource utilization lower bound for resource reservation (optional);
%   (11)node capacity constraint with resource reservation (optional)
%
%           1	  2	  3	  4	  5	  6	  7 	8 	9 	10
% InitRSV   +   +   -   -   -   +   -   -   -   +
% Redim     +   +   +   +   +   +   -   -   -   -
% RedimRSV  +   +   +   +   +   +   +  (+) (+)  -
% RedimRSV+ +   +   +   +   +   +   +  (+) (+)  +
%% Number of Variables
%   x: path; z: VNF instances assignment; v: VNF instances capacity;
%   tx;      tz;                          tv;
%   c: link capacity
%           x    z    v 	tx    tz     tv     c
% InitRSV   +    +    + 	-     -      -      +
% Redim     +    +    + 	+     +      +      +
% RedimRSV  +    +    +   +     +      +      +
% RedimRSV+ +    +    +   +     +      +      +
if options.bInitialize
	this.num_vars = [Np; Nvnf*Nsn*Np; Nvnf*Nsn];		% reconfigure x_p and z_npf, v_nf.
	theta0 = this.GLOBAL_OPTIONS.theta0;
	options = structmerge(options, ...
		getstructfields(this.options, {'Reserve'}, 'default', {theta0}), 'exclude');
	if this.b_derive_vnf || this.options.ReconfigMethod == ReconfigMethod.DimBaseline
		% Initialization
		% For the fisrt time, or for the case that no reconfigration cost is considered.
		this.invoke_method = 1;
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
	num_varxzv = sum(this.num_vars(1:3));
	num_varxz = sum(this.num_vars(1:2));
	num_varz = this.num_vars(2);
	num_varv = this.num_vars(3);
	num_lcon_res = size(this.As_res,1);
	num_lcon = num_lcon_res + num_varv;
	nnz_As = nnz(this.As_res) + (nnz(this.Hdiag)+num_varv);
	if this.invoke_method ~= 1
		this.num_vars = [this.num_vars; Np; Nvnf*Nsn*Np; Nvnf*Nsn];
		num_lcon = num_lcon + 2*Np + 2*num_varz + 2*num_varv;
		nnz_As = nnz_As + 4*Np + 4*num_varz+ 4*num_varv;
	end
	this.num_vars = [this.num_vars; Nl];
	num_lcon = num_lcon + Nl;
	nnz_As = nnz_As + nnz(this.I_edge_path) + Nl;
	num_vars = sum(this.num_vars);		% Total number of variables including temp ones.
	
	if this.invoke_method == 3 || this.invoke_method == 4
		if isfield(this.lower_bounds, 'node')
			if isfield(this.options, 'ResidualCapacity')
				node_idx = this.options.ResidualCapacity.Node-this.lower_bounds.node>1;
			else
				node_idx = true(Nsn,1);
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
		%     nnz_As = nnz_As + num_varz + num_varv + NP + NL;
		nnz_As = nnz_As + nnz(this.I_dc_path)*Nvnf + num_varv + Np + Nl;
	end
	if this.invoke_method == 1 || this.invoke_method == 4
		switch this.options.bReserve
			case 1 % only resource reservation to achieve overall utilization ratio
				theta1 = theta0;
			case 2 % resource reservation for individual resource (per node)
				theta1 = min((1+theta)/2, theta0);
				num_lcon = num_lcon + Nsn;
				nnz_As = nnz_As + nnz(this.Hrep) + num_varv;
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
	
	%% Initialize Problem
	this.problem = Dictionary(); prbm = this.problem;
	prbm.A = spalloc(num_lcon, num_vars, nnz_As);
	prbm.b = sparse(num_lcon,1);
	prbm.lb = sparse(num_vars, 1);
	%% (1) flow processing requirement
	prbm.A(1:num_lcon_res,1:num_varxz) = this.As_res;
	% b(1:num_lcon_res+(1:NL)) = sparse(num_lcon_res,1);
	row_offset = num_lcon_res;
	%% (2) VNF instance capacity constraint
	prbm.A(row_offset+(1:num_varv), Np+(1:num_varz)) = this.Hdiag;
	prbm.A(row_offset+(1:num_varv), num_varxz+(1:num_varv)) = -theta0*speye(num_varv);
	% b(row_offset+(1:num_varv)) = sparse(num_varv,1);
	row_offset = row_offset + num_varv;
	if this.invoke_method~=1
		%% (3) link reconfiguration cost constraint
		prbm.A(row_offset+(1:Np), 1:Np) = speye(Np);
		prbm.A(row_offset+(1:Np), num_varxzv+(1:Np)) = -speye(Np);
		prbm.b(row_offset+(1:Np)) = this.topts.old_variables.x; % which have been pre-processed, so it can be compared with the current states.
		row_offset = row_offset + Np;
		prbm.A(row_offset+(1:Np), 1:Np) = -speye(Np);
		prbm.A(row_offset+(1:Np), num_varxzv+(1:Np)) = -speye(Np);
		prbm.b(row_offset+(1:Np)) = -this.topts.old_variables.x;
		row_offset = row_offset + Np;
		%% (4) VNF allocation reconfiguration cost constraint
		prbm.A(row_offset+(1:num_varz), (Np+1):num_varxz) = speye(num_varz);
		prbm.A(row_offset+(1:num_varz), (num_varxzv+Np)+(1:num_varz)) = -speye(num_varz);
		prbm.b(row_offset+(1:num_varz)) = this.topts.old_variables.z;
		row_offset = row_offset + num_varz;
		prbm.A(row_offset+(1:num_varz), (Np+1):num_varxz) = -speye(num_varz);
		prbm.A(row_offset+(1:num_varz), (num_varxzv+Np)+(1:num_varz)) = -speye(num_varz);
		prbm.b(row_offset+(1:num_varz)) = -this.topts.old_variables.z;
		row_offset = row_offset + num_varz;
		%% (5) VNF instance capacity reconfiguration cost constraint
		prbm.A(row_offset+(1:num_varv), (num_varxz+1):num_varxzv) = speye(num_varv);
		prbm.A(row_offset+(1:num_varv), (num_varxzv+num_varxz+1):num_varxzv*2) = -speye(num_varv);
		prbm.b(row_offset+(1:num_varv)) = this.old_variables.v;
		row_offset = row_offset + num_varv;
		prbm.A(row_offset+(1:num_varv), (num_varxz+1):num_varxzv) = -speye(num_varv);
		prbm.A(row_offset+(1:num_varv), (num_varxzv+num_varxz+1):num_varxzv*2) = -speye(num_varv);
		prbm.b(row_offset+(1:num_varv)) = -this.old_variables.v;
		row_offset = row_offset + num_varv;
	end
	%% (6) virtual link capacity constraint
	% (optional) resource reservation for individual links
	prbm.A(row_offset+(1:Nl), 1:Np) = this.I_edge_path;
	if this.invoke_method == 1
		prbm.A(row_offset+(1:Nl), num_varxzv+(1:Nl)) = -theta1*speye(Nl);
	else
		prbm.A(row_offset+(1:Nl), 2*num_varxzv+(1:Nl)) = -theta1*speye(Nl);
	end
	% b(row_offset+(1:NL)) = sparse(NL,1);
	row_offset = row_offset + Nl;
	%% (8) lower bound of node capacity
	if this.invoke_method == 3 || this.invoke_method == 4
		if isfield(this.lower_bounds, 'node')
			reduced_eye = speye(Nsn);
			reduced_eye = reduced_eye(node_idx,:);
			prbm.A(row_offset+find(node_idx), num_varxz+(1:num_varv)) = -repmat(reduced_eye,1,Nvnf);
			prbm.b(row_offset+find(node_idx)) = -this.lower_bounds.node(node_idx);
			row_offset = row_offset + nn_bounds;
		end
	end
	
	%% (10) resource utilization lower bound for resource reservation
	if this.invoke_method == 1 || this.invoke_method == 4
		if this.options.bReserve  % overall resource utilization constraint
			% A(row_offset+1, [NP+(1:num_varz), num_varxz+(1:num_varv)]) = ...
			% 	[ones(1,num_varz), -theta*ones(1,num_varv)];
			prbm.A(row_offset+1, [Np+(1:num_varz), num_varxz+(1:num_varv)]) = ...
				[repmat((this.I_dc_path(:))', 1, Nvnf), -theta*ones(1,num_varv)];  % z_npf => (node, path, function)
			prbm.A(row_offset+2, [1:Np, (num_vars-Nl+1):num_vars]) = ...
				[sum(this.I_edge_path,1), -theta*ones(1, Nl)];
			% b(row_offset+(1:2)) = [0; 0];
			row_offset = row_offset + 2;
		end
		if this.options.bReserve == 2 % per node
			%% (11) virtual node capacity constraint for resource reservation
			% (optional) resource reservation for individual nodes
			prbm.A(row_offset+(1:Nsn), Np+(1:num_varz)) = this.Hrep;
			prbm.A(row_offset+(1:Nsn), num_varxz+(1:num_varv)) = -theta1*repmat(speye(Nsn),1,Nvnf);
			% b(row_offset+(1:NC)) = sparse(NC,1);
			row_offset = row_offset + Nsn; %#ok<NASGU>
		elseif this.options.bReserve == 3
			
		end
		% (per VNF instance is appied to constriant (2))
	end
	%% (12)
	% row_offset = row_offset + NC;
	% if this.invoke_method == 4
	%     if ~isempty(idx_up_node)
	%         d = repmat(speye(NC),1,NV);
	%         A(row_offset+(1:numel(idx_up_node)), num_varxz+(1:num_varv)) = d(idx_up_node,:);
	%         row_offset = row_offset + numel(idx_up_node);
	%         b = [b; this.upper_bounds.node(idx_up_node)];
	%     end
	%     if ~isempty(idx_up_link)
	%         A(row_offset+(1:numel(idx_up_link)), num_varxz+(1:NP)) = this.I_edge_path(idx_up_link,:);
	%         b = [b; this.upper_bounds.link(idx_up_link)];
	%     end
	% end
	%% (7)(9)
	if (this.invoke_method == 3 || this.invoke_method == 4)
		if isfield(this.options, 'ResidualCapacity')
			link_idx = this.options.ResidualCapacity.Link-this.lower_bounds.link>1;
		else
			link_idx = true(Nl,1);
		end
		prbm.lb(2*num_varxzv+find(link_idx)) = this.lower_bounds.link(link_idx);
		if isfield(this.lower_bounds, 'VNF')
			vnf_idx = repmat(node_idx, Nvnf, 1);
			prbm.lb(num_varxz+find(vnf_idx)) = this.lower_bounds.VNF(vnf_idx);
		end
	end
	
	if this.invoke_method == 1
		prbm.x0 = sparse(num_vars,1);
	else
		prbm.x0= [this.topts.old_variables.x/4;
			this.topts.old_variables.z/2;
			this.old_variables.v;
			sparse(num_varxzv,1);
			%         this.topts.old_variables.x;
			%         this.topts.old_variables.z;
			%         this.old_variables.v;
			this.old_variables.c];
	end
	assert(this.checkFeasible(prbm.x0), 'error: infeasible start point.');
	
	%% Compress
	if strcmpi(this.options.Form, 'compact')
		% column/variables: isequal(this.I_active_variables', sum(this.As_res,1)~=0)
		% row/constraints:  isequal(active_rows, find(sum(this.As_res,2)~=0))
		z_filter = sparse(repmat(logical(this.I_dc_path(:)), Nvnf, 1));
		this.I_active_variables = [true(Np,1);  z_filter;  true(num_varv,1)];
		if this.invoke_method ~= 1
			%         tx_filter = this.topts.x_reconfig_cost~=0;
			%         tz_filter = z_filter&(this.topts.z_reconfig_cost~=0);
			tx_filter = true(numel(this.topts.x_reconfig_cost),1);
			tz_filter = z_filter;
			this.I_active_variables = [this.I_active_variables;...
				tx_filter; tz_filter; true(num_varv,1)];
		end
		this.I_active_variables = [this.I_active_variables; true(Nl,1)];
		if this.invoke_method ~= 1
			row_offset = num_lcon_res + num_varv;
			active_rows = [true(row_offset,1); tx_filter; tx_filter; tz_filter; tz_filter;
				true(2*num_varv,1); true(Nl,1)];
			if (this.invoke_method == 4 || this.invoke_method == 3) && ...
					isfield(this.lower_bounds, 'node')
				active_rows = [active_rows; true(nn_bounds,1)];
			end
			if this.invoke_method == 4
				active_rows = [active_rows; true(2,1)];
				if this.options.bReserve == 2
					active_rows = [active_rows; true(Nsn,1)];
				end
			end
		else
			active_rows = true(num_lcon,1);
		end
		% NOTE that the |reconfig_cost| cosefficient is not filtered like the
		% <fastReconfigure>. The objective function <fcnProfitReconfigureSlicing> will use the
		% orginal cost coefficient.
		prbm.A = prbm.A(active_rows, this.I_active_variables);
		prbm.x0 = prbm.x0(this.I_active_variables);
		prbm.lb = prbm.lb(this.I_active_variables);
		prbm.b = prbm.b(active_rows);
		prbm.bCompact = true;
	else
		prbm.bCompact = false;
	end	
	%% Optimization options
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
	prbm.minopts = optimoptions(@fmincon);
	prbm.minopts.Algorithm = 'interior-point';
	prbm.minopts.SpecifyObjectiveGradient = true;
	prbm.minopts.OptimalityTolerance = 10^-5;
	prbm.minopts.MaxIterations = 60;
	prbm.minopts.MaxFunctionEvaluations = 180;
	prbm.minopts.Display = 'iter';   %'notify-detailed'; %'iter'; 'notify'
	% minopts.CheckGradients = true;
	% minopts.FiniteDifferenceType = 'central';
	% minopts.Diagnostics = 'on';
	% options.Form = 'normal';
	%% Save information
	prbm.Reserve = options.Reserve;
	if options.bParallel % update to optimizer
		prbm.invoke_method = this.invoke_method;
		prbm.num_full_vars = this.num_vars;
		if prbm.bCompact
			prbm.I_active_variables = this.I_active_variables; 
		end
	end
else
	prbm = this.problem(1);
	num_vars = sum(this.num_vars);
	options.Reserve = prbm.Reserve;
end
options = structupdate(options, prbm, {'bCompact'}, 'ignore');
if options.bCompact  % original number of variables in the problem
	options.num_orig_vars = sum(this.num_vars);
end
options = setdefault(options, this.pardata, {'PricingPolicy'});
options.bFinal = false;

if nargin >= 2 && ~isempty(x0)
	prbm.x0 = x0;
	assert(this.checkFeasible(prbm.x0, options), 'error: infeasible start point.');
	if options.bCompact
		prbm.x0 = prbm.x0(this.I_active_variables);
	end
end

if this.invoke_method == 1
	prbm.minopts.HessianFcn = ...
		@(x,lambda)SimpleDynamicSliceOptimizer.hessInitialSlicing(x, lambda, this, options);
	[xs, fval, exitflag, foutput] = ...
		fmincon(@(x)SimpleDynamicSliceOptimizer.fcnProfitReserveSlicing(x, this, options), ...
		prbm.x0, prbm.A, prbm.b, [], [], prbm.lb, [], [], prbm.minopts);
else
	prbm.minopts.HessianFcn = ...
		@(x,lambda)SimpleDynamicSliceOptimizer.hessSlicing(x, lambda, this, options);
	[xs, fval, exitflag, foutput] = ...
		fmincon(@(x)SimpleDynamicSliceOptimizer.fcnProfitReconfigureSlicing(x, this, options), ...
		prbm.x0, prbm.A, prbm.b, [], [], prbm.lb, [], [], prbm.minopts);
end
if options.bCompact
	x = sparse(num_vars, 1);
	x(this.I_active_variables) = xs;
else
	x = xs;
end
this.interpretExitflag(exitflag, foutput.message);
if (~isempty(DEBUG) && DEBUG) || (~isempty(INFO) && INFO)
	fprintf('\tThe optimal net profit of the slice: %G.\n', -fval);
end

%% Output solution
options.ConstraintTolerance = prbm.minopts.ConstraintTolerance;
assert(this.checkFeasible(x, options), 'error: infeasible solution.');
output = Dictionary();
output.iterations = foutput.iterations;
output.funcCount = foutput.funcCount;
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
	if options.bInitialize
		prbm = Dictionary();
		theta = options.Reserve;
		loads.Node = this.getNodeLoad(false);
		loads.Link = this.getLinkLoad(false);
		prbm.x0r = [this.temp_vars.c;
			this.old_variables.v;
			sparse(Nsn*Nvnf,1)];
		prbm.lbr = [max(this.lower_bounds.link, loads.Link);      % options.bDisableLowerBound
			this.Hdiag*this.temp_vars.z;
			sparse(Nsn*Nvnf,1)];
		%%
		num_lcons = Nsn + 2 + 2*Nsn*Nvnf;
		row_offset = 0;
		prbm.Ar = sparse(num_lcons, num_vars);
		prbm.br = sparse(num_lcons, 1);
		prbm.br(1:Nsn) = -max(this.lower_bounds.node, loads.Node);    % options.bDisableLowerBound
		prbm.Ar(1:Nsn,Nl+(1:Nsn*Nvnf)) = -repmat(speye(Nsn),1,Nvnf);
		row_offset = row_offset + Nsn;
		prbm.br(row_offset+1) = -sum(loads.Link)/theta;
		prbm.Ar(row_offset+1,1:Nl) = -ones(1,Nl);
		row_offset = row_offset + 1;
		prbm.br(row_offset+1) = -sum(loads.Node)/theta;
		prbm.Ar(row_offset+1,Nl+(1:Nsn*Nvnf)) = -ones(1,Nsn*Nvnf);
		row_offset = row_offset + 1;
		prbm.Ar(row_offset+(1:Nsn*Nvnf), Nl+(1:2*Nsn*Nvnf)) = [speye(Nsn*Nvnf), -speye(Nsn*Nvnf)];
		prbm.br(row_offset+(1:Nsn*Nvnf)) = this.old_variables.v;
		row_offset = row_offset + Nsn*Nvnf;
		prbm.Ar(row_offset+(1:Nsn*Nvnf), Nl+(1:2*Nsn*Nvnf)) = [-speye(Nsn*Nvnf), -speye(Nsn*Nvnf)];
		prbm.br(row_offset+(1:Nsn*Nvnf)) = -this.old_variables.v;
		%     row_offset = row_offset + Nsn*Nvnf;
		%%
		% we might add a resource utilization ratio for individual resources.
		prbm.minopts = optimoptions(@fmincon);
		prbm.minopts.Algorithm = 'interior-point';
		prbm.minopts.SpecifyObjectiveGradient = true;
		prbm.minopts.Display = 'iter';   %'notify-detailed'; %'iter'; 'notify'
		%     fr_opt.CheckGradients = true;
		%     fr_opt.FiniteDifferenceType = 'central';
		
		%%
		% we might add a resource utilization ratio for individual resources.
		this.problem(2) = prbm;
	else
		prbm = this.problem(2);
	end
	prbm.minopts.HessianFcn = @(x,lambda)hessReconfigCost(x, lambda, this, options);
	[x, ~, exitflag, foutput] = ...
		fmincon(@(x)fcnReconfigCost(x, this, options), ...
		prbm.x0r, prbm.Ar, prbm.br, [], [], prbm.lbr, [], [], prbm.minopts);
	this.interpretExitflag(exitflag, foutput.message);
	this.temp_vars.c = x(1:Nl);
	this.temp_vars.v = x(Nl+(1:Nsn*Nvnf));
	this.temp_vars.tv = x(Nl+Nsn*Nvnf+(1:Nsn*Nvnf));
	output.iterations = output.iterations + foutput.iterations;
	output.funcCount = output.funcCount + foutput.funcCount;
end

if options.bParallel
	output.x0 = x;
	output.temp_vars = Dictionary(this.temp_vars);
	if options.bInitialize
		output.problem = prbm;
	end
end
output.net_profit = -fval;
end

%% TEST
%{
tx_index = (1:NP) + num_base_vars;
tx = x(tx_index);
load.Node = this.getNodeLoad(false);
node_capacity = this.getNodeCapacity(false);
disp([load.Node, node_capacity])
disp(load.Node./node_capacity)
load.Link = this.getLinkLoad(false);
link_capacity = this.getLinkCapacity(false);
disp('link load | link capacity |  old link capacity');
disp([load.Link, link_capacity, this.old_state.link_capacity, ...
    this.lower_bounds.link, this.upper_bounds.link]);
sum(load.Node)./sum(node_capacity)
sum(load.Link)./sum(link_capacity)
%}

function [f,g] = fcnReconfigCost(vars, this, options)
Nsn = this.pardata.NumberServiceNodes;
Nl = this.pardata.NumberLinks;
Nvnf = this.pardata.NumberVNFs;
link_cap = vars(1:Nl);
vnf_cap = vars(Nl+(1:Nsn*Nvnf));
vnf_cap_diff = vars(Nl+Nsn*Nvnf+(1:Nsn*Nvnf));
load.Node = sum(reshape(vnf_cap, Nsn, Nvnf),2);
switch options.PricingPolicy
	case {'quadratic-price', 'quadratic'}
		[link_payment,link_price_grad] = this.fcnLinkPricing(this.prices.Link, link_cap);
		[node_payment,node_price_grad] = this.fcnNodePricing(this.prices.Node, load.Node);
		f = link_payment + node_payment;
	otherwise
		error('%s: invalid pricing policy', calledby);
end
f = f + dot(vnf_cap_diff, this.topts.vnf_reconfig_cost);

g = [link_price_grad; repmat(node_price_grad, Nvnf, 1); this.topts.vnf_reconfig_cost];

end

function h = hessReconfigCost(vars, lambda, this, options) %#ok<INUSL>
Nsn = this.pardata.NumberServiceNodes;
Nl = this.pardata.NumberLinks;
Nvnf = this.pardata.NumberVNFs;
link_cap = vars(1:Nl);
vnf_cap = vars(Nl+(1:Nsn*Nvnf));
load.Node = sum(reshape(vnf_cap, Nsn, Nvnf),2);
num_vars = length(vars);
h = spalloc(num_vars,num_vars, num_vars);
switch options.PricingPolicy
	case {'quadratic-price', 'quadratic'}
		[~,~,lph] = this.fcnLinkPricing(this.prices.Link, link_cap);
		[~,~,nph] = this.fcnNodePricing(this.prices.Node, load.Node);
		h(1:Nl,1:Nl) = diag(lph);
		h(Nl+(1:Nsn*Nvnf),Nl+(1:Nsn*Nvnf)) = block_diag(diag(nph), Nvnf);
	otherwise
		error('%s: invalid pricing policy', calledby);
end
end