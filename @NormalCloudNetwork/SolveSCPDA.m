%% Subroutine: SCP
% Solve the SCP with Dual-ADMM.
function results = SolveSCPDA(this, slices, prices, options)
%% Global Declaration
global ITER_LIMIT; 
if isempty(ITER_LIMIT)
	ITER_LIMIT = inf;
end
%% Initialization parameters
num_process = length(slices);
r = this.options.InterSlicePenalty;
eps_rel = this.options.RelativeTolerance;
eps_abs = this.options.AbsoluteTolerance;
opt_order = this.options.OptimizeOrder;
M = getParallelInfo();
%% Initialization variables
Ns = length(slices);
num_duals = this.NumberLinks + this.NumberDataCenters;
num_links = this.NumberLinks;
for si = 1:num_process
	sl = slices(si);
	dc_id = sl.getDCPI;
	sl.prices.Link = link_price;
	sl.prices.Link = sl.prices.Link(sl.VirtualLinks.PhysicalLink);
	sl.prices.Node = node_price;
	sl.prices.Node = sl.prices.Node(dc_id);
end
gamma_k = zeros(num_duals, num_process); % auxiliary variables to lambda.
q_k = zeros(num_duals, num_process); % Dual variables for the dual ADMM formulation.
% 		if ~isempty(this.init_gamma_k)
% 			gamma_k = repmat(this.init_gamma_k, 1, num_process);
% 			q_k = repmat(this.init_q_k, 1, num_process);
% 		else
% 			gamma_k = ones(num_duals, num_process); % auxiliary variables to lambda.
% 			q_k = zeros(num_duals, num_process); % Dual variables for the dual ADMM formulation.
% 		end
fval_lambda_k = 0;
fval_gamma_k = zeros(num_process,1);
output_k = cell(num_process,1);

%% Solve the dual problem by ADMM
if opt_order~=1
	lambda_k = (sum(gamma_k,2) - sum(q_k,2)/r)/num_process;
end
k = 1;
t1 = tic;
while true
	%% Step 1: update the dual-variables
	fval_gamma_kminus = fval_gamma_k;
	fval_lambda_kminus = fval_lambda_k;
	if opt_order == 1
		gamma_kminus = gamma_k;
		lambda_k = (sum(gamma_k,2) - sum(q_k,2)/r)/num_process;
	end
	%% Step 2: solve the sub-problems, returning the primal variables
	parfor (sj = 1:num_process,M)
		% 			for sj = 1:num_process
		sl = slices(sj);
		dc_id = sl.getDCPI;
		idx = [sl.VirtualLinks.PhysicalLink; num_links+dc_id];
		lambda = lambda_k;
		lambda = lambda(idx);
		q = q_k(:,sj);
		q = q(idx);
		[gamma_k(:,sj), fval_gamma_k(sj), output_k{sj}] = ...
			sl.priceOptimalFlowRate0([], lambda, q, options);
	end
	if opt_order ~= 1
		lambda_kminus = lambda_k;
		lambda_k = (sum(gamma_k,2) - sum(q_k,2)/r)/num_process;
	end
	fval_lambda_k = sum(sum(q_k.*gamma_k));
	%% Step 3: update the auxiliary variables
	q_k = q_k + r*(lambda_k-gamma_k);
	%% Step 4: test stop condition
	re = gamma_k - lambda_k;
	re_norm = norm(re(:));
	if opt_order == 1
		se = r*sum(gamma_k-gamma_kminus,2);
		se_norm = norm(se);
		tol_dual = eps_abs*sqrt(num_duals) + eps_rel*norm(sum(q_k,2));       % eps_rel*|A'*y|_2
	else
		se = r*(lambda_k-lambda_kminus);
		se_norm = sqrt(num_process)*norm(se);
		tol_dual = eps_abs*sqrt(num_duals*num_process) + eps_rel*norm(q_k);       % eps_rel*|A'*y|_2
	end
	tol_primal = eps_abs*sqrt(num_duals*num_process) + ...
		eps_rel*max(sqrt(num_process)*norm(lambda_k), norm(gamma_k(:)));
	fval_change = (sum(fval_gamma_k)+fval_lambda_k)-(sum(fval_gamma_kminus)+fval_lambda_kminus);
	num_iters = 0;
	num_funccount = 0;
	for si = 1:num_process
		num_iters = num_iters + output_k{si}.iterations;
		num_funccount = num_funccount + output_k{si}.funcCount;
	end
	num_iters = round(num_iters/num_process);
	num_funccount = round(num_funccount/num_process);
	if mod(k,20) == 1
		fprintf('                                     Primal-                  Dual-                   Sub-       Sub- \n');
		fprintf('Iteration Step-length  Dual-change   optimality   Tolerance   optimality   Tolerance  Iterations Evaluations\n');
		cprintf('*text', ...
			      '————————— ——————————— ————————————— ———————————— ——————————— ———————————— ——————————— —————————— ———————————\n');
	end
	fprintf('%8d  %10.2f   %11.4g   %10.4g  %10.4g   %10.4g  %10.4g  %8d    %8d  \n',...
		k, r, fval_change, re_norm, tol_primal, se_norm, tol_dual, num_iters, num_funccount);
	if mod(k,20) == 0
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
%%

loads = zeros(num_duals, num_process);
for si = 1:num_process
	sl = slices(si);
	sl.op.saveTempResults(output_k{si});
	loads(:,si) = output_k{si}.loads;
end
this.init_gamma_k = mean(gamma_k, 2);
this.init_q_k = mean(q_k, 2);
results.LinkPrice = lambda_k(1:this.NumberLinks) + link_price;
results.NodePrice = lambda_k(this.NumberLinks+(1:this.NumberDataCenters)) + node_price;
results.lambda = max([gamma_k, lambda_k], [], 2);
% results.lambda = max(lambda_k,0);
results.loads = loads;
results.numiters = k;
%% Capacity Distribution
% (a) Distribute the capacity according to the value of 'q';
% (b) Distribute the capacity accroding to the load;
if isfield(options, 'Capacities')
	capacities = options.Capacities/Ns;
	delta_capacity = zeros(size(q_k));
	for i = 1:num_duals
		delta_capacity(i,:) = min(q_k(i,:), capacities(i));
		spare_capacity = sum(delta_capacity(i,:));
		if spare_capacity < 0
			b_transfer = delta_capacity(i,:)>=0;
			delta_capacity(i,~b_transfer) = delta_capacity(i,~b_transfer)/...
				(-sum(delta_capacity(i,~b_transfer))/sum(delta_capacity(i,b_transfer)));
		end
	end
	results.Capacities = capacities - delta_capacity;
end
end
