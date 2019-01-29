%% Subroutine: SCP
% Solve the SCP with Dual-ADMM.
function [sp_profit, b_violate, results] = SolveSCPDA(this, slices, prices, options)
%% Global Declaration
global ITER_LIMIT;
if isempty(ITER_LIMIT)
	ITER_LIMIT = inf;
end
assert(isa(options, 'Dictionary'), ...
	'''options'' should be an output argument (type of <Dictionary>).');
prices = this.convertParameter(prices, 'struct');
if options.bCountTime
	t0 = tic; prt = 0; srt = 0; 
end

%% Initialization parameters
Ne = this.NumberLinks;
Ndc = this.NumberDataCenters;
M = getParallelInfo();
Ns = length(slices);
op = this.Optimizer;
r = max(1,mean([prices.Link; prices.Node]));		% op.options.InterSlicePenalty | {0.1|0.8|1} | mean(prices)
options.InterSlicePenalty = r;
eps_rel = op.options.RelativeTolerance;
eps_abs = op.options.AbsoluteTolerance;
opt_order = op.options.OptimizeOrder;
%% Initialization variables
output_k = cell(Ns,1);
num_duals = Ne + Ndc;
par_ops = SliceOptimizer.empty(Ns,0);
for si = 1:Ns
	sl = slices(si);
	sl.Optimizer.setProblem('Price', prices);
	sl.Optimizer.update_options(options);
	capacity = this.convertParameter(options.Capacity/Ns);
	sl.Optimizer.setProblem('Capacity', capacity);
	if M > 0
		sl.Optimizer.Host = Slice.empty();
	end
	par_ops(si) = sl.Optimizer;
end
gamma_k = zeros(num_duals, Ns); % auxiliary variables to lambda.
q_k = zeros(num_duals, Ns);     % Dual variables for the dual ADMM formulation.
% 		if ~isempty(this.init_gamma_k)
% 			gamma_k = repmat(this.init_gamma_k, 1, Ns);
% 			q_k = repmat(this.init_q_k, 1, Ns);
% 		else
% 			gamma_k = ones(num_duals, Ns); % auxiliary variables to lambda.
% 			q_k = zeros(num_duals, Ns); % Dual variables for the dual ADMM formulation.
% 		end
fval_lambda_k = 0;
fval_gamma_k = zeros(Ns,1);
if options.bCountTime
	t1 = toc(t0);	prt = prt + t1; srt = srt + t1; t0 = tic;
end
%% Solve the dual problem by ADMM
lambda_k = (sum(gamma_k,2) - sum(q_k,2)/r)/Ns;
k = 1;
while true
	%% Step 1: update the dual-variables
	fval_gamma_k__ = fval_gamma_k;
	fval_lambda_k__ = fval_lambda_k;
	gamma_k__ = gamma_k;
	lambda_k__ = lambda_k;
	%% Step 2: solve the sub-problems, returning the primal variables
	if M > 0
		options.bParallel = true;
		parfor (sj = 1:Ns,M)
			[gamma_k(:,sj), fval_gamma_k(sj), output_k{sj}] = ...
				par_ops(sj).priceOptimalFlowRateDA([], lambda_k, q_k(:,sj), options);
		end
	else
		for sj = 1:Ns
			[gamma_k(:,sj), fval_gamma_k(sj), output_k{sj}] = ...
				par_ops(sj).priceOptimalFlowRateDA([], lambda_k, q_k(:,sj), options);
		end
	end
	lambda_k = (sum(gamma_k,2) - sum(q_k,2)/r)/Ns;  % The next iteration's λ_k
	if opt_order == 1
		lambda = lambda_k__;
	else
		lambda = lambda_k;
	end
	%% Step 3: update the auxiliary variables
	fval_lambda_k = sum(sum(q_k.*gamma_k));
	q_k = q_k + r*(lambda-gamma_k);
	%% Step 4: test stop condition
	re = gamma_k - lambda;
	re_norm = norm(re(:));
	tol_primal = eps_abs*sqrt(num_duals*Ns) + ...
		eps_rel*max(sqrt(Ns)*norm(lambda), norm(gamma_k(:)));
	if opt_order == 1
		se = r*sum(gamma_k-gamma_k__,2);
		se_norm = norm(se);
		tol_dual = eps_abs*sqrt(num_duals) + eps_rel*norm(sum(q_k,2));       % eps_rel*|A'*y|_2
	else
		se = r*(lambda_k-lambda_k__);
		se_norm = sqrt(Ns)*norm(se);
		tol_dual = eps_abs*sqrt(num_duals*Ns) + eps_rel*norm(q_k);       % eps_rel*|A'*y|_2
	end
	fval_change = (sum(fval_gamma_k)+fval_lambda_k)-(sum(fval_gamma_k__)+fval_lambda_k__);
	num_iters = 0;
	num_funccount = 0;
	for si = 1:Ns
		num_iters = num_iters + output_k{si}.iterations;
		num_funccount = num_funccount + output_k{si}.funcCount;
	end
	num_iters = round(num_iters/Ns);
	num_funccount = round(num_funccount/Ns);
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
	
	if M>0 && isfield(options, 'bInitialize') && options.bInitialize
		for sj = 1:Ns
			sl = slices(sj);
			sl.Optimizer.setProblem('Problem', output_k{si}.problem);
		end
	end
	% Update this option, so that if SolveSCPDA will be called multiple times, we only
	% need initialize the problem once, providing that only prices change.
	options.bInitialize = false;
	b_stop = false;
  % 	if re_norm < tol_primal && se_norm < tol_dual %&& abs(fval_change)< 100 %
  % 		b_stop = true;
  % 	end
	if b_stop || k>=ITER_LIMIT || re_norm == 0
		break;
	end
	
	k = k + 1;
end
if options.bCountTime
	t1 = toc(t0);
	if M > 0
		prt = prt + t1*M/Ns; srt = srt + t1*M;
	else
		prt = prt + t1/Ns; srt = srt + t1;
	end
	t0 = tic;
	fprintf('Dual-ADMM: elapsed time: %d\n', t1);
end
%%

loads = zeros(num_duals, Ns);
for si = 1:Ns
	sl = slices(si);
	sl.Optimizer.setProblem('Price', []);
	if M>0
		sl.Optimizer.Host = sl; % re-connect slice with optimizer
		sl.Optimizer.saveTempResults(output_k{si});
	end
	idx = [sl.Links.PhysicalLink; Ne+sl.getDCPI()];
	loads(idx, si) = [sl.getLinkLoad(false); sl.getNodeLoad(false)];
end
% this.op.init_gamma_k = mean(gamma_k, 2);
% this.op.init_q_k = mean(q_k, 2);
results = Dictionary();
%% Determine the final price
% As the procedure is not precisely terminated, we need to approximated the dual
% variables.
% higher price:
% results.lambda = max([gamma_k, lambda_k], [], 2);
% mean price:
results.lambda = mean(gamma_k, 2);  % gamma is non-negative
% lower price:
% results.lambda = max(lambda_k,0);
results.Prices = results.lambda + [prices.Link; prices.Node];
results.loads = loads;
results.numiters = k;
results.runtime = toc(t0);
%% Capacity Distribution
% (a) Distribute the capacity according to the value of 'q';
%			'q' is a compensation to slice capacity. Negative 'q' means that capacity violation
%			is compensated by 'q', like that the capacity allocated to the slice is increased by
%			'q', vice versa. See <Fukushima1992, Sec.4>.
% (b) Distribute the capacity accroding to the load;
if isfield(options, 'Capacity')
	capacities = options.Capacity/Ns;
	delta_capacity = zeros(size(q_k));
	for i = 1:num_duals
		%% capacity that can be transfered
		% 'q' is positive and smaller than the pre-allocated capacity.
		delta_capacity(i,:) = min(q_k(i,:), capacities(i));
		% count the transfer demand (q_k<0) and the available offer (0<q_k<capacities)
		spare_capacity = sum(delta_capacity(i,:));
		if spare_capacity < 0
			% when the available resource is limited, the resource is allocated propotinoally.
			b_transfer = delta_capacity(i,:)>=0;
			delta_capacity(i,~b_transfer) = delta_capacity(i,~b_transfer)/...
				(-sum(delta_capacity(i,~b_transfer))/sum(delta_capacity(i,b_transfer)));
		end
	end
	results.Capacities = capacities - delta_capacity;
end
sp_profit = this.getSliceProviderProfit(slices, prices, ...
  getstructfields(options, {'PricingPolicy'}, 'error'));
if ~options.CapacityConstrained
	results.loads = this.convertParameter(this.getNetworkLoad(slices, options));
	results.loads(results.loads<this.op.options.NonzeroTolerance) = 0;
	results.violates = results.loads > [options.ResidualCapacity.Link; ...
		options.ResidualCapacity.Node];
	b_violate = ~isempty(find(results.violates,1));
else
	b_violate = false;
end
if options.bCountTime
	t1 = toc(t0); prt = prt + t1; srt = srt + t1;
	results.runtime = struct('Parallel', prt, 'Serial', srt);
end
end
