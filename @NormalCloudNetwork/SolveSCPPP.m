%% Proximal point method
% [Spingarn]
function results = SolveSCPPP(this, slices, prices, options)
%% Global Declaration
global ITER_LIMIT; 
if isempty(ITER_LIMIT)
	ITER_LIMIT = inf;
end
assert(isa(options, 'Dictionary'), ...
	'''options'' should be an output argument (type of <Dictionary>).');
if ~isfield(options, 'InterSlicePenalty')
	% 	options.InterSlicePenalty = 0.1*mean([prices.Node; prices.Link]);
	options.InterSlicePenalty = 32/length(slices)/mean([prices.Node; prices.Link]);
end
r = options.InterSlicePenalty;

%% Initialization parameters
M = getParallelInfo();
Ns = length(slices);
options.NumberSlices = Ns;
Ne = this.NumberLinks;
Ndc = this.NumberDataCenters;
num_duals = Ne + Ndc;
gamma_k = zeros(num_duals, Ns); % auxiliary variables to lambda.
q_k = zeros(num_duals, Ns); % Dual variables for the dual ADMM formulation.
x_prox = cell(Ns,1);
fval_k = zeros(Ns,1);
fval__ = 0;
output_k = cell(Ns,1);
loads = zeros(num_duals, Ns);
par_ops = SliceOptimizer.empty(Ns,0);
for si = 1:Ns
	sl = slices(si);
	sl.Optimizer.setProblem('Price', prices);
	sl.Optimizer.update_options(options);
	capacity = struct('Link', options.Capacity.Link/Ns, ...
		'Node', options.Capacity.Node/Ns);
	sl.Optimizer.setProblem('Capacity', capacity);
	if M > 0
		sl.Optimizer.Host = Slice.empty();
	end
	par_ops(si) = sl.Optimizer;
	x_prox{si} = zeros(sl.num_flow_vars,1);
end
k = 1;
t1 = tic;
while true
	%% Step 1: update the dual-variables
	lambda_k = 1/Ns*sum(gamma_k,2);
	x_prox__ = x_prox;
	
	%% Step-2: solve the sup-problems, return the resource load
	if M > 0
		options.bParallel = true;
		parfor (sj = 1:Ns,M)
			[gamma_k(:,sj), fval_k(sj), x_prox{sj}, output_k{sj}] = ...
				par_ops(sj).priceOptimalFlowRatePP([], lambda_k, q_k(:,sj), x_prox{sj}, options);
		end
	else
		for sj = 1:Ns
			[gamma_k(:,sj), fval_k(sj), x_prox{sj}, output_k{sj}] = ...
				par_ops(sj).priceOptimalFlowRatePP([], lambda_k, q_k(:,sj), x_prox{sj}, options);
		end
	end
	%% Step 3: update the auxiliary variables
	q_k0 = q_k + r*(lambda_k - gamma_k);
	q_k = q_k0 - 1/Ns*sum(q_k0,2);
	%% Step-4: stop condition test
	re_norm = norm(mean(gamma_k,2)-lambda_k);
	xe_norm = 0;
	xsum = 0;
	for sj = 1:Ns
		xe_norm = xe_norm + sum((x_prox{sj} - x_prox__{sj}).^2);
		xsum = xsum + sum((x_prox{sj}).^2);
	end
	xe_norm = sqrt(xe_norm);
	tol_primal = 10^-3 * sqrt(xsum);
	tol_dual = 10^-3 * norm(lambda_k);
	fval = sum(fval_k);
	fval_change = fval - fval__;
	num_iters = 0;
	num_funccount = 0;
	for si = 1:Ns
		num_iters = num_iters + output_k{si}.iterations;
		num_funccount = num_funccount + output_k{si}.funcCount;
	end
	num_iters = round(num_iters/Ns);
	num_funccount = round(num_funccount/Ns);
	if mod(k,20) == 1
		fprintf('                         Primal-                  Dual-                   Sub-       Sub- \n');
		fprintf('Iteration  Fval-change   optimality   Tolerance   optimality   Tolerance  Iterations Evaluations\n');
		cprintf('*text', ...
			      '————————— ————————————— ———————————— ——————————— ———————————— ——————————— —————————— ———————————\n');
	end
	fprintf('%8d   %11.4g   %10.4f   %9.4g   %10.4f   %9.4f   %9d  %9d\n',...
		k, fval_change, xe_norm, tol_primal, re_norm, tol_dual, num_iters, num_funccount);
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
	if xe_norm <= tol_primal && re_norm <= tol_dual
		b_stop = true;
	end
	if b_stop || k>=ITER_LIMIT
		break;
	end
	fval__ = fval;
	k = k + 1;
end
t2 = toc(t1);
fprintf('Partial-Inverse: elapsed time: %d\n', t2);

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
results = Dictionary();
results.lambda = mean(gamma_k, 2);  % gamma is non-negative, see also <SolveSCPDA.m>
results.LinkPrice = results.lambda(1:Ne) + prices.Link;
results.NodePrice = results.lambda(Ne+(1:Ndc)) + prices.Node;
results.x_prox = x_prox;
results.loads = loads;
results.numiters = k;
%% Capacity Distribution
% (a) Distribute the capacity according to the value of 'q';
% (b) Distribute the capacity accroding to the load;
if isfield(options, 'Capacity')
	capacities = [options.Capacity.Link; options.Capacity.Node]/Ns;
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