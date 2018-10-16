%% Proximal point method
% [Spingarn]
function results = SolveSCPPP(this, slices, node_price, link_price, options)
%% Global Declaration
global ITER_LIMIT; 
if isempty(ITER_LIMIT)
	ITER_LIMIT = inf;
end
%%
if ~isfield(options, 'PenaltyFactor')
	% 	options.PenaltyFactor = 0.1*mean([node_price; link_price]);
	options.PenaltyFactor = 32/length(slices)/mean([node_price; link_price]);
end
M = getParallelInfo();

%% Initialization parameters
r = options.PenaltyFactor;
Ns = length(slices);
num_process = Ns;
options.NumberSlices = Ns;
num_links = this.NumberLinks;
num_duals = this.NumberLinks + this.NumberDataCenters;
gamma_k = zeros(num_duals, num_process); % auxiliary variables to lambda.
q_k = zeros(num_duals, num_process); % Dual variables for the dual ADMM formulation.
x_prox = cell(num_process,1);
fval_k = zeros(num_process,1);
fval__ = 0;
output_k = cell(num_process,1);
loads = zeros(num_duals, num_process);
for si = 1:num_process
	sl = slices{si};
	dc_id = sl.getDCPI;
	sl.prices.Link = link_price;
	sl.prices.Link = sl.prices.Link(sl.VirtualLinks.PhysicalLink);
	sl.prices.Node = node_price;
	sl.prices.Node = sl.prices.Node(dc_id);
	x_prox{si} = zeros(sl.num_flow_vars,1);
end
k = 1;
t1 = tic;
while true
	%% Step 1: update the dual-variables
	lambda_k = 1/num_process*sum(gamma_k,2);
	x_prox__ = x_prox;
	
	%% Step-2: solve the sup-problems, return the resource load
	parfor (sj = 1:num_process,M)
% 	for sj = 1:num_process
		sl = slices{sj};
		dc_id = sl.getDCPI;
		idx = [sl.VirtualLinks.PhysicalLink; num_links+dc_id];
		lambda = lambda_k;
		lambda = lambda(idx);
		q = q_k(:,sj);
		q = q(idx);
		[gamma_k(:,sj), fval_k(sj), x_prox{sj}, output_k{sj}] = ...
			sl.priceOptimalFlowRatePP([], lambda, q, x_prox{sj}, options);
	end
	%% Step 3: update the auxiliary variables
	q_k0 = q_k + r*(lambda_k - gamma_k);
	q_k = q_k0 - 1/num_process*sum(q_k0,2);
	%% Step-4: stop condition test
	re_norm = norm(mean(gamma_k,2)-lambda_k);
	xe_norm = 0 ;
	xsum = 0;
	for sj = 1:num_process
		xe_norm = xe_norm + sum((x_prox{sj} - x_prox__{sj}).^2);
		xsum = xsum + sum((x_prox{sj}).^2);
		loads(:,sj) = output_k{sj}.loads;
	end
	xe_norm = sqrt(xe_norm);
	tol_primal = 10^-3 * sqrt(xsum);
	tol_dual = 10^-3 * norm(lambda_k);
	fval = sum(fval_k);
	fval_change = fval - fval__;
	num_iters = 0;
	num_funccount = 0;
	for si = 1:num_process
		num_iters = num_iters + output_k{si}.iterations;
		num_funccount = num_funccount + output_k{si}.funcCount;
	end
	num_iters = round(num_iters/num_process);
	num_funccount = round(num_funccount/num_process);
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
	sl = slices{si};
	sl.saveResults(output_k{si});
	sl.prices.Link = [];
	sl.prices.Node = [];
end
results.LinkPrice = lambda_k(1:num_links) + link_price;
results.NodePrice = lambda_k(num_links+(1:this.NumberDataCenters)) + node_price;
results.lambda = lambda_k;
results.x_prox = x_prox;
results.loads = loads;
results.numiters = k;
%% Capacity Distribution
% (a) Distribute the capacity according to the value of 'q';
% (b) Distribute the capacity accroding to the load;
if isfield(options, 'capacities')
	capacities = options.capacities/Ns;
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
	results.capacities = capacities - delta_capacity;
end

end