function results = SolveSCPDD(this, slices, prices, options)		% Dual Decomposition
M = getParallelInfo();
%% Initialization parameters
step_length = 0.001;
num_process = length(slices);
num_links = this.NumberLinks;
num_duals = this.NumberLinks + this.NumberDataCenters;
lambda_k = zeros(num_duals, 1);
lambda_k__ = lambda_k;
loads = zeros(num_duals, num_process);
capacities = [link_capacity; node_capacity];
fval_k = zeros(num_process,1);
fval__ = 0;
output_k = cell(num_process,1);
for si = 1:num_process
	sl = slices(si);
	dc_id = sl.getDCPI;
	sl.prices.Link = link_price;
	sl.prices.Link = sl.prices.Link(sl.VirtualLinks.PhysicalLink);
	sl.prices.Node = node_price;
	sl.prices.Node = sl.prices.Node(dc_id);
end
k = 1;
t1 = tic;
while true
	%% Step-1: solve the sup-problems, return the resource load
	parfor (sj = 1:num_process,M)
		%for sj = 1:num_process
		sl = slices(sj);
		dc_id = sl.getDCPI;
		idx = [sl.VirtualLinks.PhysicalLink; num_links+dc_id];
		lambda = lambda_k;
		lambda = lambda(idx);
		[loads(:,sj), fval_k(sj), output_k{sj}] = ...
			sl.priceOptimalFlowRateDD([], lambda, options);
	end
	%% Step-2: according to the load, update the dual varialbes
	lambda_k = max(0,lambda_k + (step_length/k)*(sum(loads,2)-capacities));
	%% Step-3: stop condition test
	re_norm = norm(lambda_k__ - lambda_k);
	if re_norm > 0
		re_norm = re_norm/norm(lambda_k);
	end
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
		fprintf('\nIteration Step-length  Dual-change   Optimality  Iterations Evaluations\n');
		cprintf('*text', ...
			'————————— ——————————— ————————————— ———————————— —————————— ———————————\n');
	end
	fprintf('%8d  %10.4g   %11.4f   %10.4g  %9d   %9d \n',...
		k, step_length/k, fval_change, re_norm, num_iters, num_funccount);
	if mod(k,20) == 0
		fprintf('\n');
	end
	
	if re_norm <= 10^-3
		break;
	else
		lambda_k__ = lambda_k;
		fval__ = fval;
	end
	k = k + 1;
end
t2 = toc(t1);
fprintf('Dual-Decomposition: elapsed time: %d\n', t2);

for si = 1:Ns
	sl = slices(si);
	sl.op.saveTempResults(output_k{si});
end
results.LinkPrice = lambda_k(1:this.NumberLinks) + link_price;
results.NodePrice = lambda_k(this.NumberLinks+(1:this.NumberDataCenters)) + node_price;
end