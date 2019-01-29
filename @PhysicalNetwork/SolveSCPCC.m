function [sp_profit, b_violate, output] = SolveSCPCC(this, slices, prices, options)
if options.bCountTime
	prt = 0; srt = 0; t0 = tic;
end
prices = this.convertParameter(prices, 'struct');

Ns = length(slices);
Ne = this.NumberLinks;
Ndc = this.NumberDataCenters;
M = getParallelInfo();
par_ops = SliceOptimizer.empty(Ns,0);
for si = 1:Ns
  sl = slices(si);
	sl.Optimizer.setProblem('Price', prices);
	sl.Optimizer.update_options(options);
	if options.CapacityConstrained
		caps = struct('Link', options.Capacities(1:Ne,si), ...
			'Node', options.Capacities(Ne+(1:Ndc),si));
		sl.Optimizer.setProblem('Capacity', caps);
	end
	if M > 0
		sl.Optimizer.Host = Slice.empty();
	end
	par_ops(si) = sl.Optimizer;
end

output_k = cell(Ns,1);
if options.bCountTime
	t1 = toc(t0);	prt = prt + t1; srt = srt + t1; t0 = tic;
end
if M > 0
	options.bParallel = true;
	% check handle objects in options and remove them.
	fields = fieldnames(options);
	for i = 1:length(fields)
		if ishandle(options.(fields{i}))
			options = rmfield(options, fields{i});
		end
	end
	options.Display = 'off';
	options.Warning = 'off';
	parfor (sj = 1:Ns,M)    % test if op.options contains handle to other objects.
		output_k{sj} = par_ops(sj).priceOptimalFlowRate([], options);
		% slices(sj).Optimizer = par_op;  % no data send back to the optimizer.
	end
	if options.bCountTime
		% not procise for serial running time, due to the heterogenuousity of each slice.
		t1 = toc(t0); prt = prt + t1*M/Ns; srt = srt + t1*M; t0 = tic;
	end
else
	for sj = 1:Ns
		output_k{sj} = par_ops(sj).priceOptimalFlowRate([], options);
	end
	if options.bCountTime
		% not procise for parallel running time, due to the heterogenuousity of each slice.
		t1 = toc(t0);	prt = prt + t1/Ns; srt = srt + t1; t0 = tic;
	end
end
for si = 1:Ns
  sl = slices(si);
	sl.Optimizer.setProblem('Price', []);
	if M > 0
		sl.Optimizer.Host = sl; % re-connect slice with optimizer
		sl.Optimizer.saveTempResults(output_k{si});
		if isfield(options, 'bInitialize') && options.bInitialize
			sl.Optimizer.setProblem('Problem', output_k{si}.problem);
		end
	end
	if options.CapacityConstrained
		sl.Optimizer.setProblem('Capacity', []);
	end
end

%% Output
sp_profit = this.getSliceProviderProfit(slices, prices, ...
  getstructfields(options, {'PricingPolicy'}, 'error'));
if ~options.CapacityConstrained
	output.loads = this.convertParameter(this.getNetworkLoad(slices, options));
	output.loads(output.loads<this.op.options.NonzeroTolerance) = 0;
	output.violates = output.loads > [options.ResidualCapacity.Link; ...
		options.ResidualCapacity.Node];
	b_violate = ~isempty(find(output.violates,1));
else
	b_violate = false;
end
if options.bCountTime
	t1 = toc(t0); prt = prt + t1; srt = srt + t1;
	output.runtime = struct('Parallel', prt, 'Serial', srt);
end
end
