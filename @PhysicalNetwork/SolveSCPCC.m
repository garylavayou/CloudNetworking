function [sp_profit, b_violate, violates] = SolveSCPCC(this, slices, prices, options)
Ns = length(slices);
Ne = this.NumberLinks;
Ndc = this.NumberDataCenters;
num_process = Ns;
output_k = cell(num_process,1);
nout = nargout;
if nout >= 2
	loads = zeros(Ne+Ndc, num_process);
end
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
if M > 0
	options.bParallel = true;
	parfor (sj = 1:Ns,M)    % test if op.options contains handle to other objects.
		output_k{sj} = par_ops(sj).priceOptimalFlowRate([], options);
		% slices(sj).Optimizer = par_op;  % no data send back to the optimizer.
	end
else
	for sj = 1:num_process
		output_k{sj} = par_ops(sj).priceOptimalFlowRate([], options);
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
	if nargout >= 2
		idx = [sl.Links.PhysicalLink; this.NumberLinks+sl.getDCPI()];
		loads(idx, si) = [sl.getLinkLoad(false); sl.getNodeLoad(false)];
	end
end

%% Output
sp_profit = this.getSliceProviderProfit(slices, prices, ...
  getstructfields(options, {'Slices','PricingPolicy'}, 'error'));
if nout >= 2
  if ~options.CapacityConstrained
    violates = sum(loads,2) > [options.ResidualCapacity.Link; options.ResidualCapacity.Node];
    b_violate = ~isempty(find(violates,1));
  else
    violates = [];
    b_violate = false;
  end
end
end
