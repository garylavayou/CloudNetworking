%% Allocate resources with fixed pirces
% DATE: 2019-01-22
%		[stat, slice_stat, results] = fixResourcePricing(this, slices, options)
function varargout = fixResourcePricing(this, varargin)
[slices, capacities, ~, fields, options] = ...
	this.preOptimizeResourcePrice(varargin{:});
defaultopts = Dictionary('PricingPolicy', 'linear', 'unit', 1); 
options = structmerge(defaultopts, options, 'exclude');
if options.bCountTime
	t_start = tic; prt = 0; srt = 0;
end
Ns = length(slices);
for i = 1:Ns
  slices(i).initialize();
	slices(i).Optimizer.initializeParallel('priceOptimalFlowRateDA', options);
end
if options.bCountTime
	t_stop = toc(t_start); prt = prt + t_stop; srt = srt + t_stop; t_start = tic;
end

%% allocate resources with fixed prices
trial_price = [this.Links.Price; this.DataCenters.Price];
[sp_profit, b_violate, output] = this.SolveSCPCC(slices, trial_price, options);  %#ok<ASGLU>
RECORD_SUBROUTINE_TIME;
num_iters = 1;

if b_violate
	%% Partition resources for slices
	part_loads = zeros(this.NumberLinks+this.NumberDataCenters, Ns);
	for si = 1:Ns
		sl = slices(si);
		idx = [sl.Links.PhysicalLink; this.NumberLinks+sl.getDCPI()];
		part_loads(idx, si) = [sl.getLinkLoad(false); sl.getNodeLoad(false)];
	end
	ratio = output.loads./capacities;
	part_loads(output.violates,:) = part_loads(output.violates,:)./ratio(output.violates);
	options.CapacityConstrained = true;
	options.Capacities = part_loads;
	options.bInitialize = true;
	[sp_profit, b_violate, output] = this.SolveSCPCC(slices, trial_price, options);  %#ok<ASGLU>
	RECORD_SUBROUTINE_TIME;
	num_iters = num_iters + 1;
end
if options.bCountTime
	t_start = tic;
end
this.finalize(this.convertParameter(trial_price), slices);
if options.bCountTime
	t_stop = toc(t_start); prt = prt + t_stop/Ns; srt = srt + t_stop;
	this.op.runtime = struct('Parallel', prt, 'Serial', srt);
	this.op.iterations = num_iters;
end
[varargout{1:nargout}] = this.calculateOutput(slices, fields, ...
	struct('PricingPolicy', 'linear'));
end