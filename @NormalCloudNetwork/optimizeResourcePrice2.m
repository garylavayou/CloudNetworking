% Optimize resource price
% Firstly, search a price 'αρ' that close to the optimal price.
% Then use "Dual-ADMM" or "Partial Inverse" method to find the final solution.
% * Resource Cost Model: linear
% DATE: 2018-12-20
function [output, runtime] = optimizeResourcePrice2(this, slices, options)
global DEBUG ITER_LIMIT; %#ok<NUSED>
if isempty(ITER_LIMIT)
	ITER_LIMIT = inf;
end

%% Initialization
if nargin <= 1 || isempty(slices)
	% all slices are involved in slice dimensioning
	slices = this.slices;
  capacities.Node = this.readDataCenter('Capacity');
  capacities.Link = this.readLink('Capacity');
else
	% residual capacity + reallocatable capacity.
	capacities.Node = this.readDataCenter('ResidualCapacity');
  capacities.Link = this.readLink('ResidualCapacity');
	for i = 1:length(slices)
    sl = slices(i);
    capacities.Node(sl.getDCPI) = capacities.Node(sl.getDCPI) + ...
      sl.ServiceNodes.Capacity;
    capacities.Link(sl.Links.PhysicalLink) = ...
      capacities.Link(sl.Links.PhysicalLink) + sl.Links.Capacity;
  end	
end

Ndc = this.NumberDataCenters;
Ns = length(slices);
Ne = this.NumberLinks;
if nargin <= 2 || ~isfield(options, 'PricingMethod')
	cost_init_method = 'RandomizeCost';
else
	cost_init_method = options.PricingMethod;
	refield(options, 'PricingMethod');
end
switch cost_init_method
	case 'RandomizeCost'
		st = rng;
		rng(20180909);
		link_uc = this.getLinkCost().*((rand(Ne,1)-0.5)/100+1);
		node_uc = this.getNodeCost().*((rand(Ndc,1)-0.5)/100+1);
		rng(st);
	case 'UniformCost'
		link_uc = ones(Ne,1)*mean(this.getLinkCost());
		node_uc = ones(Ndc,1)*mean(this.getNodeCost());
	case 'OriginCost'
		link_uc = this.getLinkCost();
		node_uc = this.getNodeCost();
	otherwise
		error('error: %s.',cost_init_method);
end
if nargin <= 2 || ~isfield(options, 'TuningMethod')
	iter_method = 'DualADMM';		% {DualDecomposition, DualADMM, PartialInverse}
else
	iter_method = options.TuningMethod;
	refield(options, 'TuningMethod');
end
defaultopts = Dictionary('ResidualCapacity', capacities, ...
	'Slices', slices, ...			% slices is an HetArray
	'PricingPolicy', 'linear', ...
	'Stage', 'temp', ...
	'bFinal', false, ...
	'OptimizeOrder', 0, ...
	'CapacityConstrained', false ...
);
if nargin <= 2
	options = defaultopts;
else
	options = structmerge(defaultopts, options);
end
if nargout == 2 
	options.CountTime = true;
	t_start = tic;
else
	options.CountTime = false;
end
options.bInitialize = true;

%% Parallel initialization
parfor i = 1:Ns
  slices(i).initialize();
	slices(i).Optimizer.initializeParallel('priceOptimalFlowRate', options);
end

%% Trail Prices
% Find the price 'αρ' that possibly maximize the profit of SP.
sp_profit = -inf*ones(3,1);
b_violate = true(3,1);
k = 0;
t = 0;
while true
	trial_price.Node = node_uc * 2^t;
	trial_price.Link = link_uc * 2^t;
	b_violate(1:2) = b_violate(2:3);
	sp_profit(1:2) = sp_profit(2:3);
	[sp_profit(3), b_violate(3)] = this.SolveSCPCC(slices, trial_price, options);  % with price model linear.
	k = k + 1;
	if k == 1
		options.bInitialize = false;
	end
	if sp_profit(3) <= sp_profit(2)
		break;
	end
	t = t+1;
end

beta_L = 2^(t-2);
beta_R = 2^t;
epsilon = 10^-3;
sp_profit_N = zeros(2,1);
b_violate_N = zeros(2,1);
while abs(beta_L-beta_R)/beta_L > epsilon
	beta_N = [beta_L*2+beta_R; beta_L+beta_R*2]/3;
	for i = 1:2
		trial_price.Node = node_uc * beta_N(i);
		trial_price.Link = link_uc * beta_N(i);
		[sp_profit_N(i), b_violate_N(i)]= this.SolveSCPCC(slices, trial_price, options);
		k = k + 2;
	end
	if sp_profit_N(1)<=sp_profit_N(2)
		beta_L = beta_N(1);
		sp_profit(1) = sp_profit_N(1);
		sp_profit(2) = sp_profit_N(2);
	else
		beta_R = beta_N(2);
		sp_profit(3) = sp_profit_N(2);
		sp_profit(2) = sp_profit_N(1);
	end
end
if ~b_violate_N(1)
	output.beta = beta_N(1);
	output.sp_profit = sp_profit_N(1);
	output.LinkPrice = link_uc * output.beta;
	output.NodePrice = node_uc * output.beta;
	output.procedure = 'max';	
else
	%%	Serach the final price
	output.beta = (beta_L+beta_R)/2;
	trial_price.Node = node_uc * output.beta;
	trial_price.Link = link_uc * output.beta;
	for j = 1:Ns
		slices(j).setProblem([], [], [], link_capacity/Ns, node_capacity/Ns);
	end
	switch iter_method
		case 'PartialInverse'
			tmp_output = this.SolveSCPPP(slices, trial_price, options);
		case 'DualDecomposition'
			tmp_output = this.SolveSCPDD(slices, trial_price, options);
		case 'DualADMM'
			tmp_output = this.SolveSCP0(slices, trial_price, options);			% return results
		otherwise
			error('error: un-recognized method.');
	end
	k = k + tmp_output.numiters;
	options.CapacityConstrained = true;
  options.Capacities = tmp_output.Capacities;
  output.Prices = struct('Link', tmp_output.LinkPrice, 'Node',  tmp_output.NodePrice);
	[output.sp_profit, b] = this.SolveSCPCC(slices, output.Prices, options);
	k = k + 1;
	if b
		warning('%s: Capacity violation.', calledby);
	end
	output.procedure = 'dual-iter';
end
output.numiters = k;

%% Finalize substrate network
this.finalize(output.Prices, slices);
output = this.calculateOutput(output, struct('Slices', {slices}, 'PricingPolicy', 'linear'));
output.utilization = this.utilizationRatio;
if nargout == 2
	runtime = toc(t_start);
end

end

