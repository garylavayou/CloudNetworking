function [output, runtime] = optimizeResourcePrice2(this, slices, options)
global DEBUG ITER_LIMIT; %#ok<NUSED>
if isempty(ITER_LIMIT)
	ITER_LIMIT = inf;
end

%% Initialization
if nargin <= 1 || isempty(slices)
	slices = this.slices;       % all slices are involved in slice dimensioning
	node_capacity = this.getDataCenterField('Capacity');
	link_capacity = this.getLinkField('Capacity');
else
	% residual capacity + reallocatable capacity.
	node_capacity = this.getDataCenterField('ResidualCapacity');
	link_capacity = this.getLinkField('ResidualCapacity');
	for i = 1:length(slices)
		sl = slices{i};
		node_capacity(sl.getDCPI) = node_capacity(sl.getDCPI) + ...
			sl.VirtualDataCenters.Capacity;
		link_capacity(sl.VirtualLinks.PhysicalLink) = ...
			link_capacity(sl.VirtualLinks.PhysicalLink) + sl.VirtualLinks.Capacity;
	end	
end
options.capacities = [link_capacity; node_capacity];
Ns = length(slices);
Ne = this.NumberLinks;
Ndc = this.NumberDataCenters;
if nargin <= 2 || ~isfield(options, 'PricingMethod')
	iter_method = 'RandomizeCost';
else
	iter_method = options.PricingMethod;
end
switch iter_method
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
		error('error: %s.',iter_method);
end
if nargin <= 2 || ~isfield(options, 'TuningMethod')
	iter_method = 'DualADMM';		% {DualDecomposition, DualADMM, PartialInverse}
else
	iter_method = options.TuningMethod;
end
if strcmpi(this.options.Form, 'compact')
	options.bCompact = true;
end
if nargout == 2 
	options.CountTime = true;
	t_start = tic;
else
	options.CountTime = false;
end
options.OptimizeOrder = 0;
options.PricingPolicy = 'linear';
options.CapacityConstrained = false;
sp_profit = -inf*ones(3,1);
t = 0;

%% Parallel initialization
M = getParallelInfo();
array = cell(Ns, 1);
problem = cell(Ns, 1);
indices = cell(Ns, 1);
parfor (pj = 1:Ns,M)
	sl = slices{pj};
	[array{pj}, problem{pj}, indices{pj}] = sl.initializeProblem(options);
end
for j = 1:Ns
	slices{j}.setProblem(array{j}, problem{j}, indices{j});
end

%% Trail Prices
% Find the price 'αp' that possibly maximize the profit of SP.
b_violate = true(3,1);
k = 0;
while true
	trial_price_node = node_uc * 2^t;
	trial_price_link = link_uc * 2^t;
	b_violate(1:2) = b_violate(2:3);
	sp_profit(1:2) = sp_profit(2:3);
	[sp_profit(3), b_violate(3)] = this.SolveSCPCC(slices, trial_price_node, trial_price_link, options);
	k = k + 1;
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
		trial_price_node = node_uc * beta_N(i);
		trial_price_link = link_uc * beta_N(i);
		[sp_profit_N(i), b_violate_N(i)]= this.SolveSCPCC(slices, trial_price_node, trial_price_link, options);
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
	trial_price_node = node_uc * output.beta;
	trial_price_link = link_uc * output.beta;
	for j = 1:Ns
		slices{j}.setProblem([], [], [], link_capacity/Ns, node_capacity/Ns);
	end
	switch iter_method
		case 'PartialInverse'
			tmp_output = this.SolveSCPPP(slices, trial_price_node, trial_price_link, options);
		case 'DualDecomposition'
			tmp_output = this.SolveSCPDD(slices, trial_price_node, trial_price_link, options);
		case 'DualADMM'
			tmp_output = this.SolveSCP(slices, trial_price_node, trial_price_link, options);			% return results
		otherwise
			error('error: un-recognized method.');
	end
	k = k + tmp_output.numiters;
	options.CapacityConstrained = true;
	options.capacities = tmp_output.capacities;
	[output.sp_profit, b] = this.SolveSCPCC(slices, tmp_output.NodePrice, ...
		tmp_output.LinkPrice, options);
	k = k + 1;
	if b
		warning('%s: Capacity violation.', calledby);
	end
	output.LinkPrice = tmp_output.LinkPrice;
	output.NodePrice = tmp_output.NodePrice;
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

