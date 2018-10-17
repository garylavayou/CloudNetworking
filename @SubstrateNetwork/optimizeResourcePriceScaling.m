function [output, runtime] = optimizeResourcePriceScaling(this, slices, options)
global DEBUG; %#ok<NUSED>

%% Initialization
if nargin <= 1 || isempty(slices)
	slices = this.slices;       % all slices are involved in slice dimensioning
	node_capacity = this.readDataCenter('Capacity');
	link_capacity = this.readLink('Capacity');
else
	% residual capacity + reallocatable capacity.
	node_capacity = this.readDataCenter('ResidualCapacity');
	link_capacity = this.readLink('ResidualCapacity');
	for i = 1:length(slices)
		sl = slices{i};
		node_capacity(sl.getDCPI) = node_capacity(sl.getDCPI) + ...
			sl.VirtualDataCenters.Capacity;
		link_capacity(sl.VirtualLinks.PhysicalLink) = ...
			link_capacity(sl.VirtualLinks.PhysicalLink) + sl.VirtualLinks.Capacity;
	end	
end
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
if strcmpi(this.options.Form, 'compact')
	options.bCompact = true;
end
if nargout >= 2 
	options.CountTime = true;
	t_start = tic;
else
	options.CountTime = false;
end
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

%% Trial Prices
% Find the price 'αp' that possibly maximize the profit of SP.
b_violate = true;
k = 0;
while true
	trial_price_node = node_uc * 2^t;
	trial_price_link = link_uc * 2^t;
	b_violate__ = b_violate;
	sp_profit(1:2) = sp_profit(2:3);
	[sp_profit(3), b_violate] = SolveSCPCC(trial_price_node, trial_price_link);
	k = k + 1;
	if sp_profit(3) <= sp_profit(2) && ~b_violate
		break;
	end
	t = t+1;
end
% statisfy condition: sp_profit(3) <= sp_profit(2) && ~b_violate
beta_R = 2^t;
epsilon = 10^-3;
sp_profit_N = zeros(2,1);
if ~b_violate__
	beta_L = 2^(t-2);
	while abs(beta_L-beta_R)/beta_L > epsilon
		beta_N = [beta_L*2+beta_R; beta_L+beta_R*2]/3;
		for i = 2:-1:1
			trial_price_node = node_uc * beta_N(i);
			trial_price_link = link_uc * beta_N(i);
			sp_profit_N(i) = SolveSCPCC(trial_price_node, trial_price_link);
			k = k + 1;
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
	output.sp_profit = sp_profit_N(1);
	output.beta = beta_N(1);
	output.procedure = 'max';
else
	beta_L = 2^(t-1);
	while abs(beta_R-beta_L)/beta_L > epsilon
		beta_N = [beta_L*2+beta_R; beta_L+beta_R*2]/3;
		for i = 1:2
			trial_price_node = node_uc * beta_N(i);
			trial_price_link = link_uc * beta_N(i);
			[sp_profit_N(i), b_violate(i)] = SolveSCPCC(trial_price_node, trial_price_link);
			k = k+ 1;
		end
		if b_violate(2)
			beta_L = beta_N(2); sp_profit(2) = sp_profit_N(2);
		elseif b_violate(1)
			beta_L = beta_N(1); sp_profit(2) = sp_profit_N(1);
			if sp_profit_N(1) > sp_profit_N(2)
				beta_R = beta_N(2); sp_profit(3) = sp_profit_N(2);
			end
		else
			if sp_profit_N(1) > sp_profit_N(2)
				beta_R = beta_N(2); sp_profit(3) = sp_profit_N(2);
			else
				beta_L = beta_N(1); sp_profit(2) = sp_profit_N(1);
			end
		end
	end
	output.sp_profit = sp_profit(3);
	output.beta = beta_R;
	output.procedure = 'feas-max';
end
output.LinkPrice = trial_price_link;
output.NodePrice = trial_price_node;
output.numiters = k;	
%% Finalize substrate network
this.finalize(output.Prices, slices);
output = this.calculateOutput(output, struct('Slices', {slices}, 'PricingPolicy', 'linear'));
output.utilization = this.utilizationRatio;
if nargout == 2
	runtime = toc(t_start);
end


	function [sp_profit, b_violate, violates] = SolveSCPCC(prices, capacities)
		num_process = Ns;
		output_k = cell(num_process,1);
		loads = zeros(Ne+Ndc, num_process);
		nin = nargin;
		for si = 1:num_process
			sl = slices{si};
			dc_id = sl.getDCPI;
			link_id = sl.VirtualLinks.PhysicalLink;
			sl.prices.Link = prices.Link(link_id);
			sl.prices.Node = prices.Node(dc_id);
			if nin >= 3
				sl.setProblem([],[],[],capacities(1:Ne,si), capacities(Ne+(1:Ndc),si));
			end
		end
		nout = nargout;
		opts = options;
		if nargin >= 3 && ~isempty(capacities)
			opts.CapacityConstrained = true;
		end
		parfor (sj = 1:Ns,M)
			%for sj = 1:num_process
			sl = slices{sj};
			if nout >= 2
				[output_k{sj}, loads(:,sj)] = sl.priceOptimalFlowRateCC([], opts);
			else
				output_k{sj} = sl.priceOptimalFlowRateCC([], opts);
			end
		end
		for si = 1:Ns
			sl = slices{si};
			sl.saveResults(output_k{si});
			sl.prices.Link = [];
			sl.prices.Node = [];
			sl.capacities = [];
		end
		
		%% Output
		sp_profit = this.getSliceProviderProfit(prices, struct('PricingPolicy', 'linear'));
		if nargout >= 2
			violates = sum(loads,2)>[link_capacity; node_capacity];
			b_violate = ~isempty(find(violates,1)); 
		end
	end
end

