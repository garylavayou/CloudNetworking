% Optimize resource price
% * Resource Cost Model: linear, convex (quatratic)
% DATE: 2017-04-23
function [output, runtime] = optimizeResourcePriceNew(this, init_price, slices)
global DEBUG;

options = getstructfields(this.options, {'Threshold','Form'});
options.PricingPolicy = 'quadratic-price';
switch options.Threshold
	case {'min', 'average', 'max'}
		b_profit_ratio = true;
	otherwise
		b_profit_ratio = false;
end
% this.clearStates;
if nargout == 2
	slice_runtime = 0;
	runtime.Serial = 0;
	runtime.Parallel = 0;
	options.CountTime = true;
else
	options.CountTime = false;
end
if nargin <= 2
	slices = this.slices;       % all slices are involved in slice dimensioning
end
options.Slices = slices;
options.Stage = 'temp';

% network data
Nc = this.NumberDataCenters;
Ns = length(slices);
Nl = this.NumberLinks;
if nargin <= 2
	node_capacity = this.readDataCenter('Capacity');
	link_capacity = this.readLink('Capacity');
else
	%%
	% residual capacity + reallocatable capacity.
	node_capacity = this.readDataCenter('ResidualCapacity');
	link_capacity = this.readLink('ResidualCapacity');
	for i = 1:Ns
		sl = slices{i};
		node_capacity(sl.getDCPI) = node_capacity(sl.getDCPI) + ...
			sl.VirtualDataCenters.Capacity;
		link_capacity(sl.VirtualLinks.PhysicalLink) = ...
			link_capacity(sl.VirtualLinks.PhysicalLink) + sl.VirtualLinks.Capacity;
	end
end
for i = 1:Ns
	slices{i}.initialize();
end
link_uc = this.getLinkCost;
node_uc = this.getNodeCost;

%% Social-welfare aware price adjustment
% Initial Price
t1 = 1;           % {0.1|0.8|1}
if nargin >=2 && ~isempty(init_price)
	link_usage = false(this.NumberLinks,1);
	node_usage = false(this.NumberDataCenters,1);
	for i = 1:Ns
		link_usage(slices{i}.VirtualLinks.PhysicalLink) = true;
		node_usage(slices{i}.getDCPI) = true;
	end
	if find([init_price.Link(link_usage)==0;init_price.Node(node_usage)==0],1)
		prices.Link = t1*link_uc;
		prices.Node = t1*node_uc;
		%         init_price.Link = prices.Link
		%         init_price.Node = prices.Node;
	else
		prices.Link = init_price.Link;
		prices.Node = init_price.Node;
	end
else
	prices.Link = t1* link_uc;
	price.Node = t1* node_uc;
end
for i = 1:Ns
	slices{i}.VirtualLinks{:,'Price'} = 0;      % Reset it to make <isFinal()=false>.
	slices{i}.VirtualDataCenters{:,'Price'} = 0;
end
t0 = 10^-1;     % {1|0.1|0.01}
delta_price.Link = t0 * link_uc;  % init_price.link
delta_price.Node = t0 * node_uc;

number_iter = 1;
%% Initial Price
% If we specify the initial price and after the first iteration, we find that the profit
% of SP is decreased, we guess that the optimal price is between the |t1* link_uc| and
% |init_link_price| for link.
% If we do not specify the initial price, then we start from a relatively small prices
% (link_uc), usually the SP's profit will increase with the increase of the price.
% However, if the profit is decreasing, then the initial guess of price (link_uc) is still
% to large, and we set the initial price to zero.
% guess_init_link_price = t1*link_uc;
% guess_init_node_price = t1*node_uc;
% if find(guess_init_link_price>prices.Link,1)
%     guess_init_link_price = prices.Link;
% end
% if find(guess_init_node_price>prices.Node,1)
%     guess_init_node_price = prices.Node;
% end
% price_prev.Link = [guess_init_link_price, prices.Link];
% price_prev.Node = [guess_init_node_price, prices.Node];
price_prev.Link = [zeros(this.NumberLinks,1), prices.Link];
price_prev.Node = [zeros(this.NumberDataCenters,1), prices.Node];
if b_profit_ratio
	b_forced_break = false;
end
SolveSCP(prices);
sp_profit = this.getSliceProviderProfit(prices, options);

b_initial_trial = true;
while true
	number_iter = number_iter + 1;
	%%%
	% Adjust the step according to how much the capacity constraints have been violated.
	% If the resource is over provisioned, the multiplier is larger than 2.
	% Resources with high utilization ratio are likely to be bottleneck. Therefore
	% the increase amount of those resources is larger. Thus we let the increase amount
	% associated with the utilization ratio.
	%
	% If the capacity tends to infinity, |delta_price.Link| and |delta_price.Node| stay
	% the same, while |prices.Link| and |prices.Node| still increases in a constant rate.
	%
	% we only increase the price of those resources that are utilized (resource
	% utilization θ>0), since increasing the price of idle resources will not increase the
	% profit of SP.
	load = this.getNetworkLoad(slices, options);
	delta_price.Link = delta_price.Link.*(1+min(1,load.Link./link_capacity));
	delta_price.Node = delta_price.Node.*(1+min(1,load.Node./node_capacity));
	node_id = load.Node>0;
	link_id = load.Link>0;
	prices.Link(link_id) = prices.Link(link_id) + delta_price.Link(link_id);
	price.Node(node_id) = prices.Node(node_id) + delta_price.Node(node_id);
	
	SolveSCP(prices);
	sp_profit_new = this.getSliceProviderProfit(prices, options);
	%% Stop condtion
	% if the profit of SP is non-increasing, or the profit ratio reaches the predefined
	% threshold, then no need to further increase the resource prices.
	if sp_profit >= sp_profit_new
		%         if ~isempty(find(price_prev.Link(:,1),1)) || ...
		%                 ~isempty(find(price_prev.Node(:,1),1))
		if b_initial_trial
			% initial price is too high
			price_prev.Link(:,2) = price_prev.Link(:,2)/2;
			price_prev.Node(:,2) = price_prev.Node(:,2)/2;
			prices.Link = price_prev.Link(:,2);
			prices.Node = price_prev.Node(:,2);
			t0 = t0/2;
			delta_price.Link = t0 * link_uc;  % init_price.link
			delta_price.Node = t0 * node_uc;
			SolveSCP(prices);
			sp_profit = this.getSliceProviderProfit(prices, options);
			continue;
		else
			price_prev.Link(:,1) = price_prev.Link(:,2);
			price_prev.Link(:,2) = price.Link;
			price_prev.Node(:,1) = price_prev.Node(:,2);
			price_prev.Node(:,2) = prices;
			break;
		end
	end
	if b_initial_trial
		b_initial_trial = false;
	end
	sp_profit = sp_profit_new;
	price_prev.Link(:,1) = price_prev.Link(:,2);
	price_prev.Link(:,2) = prices.Link;
	price_prev.Node(:,1) = price_prev.Node(:,2);
	price_prev.Node(:,2) = prices.node;
	if b_profit_ratio && this.checkProfitRatio(prices, options)
		b_forced_break = true;
		break;
	end
end
%%
% If the last step is not stopped by the profit ratio, we need to further search the
% optimal price.
if ~b_profit_ratio || ~b_forced_break
	sp_profit_new = [1 0];
	epsilon = 10^-3;
	while true  %|| (h-l) > 0.05
		number_iter = number_iter + 1;
		price_middle = struct('Node', ...
			{(2/3)*price_prev.Node(:,1)+(1/3)*price_prev.Node(:,2), (1/3)*price_prev.Node(:,1)+(2/3)*price_prev.Node(:,2)},...
			'Link', ...
			{(2/3)*price_prev.Link(:,1)+(1/3)*price_prev.Link(:,2), (1/3)*price_prev.Link(:,1)+(2/3)*price_prev.Link(:,2)});
		for i = 1:2
			SolveSCP(price_middle(i));
			sp_profit_new(i) = this.getSliceProviderProfit(price_middle, options);
		end
		if sp_profit_new(1) > sp_profit_new(2)
			price_prev.Node(:,2) = price_middle(2).Node;
			price_prev.Link(:,2) = price_middle(2).Link;
		else
			price_prev.Node(:,1) = price_middle(1).Node;
			price_prev.Link(:,1) = price_middle(1).Link;
		end
		%%%
		% the stop condition can also be set as the difference of price.
		if abs((max(sp_profit_new)-sp_profit)/sp_profit) < epsilon
			break;
		else
			sp_profit = max(sp_profit_new);
		end
	end
	prices.Node = price_prev.Node(:,1);       % temp_node_price
	prices.Link = price_prev.Link(:,1);       % temp_link_price
end
%%
% |delta_price.Link| and |delta_price.Node| of the first step can still be used, to
% improve the convergence rate. Alternatively, one can reset the two vectors as follows
%
%    delta_price.Link = t0 * link_uc;         % init_price.link
%    delta_price.Node = t0 * node_uc;
k = 1;
while true
	%%% Compute the new resource price according to the resource consumption
	load = this.getNetworkLoad(slices, options);
	b_link_violate = (link_capacity-load.Link) < 1;
	b_node_violate = (node_capacity-load.Node) < 1;
	if isempty(find(b_link_violate==1,1)) && isempty(find(b_node_violate==1,1))
		break;
	end
	prices.Link(b_link_violate) = prices.Link(b_link_violate) + delta_price.Link(b_link_violate);
	delta_price.Link(b_link_violate) = delta_price.Link(b_link_violate) .* ...
		(load.Link(b_link_violate)./link_capacity(b_link_violate));     % {2|(k+1)/k}
	prices.Node(b_node_violate) = prices.Node(b_node_violate) + delta_price.Node(b_node_violate);
	delta_price.Node(b_node_violate) = delta_price.Node(b_node_violate) .* ...
		(load.Node(b_node_violate)./node_capacity(b_node_violate));     % {2|(k+1)/k}
	% Slices solve P1 with $��_k$, return the node (link) load v(y);
	% announce the resource price and optimize each network slice
	number_iter = number_iter + 1;
	SolveSCP(prices);
	k = k+1;
end

if k>1
	delta_price.Link = t0 * prices.Link;  % 0.01 * init_price.link
	delta_price.Node = t0 * prices.Node;
	min_delta_price.Link = delta_price.Link;
	min_delta_price.Node = delta_price.Node;
	d0 = 10^-1;
	d1 = 10^-0;
	stop_cond1 = ~isempty(find(delta_price.Link > d0 * link_uc, 1));
	stop_cond2 = ~isempty(find(delta_price.Node > d0 * node_uc, 1));
	if b_profit_ratio
		stop_cond3 = this.checkProfitRatio(prices, options);
	else
		sp_profit = this.getSliceProviderProfit(prices, options);
		stop_cond3 = true;
	end
	partial_link_violate = false(Nl, 1);
	partial_node_violate = false(Nc, 1);
	b_first = true;
	while (stop_cond1 || stop_cond2) && stop_cond3
		number_iter = number_iter + 1;
		if ~isempty(DEBUG) && DEBUG
			disp('----link price    delta link price----')
			disp([prices.Link delta_price.Link]);
		end
		b_link = prices.Link > delta_price.Link;
		prices.Link(b_link) = prices.Link(b_link) - delta_price.Link(b_link);
		if ~isempty(DEBUG) && DEBUG
			disp('----node price    delta node price----')
			disp([prices.Node delta_price.Node]);
		end
		b_node = prices.Node > delta_price.Node;
		prices.Node(b_node) = prices.Node(b_node) - delta_price.Node(b_node);
		SolveSCP(prices);
		load = this.getNetworkLoad(slices, options);
		
		if b_profit_ratio
			% the profit ratio of SP should not less than the predefined threshold.
			stop_cond3 = this.checkProfitRatio(prices, options);
		else
			% we decrease the price, the profit of SP should increase.
			sp_profit_new = this.getSliceProviderProfit(prices, options);
			stop_cond3 = sp_profit_new >= sp_profit;
		end
		b_link_violate = (link_capacity - load.Link)<0;
		b_node_violate = (node_capacity - load.Node)<0;
		assert_link_1 = isempty(find(b_link_violate==1,1));			% no violate link
		assert_node_1 = isempty(find(b_node_violate==1,1));			% no violate node
		if assert_link_1 && assert_node_1 && stop_cond3
			if b_first
				delta_price.Link = delta_price.Link * 2;
				delta_price.Node = delta_price.Node * 2;
			else
				delta_price.Link = delta_price.Link + min_delta_price.Link;
				delta_price.Node = delta_price.Node + min_delta_price.Node;
			end
			partial_link_violate = false(Nl, 1);
			partial_node_violate = false(Nc, 1);
		else
			b_first = false;
			prices.Link(b_link) = prices.Link(b_link) + delta_price.Link(b_link);
			prices.Node(b_node) = prices.Node(b_node) + delta_price.Node(b_node);
			if ~stop_cond3 && assert_link_1 && assert_node_1
				SolveSCP(prices);
				break;
			end
			%%%
			%  If $\Delta_\rho$ has been smaller than the initial step, then only those
			%  resources with residual capacity will continue reduce their price, i.e. the
			%  components of step $\Delta_\rho$ corresponding to those overloaded
			%  resources is set to 0.
			assert_link_2 = isempty(find(delta_price.Link > d1 * link_uc, 1));		% the vector is less than a threshold
			assert_node_2 = isempty(find(delta_price.Node > d1 * node_uc, 1));		% the vector is less than a threshold
			if assert_link_2
				partial_link_violate = partial_link_violate | b_link_violate;
				delta_price.Link(partial_link_violate) = 0;
			else
				partial_link_violate = false(Nl, 1);
			end
			delta_price.Link = delta_price.Link / 2;
			min_delta_price.Link = min(delta_price.Link/4, min_delta_price.Link);
			if assert_node_2
				partial_node_violate = partial_node_violate | b_node_violate;
				delta_price.Node(partial_node_violate) = 0;
			else
				partial_node_violate = false(Nc, 1);
			end
			delta_price.Node = delta_price.Node / 2;
			min_delta_price.Node = min(delta_price.Node/4, min_delta_price.Node);
		end
		%     stop_cond1 = norm(delta_price.Link) > norm(10^-4 * link_uc);
		%     stop_cond2 = norm(delta_price.Node) > norm(10^-4 * node_uc);
		stop_cond1 = ~isempty(find(delta_price.Link > d0 * link_uc, 1));
		stop_cond2 = ~isempty(find(delta_price.Node > d0 * node_uc, 1));
	end
end

%% Finalize substrate network
% # The resource allocation variables, virtual node/link load, and flow rate of each
% slice.
% # After the optimization, each network slice has record the final prices.
% # Record the substrate network's node/link load, price.
this.finalize(prices, slices);

% Calculate the output
if nargout >= 1
	output = this.calculateOutput([], getstructfields(options, {'Slices','PricingPolicy'}));
end

% output the optimization results
if DEBUG
	fprintf('Optimization results:\n');
	fprintf('\tThe optimization procedure contains %d iterations.\n', number_iter);
	if nargout >= 1
		fprintf('\tOptimal objective value: %d.\n', output.Welfare);
	end
	fprintf('\tNormalized network utilization is %G.\n', this.utilizationRatio);
end

%% sub-problem
	function SolveSCP(node_price_t, link_price_t)
		for s = 1:Ns
			sl = slices{s};
			sl.prices.Link = link_price_t(sl.VirtualLinks.PhysicalLink);
			% |prices.Node| only contain the price of data center nodes.
			dc_id = sl.getDCPI;
			sl.prices.Node = node_price_t(dc_id);
			options.ResidualCapacity.Link = link_capacity(sl.VirtualLinks.PhysicalLink);
			options.ResidualCapacity.Node = node_capacity(dc_id);
			if options.CountTime
				tic;
			end
			%%%
			% dynamically decide which version of _priceOptimalFlowRate_ to call.
			% See also <Slice.priceOptimalFlowRate> and
			% <DynamciSlice.priceOptimalFlowRate>.
			sl.priceOptimalFlowRate([], options);
			if options.CountTime
				t = toc;
				slice_runtime = max(slice_runtime, t);
				runtime.Serial = runtime.Serial + t;
			end
			sl.prices.Link = [];
			sl.prices.Node = [];
		end
		if options.CountTime
			runtime.Parallel = runtime.Parallel + slice_runtime;
		end
	end
end
