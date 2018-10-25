% Optimize resource price
% * Resource Cost Model: linear, convex (quatratic)
% DATE: 2017-04-23
function [output, runtime] = optimizeResourcePriceNew(this, slices, options)
global DEBUG;

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
    sl = slices{i};
    capacities.Node(sl.getDCPI) = capacities.Node(sl.getDCPI) + ...
      sl.ServiceNodes.Capacity;
    capacities.Link(sl.Links.PhysicalLink) = ...
      capacities.Link(sl.Links.PhysicalLink) + sl.Links.Capacity;
  end
end
options.ResidualCapacity = capacities;
options.Slices = slices;
Ndc = this.NumberDataCenters;
Ns = length(slices);
Ne = this.NumberLinks;
for i = 1:Ns
  slices{i}.initialize();
end
if nargin <= 2 || ~isfield(options, 'PricingMethod')
  iter_method = 'RandomizeCost';
else
  iter_method = options.PricingMethod;
  options = rmfield('PricingMethod');
end
switch iter_method
  case 'RandomizeCost'
    st = rng;
    rng(20180909);
    unitcost.Link = this.getLinkCost().*((rand(Ne,1)-0.5)/100+1);
    unitcost.Node = this.getNodeCost().*((rand(Ndc,1)-0.5)/100+1);
    rng(st);
  case 'UniformCost'
    unitcost.Link = ones(Ne,1)*mean(this.getLinkCost());
    unitcost.Node = ones(Ndc,1)*mean(this.getNodeCost());
  case 'OriginCost'
    unitcost.Link = this.getLinkCost();
    unitcost.Node = this.getNodeCost();
  otherwise
    error('error: %s.',iter_method);
end
options = structmerge(options, getstructfields(this.options, {'Threshold','Form'}));
options.PricingPolicy = 'quadratic-price';
options.Stage = 'temp';
options.CapacityConstrained = false;
switch options.Threshold
  case {'min', 'average', 'max'}
    b_profit_ratio = true;
  otherwise
    b_profit_ratio = false;
end
if nargout>=2
  options.CountTime = true;
  t_start = tic;
else
  options.CountTime = false;
end

%% Social-welfare aware price adjustment
% Initialize Price
prices = this.op.initPrice(this, slices, unitcost, options);
t0 = 10^-1;     % {1|0.1|0.01}
delta_price.Link = t0 * unitcost.Link;  % init_price.link
delta_price.Node = t0 * unitcost.Node;

number_iter = 1;
%% Initial Price
% If we specify the initial price and after the first iteration, we find that the profit
% of SP is decreased, we guess that the optimal price is between the |t1* unitcost.Link| and
% |init_prices.Link| for link.
% If we do not specify the initial price, then we start from a relatively small prices
% (unitcost.Link), usually the SP's profit will iNdcrease with the iNdcrease of the price.
% However, if the profit is decreasing, then the initial guess of price (unitcost.Link) is still
% to large, and we set the initial price to zero.
% guess_init_prices.Link = t1*unitcost.Link;
% guess_init_prices.Node = t1*unitcost.Node;
% if find(guess_init_prices.Link>prices.Link,1)
%     guess_init_prices.Link = prices.Link;
% end
% if find(guess_init_prices.Node>prices.Node,1)
%     guess_init_prices.Node = prices.Node;
% end
% price_prev.Link = [guess_init_prices.Link, prices.Link];
% price_prev.Node = [guess_init_prices.Node, prices.Node];
price_prev.Link = [zeros(Ne,1), prices.Link];
price_prev.Node = [zeros(Ndc,1), prices.Node];
if b_profit_ratio
  b_forced_break = false;
end
sp_profit = SolveSCP(slices, prices);

b_initial_trial = true;
while true
  number_iter = number_iter + 1;
  %%%
  % Adjust the step according to how much the capacity constraints have been violated.
  % If the resource is over provisioned, the multiplier is larger than 2.
  % Resources with high utilization ratio are likely to be bottleneck. Therefore
  % the iNdcrease amount of those resources is larger. Thus we let the iNdcrease amount
  % associated with the utilization ratio.
  %
  % If the capacity tends to infinity, |delta_price.Link| and |delta_price.Node| stay
  % the same, while |prices.Link| and |prices.Node| still iNdcreases in a constant rate.
  %
  % we only increase the price of those resources that are utilized (resource
  % utilization θ>0), siNdce iNdcreasing the price of idle resources will not iNdcrease the
  % profit of SP.
  load = this.getNetworkLoad(slices, options);
  delta_price.Link = delta_price.Link.*(1+min(1,load.Link./capacities.Link));
  delta_price.Node = delta_price.Node.*(1+min(1,load.Node./capacities.Node));
  node_id = load.Node>0;
  link_id = load.Link>0;
  prices.Link(link_id) = prices.Link(link_id) + delta_price.Link(link_id);
  prices.Node(node_id) = prices.Node(node_id) + delta_price.Node(node_id);
  
  SolveSCP(slices, prices);
  sp_profit_new = this.getSliceProviderProfit(prices, options);
  %% Stop condtion
  % if the profit of SP is non-iNdcreasing, or the profit ratio reaches the predefined
  % threshold, then no need to further iNdcrease the resource prices.
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
      delta_price.Link = t0 * unitcost.Link;  % init_price.link
      delta_price.Node = t0 * unitcost.Node;
      SolveSCP(slices, prices);
      sp_profit = this.getSliceProviderProfit(prices, options);
      continue;
    else
      price_prev.Link(:,1) = price_prev.Link(:,2);
      price_prev.Link(:,2) = prices.Link;
      price_prev.Node(:,1) = price_prev.Node(:,2);
      price_prev.Node(:,2) = prices.Node;
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
  price_prev.Node(:,2) = prices.Node;
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
      SolveSCP(slices, price_middle(i));
      sp_profit_new(i) = this.getSliceProviderProfit(price_middle(i), options);
    end
    if sp_profit_new(1) > sp_profit_new(2)
      price_prev.Node(:,2) = price_middle(2).Node;
      price_prev.Link(:,2) = price_middle(2).Link;
    else
      price_prev.Node(:,1) = price_middle(1).Node;
      price_prev.Link(:,1) = price_middle(1).Link;
    end
    %%%
    % the stop condition can also be set as the differeNdce of price.
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
% improve the convergeNdce rate. Alternatively, one can reset the two vectors as follows
%
%    delta_price.Link = t0 * unitcost.Link;         % init_price.link
%    delta_price.Node = t0 * unitcost.Node;
k = 1;
while true
  %%% Compute the new resource price according to the resource consumption
  load = this.getNetworkLoad(slices, options);
  b_link_violate = (capacities.Link-load.Link) < 1;
  b_node_violate = (capacities.Node-load.Node) < 1;
  if isempty(find(b_link_violate==1,1)) && isempty(find(b_node_violate==1,1))
    break;
  end
  prices.Link(b_link_violate) = prices.Link(b_link_violate) + delta_price.Link(b_link_violate);
  delta_price.Link(b_link_violate) = delta_price.Link(b_link_violate) .* ...
    (load.Link(b_link_violate)./capacities.Link(b_link_violate));     % {2|(k+1)/k}
  prices.Node(b_node_violate) = prices.Node(b_node_violate) + delta_price.Node(b_node_violate);
  delta_price.Node(b_node_violate) = delta_price.Node(b_node_violate) .* ...
    (load.Node(b_node_violate)./capacities.Node(b_node_violate));     % {2|(k+1)/k}
  % Slices solve P1 with $��_k$, return the node (link) load v(y);
  % annouNdce the resource price and optimize each network slice
  number_iter = number_iter + 1;
  SolveSCP(slices, prices);
  k = k+1;
end

if k>1
  delta_price.Link = t0 * prices.Link;  % 0.01 * init_price.link
  delta_price.Node = t0 * prices.Node;
  min_delta_price.Link = delta_price.Link;
  min_delta_price.Node = delta_price.Node;
  d0 = 10^-1;
  d1 = 10^-0;
  stop_cond1 = ~isempty(find(delta_price.Link > d0 * unitcost.Link, 1));
  stop_cond2 = ~isempty(find(delta_price.Node > d0 * unitcost.Node, 1));
  if b_profit_ratio
    stop_cond3 = this.checkProfitRatio(prices, options);
  else
    sp_profit = this.getSliceProviderProfit(prices, options);
    stop_cond3 = true;
  end
  partial_link_violate = false(Ne, 1);
  partial_node_violate = false(Ndc, 1);
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
    SolveSCP(slices, prices);
    load = this.getNetworkLoad(slices, options);
    
    if b_profit_ratio
      % the profit ratio of SP should not less than the predefined threshold.
      stop_cond3 = this.checkProfitRatio(prices, options);
    else
      % we decrease the price, the profit of SP should iNdcrease.
      sp_profit_new = this.getSliceProviderProfit(prices, options);
      stop_cond3 = sp_profit_new >= sp_profit;
    end
    b_link_violate = (capacities.Link - load.Link)<0;
    b_node_violate = (capacities.Node - load.Node)<0;
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
      partial_link_violate = false(Ne, 1);
      partial_node_violate = false(Ndc, 1);
    else
      b_first = false;
      prices.Link(b_link) = prices.Link(b_link) + delta_price.Link(b_link);
      prices.Node(b_node) = prices.Node(b_node) + delta_price.Node(b_node);
      if ~stop_cond3 && assert_link_1 && assert_node_1
        SolveSCP(slices, prices);
        break;
      end
      %%%
      %  If $\Delta_\rho$ has been smaller than the initial step, then only those
      %  resources with residual capacity will continue reduce their price, i.e. the
      %  components of step $\Delta_\rho$ corresponding to those overloaded
      %  resources is set to 0.
      assert_link_2 = isempty(find(delta_price.Link > d1 * unitcost.Link, 1));		% the vector is less than a threshold
      assert_node_2 = isempty(find(delta_price.Node > d1 * unitcost.Node, 1));		% the vector is less than a threshold
      if assert_link_2
        partial_link_violate = partial_link_violate | b_link_violate;
        delta_price.Link(partial_link_violate) = 0;
      else
        partial_link_violate = false(Ne, 1);
      end
      delta_price.Link = delta_price.Link / 2;
      min_delta_price.Link = min(delta_price.Link/4, min_delta_price.Link);
      if assert_node_2
        partial_node_violate = partial_node_violate | b_node_violate;
        delta_price.Node(partial_node_violate) = 0;
      else
        partial_node_violate = false(Ndc, 1);
      end
      delta_price.Node = delta_price.Node / 2;
      min_delta_price.Node = min(delta_price.Node/4, min_delta_price.Node);
    end
    %     stop_cond1 = norm(delta_price.Link) > norm(10^-4 * unitcost.Link);
    %     stop_cond2 = norm(delta_price.Node) > norm(10^-4 * unitcost.Node);
    stop_cond1 = ~isempty(find(delta_price.Link > d0 * unitcost.Link, 1));
    stop_cond2 = ~isempty(find(delta_price.Node > d0 * unitcost.Node, 1));
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
if nargout == 2
  runtime = toc(t_start);
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

end
