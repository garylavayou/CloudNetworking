% Optimize resource price
% * Resource Cost Model: linear, convex (quatratic)
% DATE: 2017-04-23
%    [stat, slice_stat, results] = optimizeResourcePrice(this, slices, options)
function varargout = optimizeResourcePrice(this, varargin)
global DEBUG;

[slices, capacities, unitcosts, fields, options] =...
	this.preOptimizeResourcePrice(varargin{:});
% override the network default
defaultopts = Dictionary('bCountTime', false,	'PricingPolicy', 'quadratic-price');  
options = setdefault(options, defaultopts);   
if options.bCountTime
  t_start = tic; prt = 0; srt = 0;
end
Ns = length(slices);
for i = 1:Ns
  slices(i).initialize();
	slices(i).Optimizer.initializeParallel('priceOptimalFlowRate', options);
end
num_resource = this.NumberDataCenters + this.NumberLinks;

%% Social-welfare aware price adjustment
% Initialize Price
trial_prices = this.op.initPrice(slices, unitcosts, options);
t0 = 10^-1;     % {1|0.1|0.01}
delta_price = t0 * unitcosts; % init_price.link
switch this.op.options.Threshold
  case {'min', 'average', 'max'}
    b_profit_ratio = true;
  otherwise
    b_profit_ratio = false;
end
if options.bCountTime
	t_stop = toc(t_start); prt = prt + t_stop; srt = srt + t_stop; t_start = tic;
end

%% Initial Price
% If we specify the initial price and after the first iteration, we find that the profit
% of SP is decreased, we guess that the optimal price is between the |t1* unitcost.Link| and
% |init_prices.Link| for link.
% If we do not specify the initial price, then we start from a relatively small prices
% (unitcost.Link), usually the SP's profit will increase with the increase of the price.
% However, if the profit is decreasing, then the initial guess of price (unitcost.Link) is still
% to large, and we set the initial price to zero.
% guess_init_prices = t1*unitcost;
% if find(guess_init_prices>prices,1)
%     guess_init_prices = prices;
% end
% price_prev = [guess_init_prices, prices];
trial_prices = [zeros(num_resource,2) trial_prices];
if b_profit_ratio
  b_forced_break = false;
end
sp_profit = -inf*ones(3,1); 
[sp_profit(3), b_violate, output] = this.SolveSCPCC(slices, trial_prices(:,end), options);
RECORD_SUBROUTINE_TIME;
options.bInitialize = false;
number_iters = 0;

b_initial_trial = true;
while true
  sp_profit(1:2) = sp_profit(2:3);  trial_prices(:,1:end-1) = trial_prices(:,2:end);
  %% Update price
	% Note that we use 'quadratic' pricing policy, so that the demand can be offloaded.
	% Adjust the step (delta_price) according to how much the capacity constraints have been
	% violated. Thus, we let the increase amount associated with the utilization ratio. 
	violate_ratio = output.loads./capacities;
  if b_violate
    % (a) update those violated;
    % To not influence the demand of flows that do not traverse bottleneck resources, we
    % only update the prices of bottleneck resources.  
    %   If the resource is over provisioned, the multiplier is larger than 1.
    violates = output.violates;
    delta_price(violates) = delta_price(violates).*(1+violate_ratio(violates));
    trial_prices(violates, end) = trial_prices(violates, end) + delta_price(violates);
  else
    % (b) update all except idle ones.
    % Increase the prices of all resources, based on the resource utilization.
    %   If the capacity tends to infinity, |delta_price| stay the same, while |prices|
    %   still increases in a constant rate. 
    %   We only increase the price of those resources that are utilized (resource
    %   utilization θ>0), since increasing the price of idle resources will not increase
    %   the profit of SP. 
    res_id = output.loads>0;
    delta_price = delta_price.*(1+violate_ratio);
    trial_prices(res_id, end) = trial_prices(res_id, end) + delta_price(res_id);
  end
	[sp_profit(3), b_violate, output] = this.SolveSCPCC(slices, trial_prices(:,end), options); 
	RECORD_SUBROUTINE_TIME;
	number_iters = number_iters + 1;
  
	if b_initial_trial 
    if sp_profit(2) >= sp_profit(3)
      % initial price is too high, half the initial price
      trial_prices = [zeros(num_resource,2), trial_prices(:,end-1)/2]; sp_profit(1) = -inf;
      t0 = t0/2;
      delta_price = t0 * unitcosts;  % init_price.link
      [sp_profit(3), b_violate, output] = this.SolveSCPCC(slices, trial_prices(:,end), options);
      RECORD_SUBROUTINE_TIME;
      number_iters = number_iters + 1;
      continue;
    else
      b_initial_trial= false;
    end
  end
	%% Stop condition
	% if the profit of SP is non-increasing, or the profit ratio reaches the predefined
	% threshold, then no need to further increase the resource prices.
	if b_profit_ratio && this.checkProfitRatio(trial_prices(:,end), options)
    b_forced_break = true;
    break;
	elseif sp_profit(3) <= sp_profit(2)
		break;
  end
end
%%
% If the resource capacity is not violated with the offered price, we serach a
% near-optimal price between [ρ_{k-1}, ρ_k]. 
% If previous procedure is forced to stopped by the profit ratio, we need to further
% search a higher price beyond ρ_k. Thus the following branch is executed only when SP
% profit start declining.
if ~b_violate && ~(b_profit_ratio && b_forced_break)
  sp_profit_N = zeros(2,1);
  epsilon = 10^-3;
  % a_L = 1; a_R = norm(prices(:,2))/norm(prices(:,1));
  a_L = 0; a_R = 1; delta_price = trial_prices(:,end) - trial_prices(:,end-2);
  while abs((sp_profit(2)-sp_profit(1))/sp_profit(1)) >= epsilon || ...
		abs((sp_profit(2)-sp_profit(3))/sp_profit(3)) >= epsilon
    a = [2/3*a_L+1/3*a_R; 1/3*a_L+2/3*a_R];
    for i = 2:-1:1
      % trial_prices = prices(:,1) * a(i);
      prices = trial_prices(:,1) + a(i)*delta_price;
			[sp_profit_N(i), b_violate, output] = this.SolveSCPCC(slices, prices, options);
			RECORD_SUBROUTINE_TIME;
    end
    number_iters = number_iters + 2;
    if sp_profit_N(1) > sp_profit_N(2)
      a_R = a(2); sp_profit(3) = sp_profit_N(2); sp_profit(2) = sp_profit_N(1);
    else
      a_L = a(1); sp_profit(1) = sp_profit_N(1); sp_profit(2) = sp_profit_N(2);
    end
  end
  trial_prices = prices; % b_violate = b_violate(1);
  sp_profit = sp_profit(2:3);
end
%%
% |delta_price.Link| and |delta_price.Node| of the first step can still be used, to
% improve the convergeNdce rate. Alternatively, one can reset the two vectors as follows
%
%    delta_price = t0 * unitcosts;         % init_price.link
if b_violate
  trial_prices = trial_prices(:,end);
	while b_violate
    % In this stage, increasing price leads to decline of profit. So only increase the
    % prices of violated resources. (Try not influence those who not traverse the
    % bottleneck resources.). 
    violates = output.violates;
    violate_ratio = output.loads(output.violates)./capacities(output.violates);
		delta_price(violates) = delta_price(violates) .* violate_ratio;     % {2|(k+1)/k}
		trial_prices(violates) = trial_prices(violates) + delta_price(violates);
		% Slices solve P1 with $ρ_k$, return the node (link) load v(y);
		% announce the resource price and optimize each network slice
    sp_profit(1) = sp_profit(2);
		[sp_profit(2), b_violate, output] = this.SolveSCPCC(slices, trial_prices, options);
		RECORD_SUBROUTINE_TIME;
		number_iters = number_iters + 1;
	end

  delta_price = t0 * trial_prices;  % 0.01 * init_price.link
  min_delta_price = delta_price;
  d0 = 10^-1;
  d1 = 10^-0;
  stop_cond_delta = ~isempty(find(delta_price > d0*unitcosts, 1));
  if b_profit_ratio
    stop_cond_profit = this.checkProfitRatio(trial_prices, options);
  else
    stop_cond_profit = true;
  end
  partial_violates = false(num_resource, 1);
  b_first = true;
  while stop_cond_delta && stop_cond_profit
    if DEBUG
			disp(table(trial_prices, delta_price, 'VariableNames', {'Prices', 'Delta Prices'}));
    end
    b_res = trial_prices > delta_price;
    trial_prices(b_res) = trial_prices(b_res) - delta_price(b_res);
    sp_profit(1) = sp_profit(2);
    [sp_profit(2), b_violate, output] = this.SolveSCPCC(slices, trial_prices, options);
		RECORD_SUBROUTINE_TIME;
    number_iters = number_iters + 1;
    
    if b_profit_ratio
      % the profit ratio of SP should not less than the predefined threshold.
      stop_cond_profit = this.checkProfitRatio(trial_prices, options);
    else
      % we decrease the price, the profit of SP should increase.
      stop_cond_profit = sp_profit(2) >= sp_profit(1);
    end
    if ~b_violate && stop_cond_profit  % no violate link/node   && 
      if b_first
        delta_price = delta_price * 2;
      else
        delta_price = delta_price + min_delta_price;
      end
      partial_violates = false(num_resource, 1);
    else
      b_first = false;
      trial_prices(b_res) = trial_prices(b_res) + delta_price(b_res);  % recover
      if ~stop_cond_profit && ~b_violate
        [~, ~, output] = this.SolveSCPCC(slices, trial_prices, options); %#ok<ASGLU>
				RECORD_SUBROUTINE_TIME;
        number_iters = number_iters + 1;
        break;
      end
      %%%
      %  If $\Delta_\rho$ has been smaller than the initial step, then only those
      %  resources with residual capacity will continue reduce their price, i.e. the
      %  components of step $\Delta_\rho$ corresponding to those overloaded
      %  resources is set to 0.
      assert_delta = isempty(find(delta_price > d1*unitcosts, 1));		% the vector is less than a threshold
      if assert_delta
        partial_violates = partial_violates | output.violates;
        delta_price(partial_violates) = 0;
      else
        partial_violates = false(num_resource, 1);
      end
      delta_price = delta_price / 2;
      min_delta_price = min(delta_price/4, min_delta_price);
    end
    %     stop_cond1 = norm(delta_price) > norm(10^-4 * unitcosts);
    stop_cond_delta = ~isempty(find(delta_price >d0*unitcosts, 1));
  end
end

%% Finalize substrate network
% # The resource allocation variables, virtual node/link load, and flow rate of each
% slice.
% # After the optimization, each network slice has record the final prices.
% # Record the substrate network's node/link load, price.
if options.bCountTime
	t_start = tic;
end
this.finalize(this.convertParameter(trial_prices), slices);
if options.bCountTime
	t_stop = toc(t_start); prt = prt + t_stop/Ns; srt = srt + t_stop;
	this.op.runtime = struct('Parallel', prt, 'Serial', srt);
  this.op.iterations = number_iters;
end

% Calculate the output
[varargout{1:nargout}] = this.calculateOutput(slices, fields, ...
	getstructfields(options, {'PricingPolicy'}));

% output the optimization results
if DEBUG
  fprintf('Optimization results:\n');
  fprintf('\tThe optimization procedure contains %d iterations.\n', number_iters);
  fprintf('\tNormalized network utilization is %G.\n', this.utilizationRatio());
end

end
