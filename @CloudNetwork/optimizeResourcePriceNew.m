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

% network data
NC = this.NumberDataCenters;
NS = length(slices);
NL = this.NumberLinks;
if nargin <= 2
    node_capacity = this.getDataCenterField('Capacity');
    link_capacity = this.getLinkField('Capacity');
else
    %% 
    % residual capacity + reallocatable capacity.
    node_capacity = this.getDataCenterField('ResidualCapacity');
    link_capacity = this.getLinkField('ResidualCapacity');
    for i = 1:NS
        sl = slices{i};
        node_capacity(sl.getDCPI) = node_capacity(sl.getDCPI) + ...
            sl.VirtualDataCenters.Capacity;
        link_capacity(sl.VirtualLinks.PhysicalLink) = ...
            link_capacity(sl.VirtualLinks.PhysicalLink) + sl.VirtualLinks.Capacity;
    end
end
link_uc = this.getLinkCost;
node_uc = this.getNodeCost;

%% Social-welfare aware price adjustment
% Initial Price
t1 = 1;           % {0.1|0.8|1}
if nargin >=2 && ~isempty(init_price)
    link_usage = false(this.NumberLinks,1);
    node_usage = false(this.NumberDataCenters,1);
    for i = 1:NS
        link_usage(slices{i}.VirtualLinks.PhysicalLink) = true;
        node_usage(slices{i}.getDCPI) = true;
    end
    if find([init_price.Link(link_usage)==0;init_price.Node(node_usage)==0],1)
        link_price = t1*link_uc;
        node_price = t1*node_uc;
        %         init_price.Link = link_price
        %         init_price.Node = node_price;
    else
        link_price = init_price.Link;
        node_price = init_price.Node;
    end
else
    link_price = t1* link_uc;
    node_price = t1* node_uc;
end
for i = 1:NS
    slices{i}.VirtualLinks{:,'Price'} = 0;
    slices{i}.VirtualDataCenters{:,'Price'} = 0;
end
t0 = 10^-1;     % {1|0.1|0.01}
delta_link_price = t0 * link_uc;  % init_price.link
delta_node_price = t0 * node_uc;

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
% if find(guess_init_link_price>link_price,1)
%     guess_init_link_price = link_price;
% end
% if find(guess_init_node_price>node_price,1)
%     guess_init_node_price = node_price;
% end
% link_price_prev = [guess_init_link_price, link_price];
% node_price_prev = [guess_init_node_price, node_price];
link_price_prev = [zeros(this.NumberLinks,1), link_price];
node_price_prev = [zeros(this.NumberDataCenters,1), node_price];
if b_profit_ratio
    b_forced_break = false;
end
new_net_welfare = SolveSCP(node_price, link_price);
sp_profit = this.getSliceProviderProfit(node_price, link_price, options);

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
    % If the capacity tends to infinity, |delta_link_price| and |delta_node_price| stay
    % the same, while |link_price| and |node_price| still increases in a constant rate.
    %
    % we only increase the price of those resources that are utilized (resource
    % utilization ¦È>0), since increasing the price of idle resources will not increase the
    % profit of SP.   
    [node_load, link_load] = this.getNetworkLoad(slices, 'sum');
    delta_link_price = delta_link_price.*(1+min(1,link_load./link_capacity));
    delta_node_price = delta_node_price.*(1+min(1,node_load./node_capacity));
    node_id = node_load>0;
    link_id = link_load>0;
    link_price(link_id) = link_price(link_id) + delta_link_price(link_id);
    node_price(node_id) = node_price(node_id) + delta_node_price(node_id);
    
    new_net_welfare = SolveSCP(node_price, link_price);
    sp_profit_new = this.getSliceProviderProfit(node_price, link_price, options);
    %% Stop condtion 
    % if the profit of SP is non-increasing, or the profit ratio reaches the predefined
    % threshold, then no need to further increase the resource prices.
    if sp_profit >= sp_profit_new
        %         if ~isempty(find(link_price_prev(:,1),1)) || ...
        %                 ~isempty(find(node_price_prev(:,1),1))
        if b_initial_trial
            % initial price is too high
            link_price_prev(:,2) = link_price_prev(:,2)/2;
            node_price_prev(:,2) = node_price_prev(:,2)/2;
            link_price = link_price_prev(:,2);
            node_price = node_price_prev(:,2);
            t0 = t0/2;
            delta_link_price = t0 * link_uc;  % init_price.link
            delta_node_price = t0 * node_uc;
            new_net_welfare = SolveSCP(node_price, link_price);
            sp_profit = this.getSliceProviderProfit(node_price, link_price, options);
            continue;
        else
            link_price_prev(:,1) = link_price_prev(:,2);
            link_price_prev(:,2) = link_price;
            node_price_prev(:,1) = node_price_prev(:,2);
            node_price_prev(:,2) = node_price;
            break;
        end
    end
    if b_initial_trial
        b_initial_trial = false;
    end
    sp_profit = sp_profit_new;
    link_price_prev(:,1) = link_price_prev(:,2);
    link_price_prev(:,2) = link_price;
    node_price_prev(:,1) = node_price_prev(:,2);
    node_price_prev(:,2) = node_price;
    if b_profit_ratio && this.checkProfitRatio(node_price, link_price, options)
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
        node_price_middle = [(2/3)*node_price_prev(:,1)+(1/3)*node_price_prev(:,2), ...
            (1/3)*node_price_prev(:,1)+(2/3)*node_price_prev(:,2)];
        link_price_middle = [(2/3)*link_price_prev(:,1)+(1/3)*link_price_prev(:,2), ...
            (1/3)*link_price_prev(:,1)+(2/3)*link_price_prev(:,2)];
        for i = 1:2
            new_net_welfare = SolveSCP(node_price_middle(:,i), link_price_middle(:,i));
            sp_profit_new(i) = this.getSliceProviderProfit(...
                node_price_middle(:,i), link_price_middle(:,i), ...
                getstructfields(options, {'Slices','PricingPolicy'}));
        end
        if sp_profit_new(1) > sp_profit_new(2)
            node_price_prev(:,2) = node_price_middle(:,2);
            link_price_prev(:,2) = link_price_middle(:,2);
        else
            node_price_prev(:,1) = node_price_middle(:,1);
            link_price_prev(:,1) = link_price_middle(:,1);
        end
        %%%
        % the stop condition can also be set as the difference of price.
        if abs((max(sp_profit_new)-sp_profit)/sp_profit) < epsilon
            break;
        else
            sp_profit = max(sp_profit_new);
        end
    end
    node_price = node_price_prev(:,1);       % temp_node_price
    link_price = link_price_prev(:,1);       % temp_link_price
end
%%
% |delta_link_price| and |delta_node_price| of the first step can still be used, to
% improve the convergence rate. Alternatively, one can reset the two vectors as follows
%
%    delta_link_price = t0 * link_uc;         % init_price.link
%    delta_node_price = t0 * node_uc;
k = 1;
while true
    %%% Compute the new resource price according to the resource consumption
    [node_load, link_load] = this.getNetworkLoad(slices, 'sum');
    b_link_violate = (link_capacity-link_load) < 1;
    b_node_violate = (node_capacity-node_load) < 1;
    if isempty(find(b_link_violate==1,1)) && isempty(find(b_node_violate==1,1))
        break;
    end
    link_price(b_link_violate)  = link_price(b_link_violate) + delta_link_price(b_link_violate);
    delta_link_price(b_link_violate) = delta_link_price(b_link_violate) .* ...
        (link_load(b_link_violate)./link_capacity(b_link_violate));     % {2|(k+1)/k}
    node_price(b_node_violate) = node_price(b_node_violate) + delta_node_price(b_node_violate);
    delta_node_price(b_node_violate) = delta_node_price(b_node_violate) .* ...
        (node_load(b_node_violate)./node_capacity(b_node_violate));     % {2|(k+1)/k}
    % Slices solve P1 with $¦Ñ_k$, return the node (link) load v(y);
    % announce the resource price and optimize each network slice
    number_iter = number_iter + 1;
    new_net_welfare = SolveSCP(node_price, link_price);
    k = k+1;
end
if ~isempty(DEBUG) && DEBUG
    fprintf('\tFirst stage objective value: %d.\n', new_net_welfare);
end

if k>1
    delta_link_price = t0 * link_price;  % 0.01 * init_price.link
    delta_node_price = t0 * node_price;
    min_delta_link_price = delta_link_price;
    min_delta_node_price = delta_node_price;
    d0 = 10^-1;
    d1 = 10^-0;
    stop_cond1 = ~isempty(find(delta_link_price > d0 * link_uc, 1));
    stop_cond2 = ~isempty(find(delta_node_price > d0 * node_uc, 1));
    if b_profit_ratio
        stop_cond3 = this.checkProfitRatio(node_price, link_price, ...
        getstructfields(options, {'PricingPolicy'}));
    else
        sp_profit = this.getSliceProviderProfit(node_price, link_price, ...
            getstructfields(options, {'Slices','PricingPolicy'}));
        stop_cond3 = true;
    end
    partial_link_violate = false(NL, 1);
    partial_node_violate = false(NC, 1);
    b_first = true;
    while stop_cond1 && stop_cond2 && stop_cond3
        number_iter = number_iter + 1;
        if ~isempty(DEBUG) && DEBUG
            disp('----link price    delta link price----')
            disp([link_price delta_link_price]);
        end
        b_link = link_price > delta_link_price;
        link_price(b_link) = link_price(b_link) - delta_link_price(b_link);
        if ~isempty(DEBUG) && DEBUG
            disp('----node price    delta node price----')
            disp([node_price delta_node_price]);
        end
        b_node = node_price > delta_node_price;
        node_price(b_node) = node_price(b_node) - delta_node_price(b_node);
        SolveSCP(node_price, link_price);
        [node_load, link_load] = this.getNetworkLoad(slices, 'sum');
        
        if b_profit_ratio
            % the profit ratio of SP should not less than the predefined threshold.
            stop_cond3 = this.checkProfitRatio(node_price, link_price, ...
                getstructfields(options, 'PricingPolicy'));
        else
            % we decrease the price, the profit of SP should increase.
            sp_profit_new = this.getSliceProviderProfit(node_price, link_price, ...
                getstructfields(options, {'Slices','PricingPolicy'}));
            stop_cond3 = sp_profit_new >= sp_profit;
        end
        b_link_violate = (link_capacity - link_load)<0;
        b_node_violate = (node_capacity - node_load)<0;
        assert_link_1 = isempty(find(b_link_violate==1,1));			% no violate link
        assert_node_1 = isempty(find(b_node_violate==1,1));			% no violate node
        if assert_link_1 && assert_node_1 && stop_cond3
            if b_first
                delta_link_price = delta_link_price * 2;
                delta_node_price = delta_node_price * 2;
            else
                delta_link_price = delta_link_price + min_delta_link_price;
                delta_node_price = delta_node_price + min_delta_node_price;
            end
            partial_link_violate = false(NL, 1);
            partial_node_violate = false(NC, 1);
        else
            b_first = false;
            link_price(b_link) = link_price(b_link) + delta_link_price(b_link);
            node_price(b_node) = node_price(b_node) + delta_node_price(b_node);
            if ~stop_cond3 && assert_link_1 && assert_node_1
                SolveSCP(node_price, link_price);
                break;
            end
            %%%
            %  If $\Delta_\rho$ has been smaller than the initial step, then only those
            %  resources with residual capacity will continue reduce their price, i.e. the
            %  components of step $\Delta_\rho$ corresponding to those overloaded
            %  resources is set to 0.   
            assert_link_2 = isempty(find(delta_link_price > d1 * link_uc, 1));		% the vector is less than a threshold
            assert_node_2 = isempty(find(delta_node_price > d1 * node_uc, 1));		% the vector is less than a threshold
            if assert_link_2
                partial_link_violate = partial_link_violate | b_link_violate;
                delta_link_price(partial_link_violate) = 0;
            else
                partial_link_violate = false(NL, 1);
            end
            delta_link_price = delta_link_price / 2;
            min_delta_link_price = min(delta_link_price/4, min_delta_link_price);
            if assert_node_2
                partial_node_violate = partial_node_violate | b_node_violate;
                delta_node_price(partial_node_violate) = 0;
            else
                partial_node_violate = false(NC, 1);
            end
            delta_node_price = delta_node_price / 2;
            min_delta_node_price = min(delta_node_price/4, min_delta_node_price);
        end
        %     stop_cond1 = norm(delta_link_price) > norm(10^-4 * link_uc);
        %     stop_cond2 = norm(delta_node_price) > norm(10^-4 * node_uc);
        stop_cond1 = ~isempty(find(delta_link_price > d0 * link_uc, 1));
        stop_cond2 = ~isempty(find(delta_node_price > d0 * node_uc, 1));
    end
end

%% Finalize substrate network
% # The resource allocation variables, virtual node/link load, and flow rate of each
% slice.
% # After the optimization, each network slice has record the final prices.
% # Record the substrate network's node/link load, price.
this.finalize(node_price, link_price, slices);

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
    function netprofit = SolveSCP(node_price_t, link_price_t)
        for s = 1:NS
            sl = slices{s};
            sl.prices.Link = link_price_t(sl.VirtualLinks.PhysicalLink);
            % |node_price| only contain the price of data center nodes.
            dc_id = sl.getDCPI;
            sl.prices.Node = node_price_t(dc_id);  
            if options.CountTime
                tic;
            end
            %%%
            % dynamically decide which version of _priceOptimalFlowRate_ to call.
            % See also <Slice.priceOptimalFlowRate> and
            % <DynamciSlice.priceOptimalFlowRate>.
            netprofit = sl.priceOptimalFlowRate([], options);
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