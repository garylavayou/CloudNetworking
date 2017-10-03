% Optimize resource price
% * *TODO* _refine the price adjustment algorithm._
% * *TODO* Resource Cost Model: linear, convex (quatratic)
function [output, runtime] = optimizeResourcePrice2(this, init_price)
global InfoLevel;
options = getstructfields(this.options, {'PricingFactor', 'PercentFactor'});

this.clearStates;
if nargout == 2
    slice_runtime = 0;
    runtime.Serial = 0;
    runtime.Parallel = 0;
end

% network data
NN = this.NumberNodes;
NS = this.NumberSlices;
NL = this.NumberLinks;
node_load = zeros(NN, NS);
link_load = zeros(NL, NS);
node_capacity = this.getNodeField('Capacity');
link_capacity = this.getLinkField('Capacity');
link_uc = this.getLinkCost;
node_uc = this.getNodeCost;

%% Social-welfare aware price adjustment
% Initial Price
t1 = 0.8;           % {0.1|1}
if nargin >=2 && ~isempty(init_price)
    link_price = t1 * init_price.Link;
    node_price = t1 * init_price.Node;
else
    init_price.Link = t1* link_uc;
    link_price = init_price.Link;
    init_price.Node = t1* node_uc;
    node_price = init_price.Node;
end
t0 = 10^-1;     % {1|0.1|0.01}
delta_link_price = t0 * link_uc;  % init_price.link
delta_node_price = t0 * node_uc;

number_iter = 1;
slice_profit = zeros(NS,1);
while true
    % Slices solve P1 with ¦Ñ_k, return the node (link) load v(y);
    % announce the resource price and optimize each network slice
    number_iter = number_iter + 1;
    for s = 1:NS
        this.slices{s}.VirtualLinks.Price = ...
            link_price(this.slices{s}.VirtualLinks.PhysicalLink);
        this.slices{s}.VirtualNodes.Price = ...
            node_price(this.slices{s}.VirtualNodes.PhysicalNode);
        if nargout == 2
            tic;
        end
        [slice_profit(s), node_load(:,s), link_load(:,s)] ...
            = this.slices{s}.priceOptimalFlowRate();
        if nargout == 2
            t = toc;
            slice_runtime = max(slice_runtime, t);
            runtime.Serial = runtime.Serial + t;
        end
    end
    if nargout == 2
        runtime.Parallel = runtime.Parallel + slice_runtime;
    end
    %%% Compute the new resource price according to the resource consumption
    [node_load, link_load] = this.getNetworkLoad(utility);
    b_link_violate = (link_capacity - link_load)<0;
    b_node_violate = (node_capacity - node_load)<0;
    if isempty(find(b_link_violate==1,1)) && isempty(find(b_node_violate==1,1))
        break;
    end
    link_price(b_link_violate)  = link_price(b_link_violate) + delta_link_price(b_link_violate);
    delta_link_price(b_link_violate) = delta_link_price(b_link_violate) * 2;
    node_price(b_node_violate) = node_price(b_node_violate) + delta_node_price(b_node_violate);
    delta_node_price(b_node_violate) = delta_node_price(b_node_violate) * 2;
end
if InfoLevel.UserModelDebug >= DisplayLevel.Notify
    fprintf('\tFirst stage objective value: %d.\n', new_net_welfare);
end

delta_link_price = t0 * link_price;  % 0.01 * init_price.link
delta_node_price = t0 * node_price;
min_delta_link_price = delta_link_price;
min_delta_node_price = delta_node_price;
d0 = 10^-1;
d1 = 10^-0;
stop_cond1 = ~isempty(find(delta_link_price > d0 * link_uc, 1));
stop_cond2 = ~isempty(find(delta_node_price > d0 * node_uc, 1));
stop_cond3 = this.checkProfitRatio(node_load, link_load, node_price, link_price, options);
slice_profit = zeros(NS,1);
partial_link_violate = false(NL, 1);
partial_node_violate = false(NN, 1);
b_first = true;
while stop_cond1 && stop_cond2 && stop_cond3
    number_iter = number_iter + 1;
    if InfoLevel.UserModelDebug == DisplayLevel.Iteration
        disp('----link price    delta link price----')
        disp([link_price delta_link_price]);
    end
    b_link = link_price > delta_link_price;
    link_price(b_link) = link_price(b_link) - delta_link_price(b_link);
    if InfoLevel.UserModelDebug == DisplayLevel.Iteration
        disp('----node price    delta node price----')
        disp([node_price delta_node_price]);
    end
    b_node = node_price > delta_node_price;
    node_price(b_node) = node_price(b_node) - delta_node_price(b_node);
    for s = 1:NS
        this.slices{s}.VirtualLinks.Price = ...
            link_price(this.slices{s}.VirtualLinks.PhysicalLink);
        this.slices{s}.VirtualNodes.Price = ...
            node_price(this.slices{s}.VirtualNodes.PhysicalNode);
        if nargout == 2
            tic;
        end
        [slice_profit(s), ~, ~] = this.slices{s}.priceOptimalFlowRate();
        if nargout == 2
            t = toc;
            slice_runtime = max(slice_runtime, t);
            runtime.Serial = runtime.Serial + t;
        end
    end
    if nargout == 2
        runtime.Parallel = runtime.Parallel + slice_runtime;
    end
    [node_load, link_load] = this.getNetworkLoad(utility);
    b_link_violate = (link_capacity - link_load)<0;
    b_node_violate = (node_capacity - node_load)<0;
    assert_link_1 = isempty(find(b_link_violate==1,1));			% no violate link
    assert_node_1 = isempty(find(b_node_violate==1,1));			% no violate node
    if assert_link_1 && assert_node_1
        if b_first
            delta_link_price = delta_link_price * 2;
            delta_node_price = delta_node_price * 2;
        else
            delta_link_price = delta_link_price + min_delta_link_price;
            delta_node_price = delta_node_price + min_delta_node_price;
        end
        partial_link_violate = false(NL, 1);
        partial_node_violate = false(NN, 1);
        stop_cond3 = this.checkProfitRatio(node_load, link_load, node_price, link_price, options);
    else
        b_first = false;
        link_price(b_link) = link_price(b_link) + delta_link_price(b_link);
        node_price(b_node) = node_price(b_node) + delta_node_price(b_node);
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
            partial_node_violate = false(NN, 1);
        end
        delta_node_price = delta_node_price / 2;
        min_delta_node_price = min(delta_node_price/4, min_delta_node_price);
    end
    %     stop_cond1 = norm(delta_link_price) > norm(10^-4 * link_uc);
    %     stop_cond2 = norm(delta_node_price) > norm(10^-4 * node_uc);
    stop_cond1 = ~isempty(find(delta_link_price > d0 * link_uc, 1));
    stop_cond2 = ~isempty(find(delta_node_price > d0 * node_uc, 1));
end

%% Profit ratio aware price adjustment
delta_link_price = t1 * link_price;
delta_node_price = t1 * node_price;
stop_cond3 = ~this.checkProfitRatio(node_load, link_load, node_price, link_price, options);
% only adjust price when the network's profit lower than the threshold.
if stop_cond3
    while stop_cond3
        number_iter = number_iter + 1;
        link_price = link_price + delta_link_price;
        node_price = node_price + delta_node_price;
        for s = 1:NS
            this.slices{s}.VirtualLinks.Price = ...
                link_price(this.slices{s}.VirtualLinks.PhysicalLink);
            this.slices{s}.VirtualNodes.Price = ...
                node_price(this.slices{s}.VirtualNodes.PhysicalNode);
            if nargout == 2
                tic;
            end
            this.slices{s}.priceOptimalFlowRate();
            if nargout == 2
                t = toc;
                slice_runtime = max(slice_runtime, t);
                runtime.Serial = runtime.Serial + t;
            end
        end
        if nargout == 2
            runtime.Parallel = runtime.Parallel + slice_runtime;
        end
        [node_load, link_load] = this.getNetworkLoad(utility);
        stop_cond3 = ~this.checkProfitRatio(node_load, link_load, node_price, link_price, options);
        if ~stop_cond3
            delta_link_price = delta_link_price * 2;
            delta_node_price = delta_node_price * 2;
        end
    end
    
    l = 0.5;
    h = 1;
    new_opts = options;
    new_opts.Epsilon = 10^-3;
    while this.checkProfitRatio(node_load, link_load, node_price, link_price, new_opts)
        number_iter = number_iter + 1;
        alpha = (l+h)/2;
        link_price = link_price + delta_link_price * alpha;
        node_price = node_price + delta_node_price * alpha;
        for s = 1:NS
            this.slices{s}.VirtualLinks.Price = ...
                link_price(this.slices{s}.VirtualLinks.PhysicalLink);
            this.slices{s}.VirtualNodes.Price = ...
                node_price(this.slices{s}.VirtualNodes.PhysicalNode);
            if nargout == 2
                tic;
            end
            this.slices{s}.priceOptimalFlowRate();
            if nargout == 2
                t = toc;
                slice_runtime = max(slice_runtime, t);
                runtime.Serial = runtime.Serial + t;
            end
        end
        if nargout == 2
            runtime.Parallel = runtime.Parallel + slice_runtime;
        end
        [node_load, link_load] = this.getNetworkLoad(utility);
        if this.checkProfitRatio(node_load, link_load, node_price, link_price, options)
            h = alpha;
        else
            l = alpha;
        end
    end
end

% Finalize substrate network
this.finalize(node_price, link_price);

% Calculate the output
output = this.calculateOutput();

% output the optimization results
if InfoLevel.UserModelDebug >= DisplayLevel.Final
    fprintf('Optimization results:\n');
    fprintf('\tThe optimization procedure contains %d iterations.\n', number_iter);
    fprintf('\tOptimal objective value: %d.\n', output_optimal.welfare_accurate);
    fprintf('\tNormalized network utilization is %G.\n', this.utilizationRatio);
end
end