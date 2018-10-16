% Optimize resource price
% * *TODO* _refine the price adjustment algorithm._
% * *TODO* Resource Cost Model: linear, convex (quatratic)
function [output, runtime] = optimizeResourcePrice2(this, init_price)
global DEBUG;
options = getstructfields(this.options, {'PricingFactor', 'PercentFactor'});
options.Stage = 'temp';

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
    prices.Link = t1 * init_price.Link;
    prices.Node = t1 * init_price.Node;
else
    init_price.Link = t1* link_uc;
    prices.Link = init_price.Link;
    init_price.Node = t1* node_uc;
    prices.Node = init_price.Node;
end
t0 = 10^-1;     % {1|0.1|0.01}
delta_price.Link = t0 * link_uc;  % init_price.link
delta_price.Node = t0 * node_uc;

number_iter = 1;
slice_profit = zeros(NS,1);
while true
    % Slices solve P1 with ��_k, return the node (link) load v(y);
    % announce the resource price and optimize each network slice
    number_iter = number_iter + 1;
    for s = 1:NS
        this.slices{s}.prices.Link = ...
            prices.Link(this.slices{s}.VirtualLinks.PhysicalLink);
        this.slices{s}.prices.Node = ...
            prices.Node(this.slices{s}.VirtualNodes.PhysicalNode);
        if nargout == 2
            tic;
        end
        [slice_profit(s), node_load(:,s), link_load(:,s)] ...
            = this.slices{s}.priceOptimalFlowRate([], options);
        if nargout == 2
            t = toc;
            slice_runtime = max(slice_runtime, t);
            runtime.Serial = runtime.Serial + t;
        end
        this.slices{s}.prices.Link = [];
        this.slices{s}.prices.Node = [];
    end
    if nargout == 2
        runtime.Parallel = runtime.Parallel + slice_runtime;
    end
    %%% Compute the new resource price according to the resource consumption
    load = this.getNetworkLoad([], options);
    b_link_violate = (link_capacity - load.Link)<0;
    b_node_violate = (node_capacity - load.Node)<0;
    if isempty(find(b_link_violate==1,1)) && isempty(find(b_node_violate==1,1))
        break;
    end
    prices.Link(b_link_violate)  = prices.Link(b_link_violate) + delta_price.Link(b_link_violate);
    delta_price.Link(b_link_violate) = delta_price.Link(b_link_violate) * 2;
    prices.Node(b_node_violate) = prices.Node(b_node_violate) + delta_price.Node(b_node_violate);
    delta_price.Node(b_node_violate) = delta_price.Node(b_node_violate) * 2;
end
if ~isempty(DEBUG) && DEBUG
    fprintf('\tFirst stage objective value: %d.\n', new_net_welfare);
end

delta_price.Link = t0 * prices.Link;  % 0.01 * init_price.link
delta_price.Node = t0 * prices.Node;
min_delta_link_price = delta_price.Link;
min_delta_node_price = delta_price.Node;
d0 = 10^-1;
d1 = 10^-0;
stop_cond1 = ~isempty(find(delta_price.Link > d0 * link_uc, 1));
stop_cond2 = ~isempty(find(delta_price.Node > d0 * node_uc, 1));
stop_cond3 = this.checkProfitRatio(prices, options);
slice_profit = zeros(NS,1);
partial_link_violate = false(NL, 1);
partial_node_violate = false(NN, 1);
b_first = true;
while stop_cond1 && stop_cond2 && stop_cond3
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
    for s = 1:NS
        this.slices{s}.prices.Link = ...
            prices.Link(this.slices{s}.VirtualLinks.PhysicalLink);
        this.slices{s}.prices.Node = ...
            prices.Node(this.slices{s}.VirtualNodes.PhysicalNode);
        if nargout == 2
            tic;
        end
        [slice_profit(s), ~, ~] = this.slices{s}.priceOptimalFlowRate([], options);
        if nargout == 2
            t = toc;
            slice_runtime = max(slice_runtime, t);
            runtime.Serial = runtime.Serial + t;
        end
        this.slices{s}.prices.Link = [];
        this.slices{s}.prices.Node = [];
    end
    if nargout == 2
        runtime.Parallel = runtime.Parallel + slice_runtime;
    end
    load = this.getNetworkLoad([], options);
    b_link_violate = (link_capacity - load.Link)<0;
    b_node_violate = (node_capacity - load.Node)<0;
    assert_link_1 = isempty(find(b_link_violate==1,1));			% no violate link
    assert_node_1 = isempty(find(b_node_violate==1,1));			% no violate node
    if assert_link_1 && assert_node_1
        if b_first
            delta_price.Link = delta_price.Link * 2;
            delta_price.Node = delta_price.Node * 2;
        else
            delta_price.Link = delta_price.Link + min_delta_link_price;
            delta_price.Node = delta_price.Node + min_delta_node_price;
        end
        partial_link_violate = false(NL, 1);
        partial_node_violate = false(NN, 1);
        stop_cond3 = this.checkProfitRatio(prices, options);
    else
        b_first = false;
        prices.Link(b_link) = prices.Link(b_link) + delta_price.Link(b_link);
        prices.Node(b_node) = prices.Node(b_node) + delta_price.Node(b_node);
        assert_link_2 = isempty(find(delta_price.Link > d1 * link_uc, 1));		% the vector is less than a threshold
        assert_node_2 = isempty(find(delta_price.Node > d1 * node_uc, 1));		% the vector is less than a threshold
        if assert_link_2
            partial_link_violate = partial_link_violate | b_link_violate;
            delta_price.Link(partial_link_violate) = 0;
        else
            partial_link_violate = false(NL, 1);
        end
        delta_price.Link = delta_price.Link / 2;
        min_delta_link_price = min(delta_price.Link/4, min_delta_link_price);
        if assert_node_2
            partial_node_violate = partial_node_violate | b_node_violate;
            delta_price.Node(partial_node_violate) = 0;
        else
            partial_node_violate = false(NN, 1);
        end
        delta_price.Node = delta_price.Node / 2;
        min_delta_node_price = min(delta_price.Node/4, min_delta_node_price);
    end
    %     stop_cond1 = norm(delta_price.Link) > norm(10^-4 * link_uc);
    %     stop_cond2 = norm(delta_price.Node) > norm(10^-4 * node_uc);
    stop_cond1 = ~isempty(find(delta_price.Link > d0 * link_uc, 1));
    stop_cond2 = ~isempty(find(delta_price.Node > d0 * node_uc, 1));
end

%% Profit ratio aware price adjustment
delta_price.Link = t1 * prices.Link;
delta_price.Node = t1 * prices.Node;
stop_cond3 = ~this.checkProfitRatio(prices, options);
% only adjust price when the network's profit lower than the threshold.
if stop_cond3
    while stop_cond3
        number_iter = number_iter + 1;
        prices.Link = prices.Link + delta_price.Link;
        prices.Node = prices.Node + delta_price.Node;
        for s = 1:NS
            this.slices{s}.prices.Link = ...
                prices.Link(this.slices{s}.VirtualLinks.PhysicalLink);
            this.slices{s}.prices.Node = ...
                prices.Node(this.slices{s}.VirtualNodes.PhysicalNode);
            if nargout == 2
                tic;
            end
            this.slices{s}.priceOptimalFlowRate([], options);
            if nargout == 2
                t = toc;
                slice_runtime = max(slice_runtime, t);
                runtime.Serial = runtime.Serial + t;
            end
            this.slices{s}.prices.Link = [];
            this.slices{s}.prices.Node = [];
        end
        if nargout == 2
            runtime.Parallel = runtime.Parallel + slice_runtime;
        end
        load = this.getNetworkLoad([], options);
        stop_cond3 = ~this.checkProfitRatio(prices, options);
        if ~stop_cond3
            delta_price.Link = delta_price.Link * 2;
            delta_price.Node = delta_price.Node * 2;
        end
    end
    
    l = 0.5;
    h = 1;
    new_opts = options;
    new_opts.Epsilon = 10^-3;
    while this.checkProfitRatio(prices, new_opts)
        number_iter = number_iter + 1;
        alpha = (l+h)/2;
        prices.Link = prices.Link + delta_price.Link * alpha;
        prices.Node = prices.Node + delta_price.Node * alpha;
        for s = 1:NS
            this.slices{s}.prices.Link = ...
                prices.Link(this.slices{s}.VirtualLinks.PhysicalLink);
            this.slices{s}.prices.Node = ...
                prices.Node(this.slices{s}.VirtualNodes.PhysicalNode);
            if nargout == 2
                tic;
            end
            this.slices{s}.priceOptimalFlowRate([], options);
            if nargout == 2
                t = toc;
                slice_runtime = max(slice_runtime, t);
                runtime.Serial = runtime.Serial + t;
            end
            this.slices{s}.prices.Link = [];
            this.slices{s}.prices.Node = [];
        end
        if nargout == 2
            runtime.Parallel = runtime.Parallel + slice_runtime;
        end
        load = this.getNetworkLoad([], options);
        if this.checkProfitRatio(prices, options)
            h = alpha;
        else
            l = alpha;
        end
    end
end

% Finalize substrate network
this.finalize(prices);

% Calculate the output
output = this.calculateOutput();

% output the optimization results
if ~isempty(DEBUG) && DEBUG
    fprintf('Optimization results:\n');
    fprintf('\tThe optimization procedure contains %d iterations.\n', number_iter);
    fprintf('\tOptimal objective value: %d.\n', output_optimal.welfare_accurate);
    fprintf('\tNormalized network utilization is %G.\n', this.utilizationRatio);
end
end