% Optimize resource price
% * See also the <optimizeResourcePriceNew>, for more details in refining the price
% adjustment algorithm.
% * *TODO* Resource Cost Model: linear, convex (quatratic)
function [output, runtime] = optimizeResourcePrice(this, init_price)
global DEBUG;
options.Threshold = this.options.Threshold;
options.PricingPolicy = 'quadratic-price';
% this.clearStates;
if nargout == 2
    slice_runtime = 0;
    runtime.Serial = 0;
    runtime.Parallel = 0;
    options.CountTime = true;
else
    options.CountTime = false;
end

% network data
NC = this.NumberDataCenters;
NS = this.NumberSlices;
NL = this.NumberLinks;
node_capacity = this.getDataCenterField('Capacity');
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

%% Increase price to limit the slice's resource demand
number_iter = 1;
while true
    % Slices solve P1 with ¦Ñ_k, return the node (link) load v(y);
    % announce the resource price and optimize each network slice
    number_iter = number_iter + 1;
    SolveSCP(node_price, link_price);
    %%% Compute the new resource price according to the resource consumption
    [node_load, link_load] = this.getNetworkLoad([], 'sum');
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
if ~isempty(DEBUG) && DEBUG
    fprintf('\tFirst stage objective value: %d.\n', new_net_welfare);
end

%% decrease price to improve the net social welfare
delta_link_price = t0 * link_price;  % 0.01 * init_price.link
delta_node_price = t0 * node_price;
min_delta_link_price = delta_link_price;
min_delta_node_price = delta_node_price;
d0 = 10^-1;
d1 = 10^-0;
stop_cond1 = ~isempty(find(delta_link_price > d0 * link_uc, 1));
stop_cond2 = ~isempty(find(delta_node_price > d0 * node_uc, 1));
stop_cond3 = this.checkProfitRatio(node_price, link_price, options);
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
    
    [node_load, link_load] = this.getNetworkLoad([], 'sum');
    stop_cond3 = this.checkProfitRatio(node_price, link_price, options);
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
            stop_cond3 = this.checkProfitRatio(node_price, link_price, options);
            break;
        end
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

%% Profit ratio aware price adjustment
% only adjust price when the network's profit lower than the threshold.
if ~stop_cond3
    delta_link_price = t1 * link_price;
    delta_node_price = t1 * node_price;
    while ~stop_cond3
        number_iter = number_iter + 1;
        link_price = link_price + delta_link_price;
        node_price = node_price + delta_node_price;
        SolveSCP(node_price, link_price);
        stop_cond3 = this.checkProfitRatio(node_price, link_price, options);
        if ~stop_cond3
            delta_link_price = delta_link_price * 2;
            delta_node_price = delta_node_price * 2;
        end
    end
    
    l = -1;
    h = 0;
    new_opt = options;
    new_opt.Epsilon = 10^-3;
    while true
        number_iter = number_iter + 1;
        alpha = (l+h)/2;
        temp_link_price = link_price + delta_link_price * alpha;
        temp_node_price = node_price + delta_node_price * alpha;
        SolveSCP(temp_node_price, temp_link_price);
        [b_profit, profit_gap] = ...
            this.checkProfitRatio(temp_node_price, temp_link_price, new_opt);
        if b_profit || h-l<10^-4
            % due to precision error, the second condition is used to step out.
            link_price = temp_link_price;
            node_price = temp_node_price;
            if h-l < 10^-4
                if ~isempty(DEBUG) && DEBUG
                    warning('precision error: %.4f', profit_gap);
                end
            end
            break;
        end
        if this.checkProfitRatio(temp_node_price, temp_link_price, options)
            h = alpha;
        else
            l = alpha;
        end
    end
end

%% Finalize substrate network
% # The resource allocation variables, virtual node/link load, and flow rate of each
% slice.
% # After the optimization, each network slice has record the final prices.
% # Record the substrate network's node/link load, price.
this.finalize(node_price, link_price);

% Calculate the output
output = this.calculateOutput([], options);

% output the optimization results
if ~isempty(DEBUG) && DEBUG
    fprintf('Optimization results:\n');
    fprintf('\tThe optimization procedure contains %d iterations.\n', number_iter);
    fprintf('\tOptimal objective value: %d.\n', output_optimal.welfare_accurate);
    fprintf('\tNormalized network utilization is %G.\n', this.utilizationRatio);
end

%%
    function SolveSCP(node_price_t, link_price_t)
        for s = 1:NS
            sl = this.slices{s};
            sl.prices.Link = link_price_t(sl.VirtualLinks.PhysicalLink);
            dc_id = sl.getDCPI;
            sl.prices.Node = node_price_t(dc_id);
            if options.CountTime
                tic;
            end
            this.slices{s}.priceOptimalFlowRate([], options);
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