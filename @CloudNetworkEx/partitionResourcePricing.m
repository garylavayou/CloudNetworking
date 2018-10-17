function [output, runtime] = partitionResourcePricing(this, init_price)
global DEBUG;
options = getstructfields(this.options, ...
    {'SlicingMethod', 'ProfitType', 'WelfareType', 'PercentFactor'});

this.clearStates;
if nargout == 2
    slice_runtime = 0;
    runtime.Serial = 0;
    runtime.Parallel = 0;
end

% network data
NN = this.NumberNodes;
NL = this.NumberLinks;
NS = this.NumberSlices;
node_load = zeros(NN, NS);
link_load = zeros(NL, NS);
aggr_link_load = zeros(NL,1);
aggr_node_load = zeros(NS,1);
node_capacity = this.readNode('Capacity');
link_capacity = this.readLink('Capacity');
link_uc = this.getLinkCost;  % dynamic and static unit cost
node_uc = this.getNodeCost;

%% Social-welfare aware price adjustment
% Initial Price
t1 = 1;           % {0.1|1}
if nargin >= 2 && ~isempty(init_price)
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
part_factor = 1;

net_profit = zeros(NS,1);
b_violate = false;
while true
    number_iter = number_iter + 1;

    sliceOptimization;
    
    %% Price adjustment    
    % Slices solve P1 with ��_k, return the node (link) load v(y);
    % announce the resource price and optimize each network slice
    % Compute the new resource price according to the resource consumption
    b_link_violate = (link_capacity - aggr_link_load)<0;
    b_node_violate = (node_capacity - aggr_node_load)<0;
    if isempty(find(b_link_violate==1,1)) && isempty(find(b_node_violate==1,1))
        b_violate = false;
        break;
    else
        b_violate = true;
    end
    prices.Link(b_link_violate)  = prices.Link(b_link_violate) + delta_price.Link(b_link_violate);
    delta_price.Link(b_link_violate) = delta_price.Link(b_link_violate) * 2;
    prices.Node(b_node_violate) = prices.Node(b_node_violate) + delta_price.Node(b_node_violate);
    delta_price.Node(b_node_violate) = delta_price.Node(b_node_violate) * 2;
end

%% Stage 1-2
% delta_price.Link = t0 * prices.Link;  % 0.01 * init_price.link
% delta_price.Node = t0 * prices.Node;
step_weight = 1./prices.Link / min(1./prices.Link);
delta_price.Link = step_weight.*link_uc;  
step_weight = 1./prices.Node / min(1./prices.Node);
delta_price.Node = step_weight.*node_uc;
min_delta_link_price = delta_price.Link;
min_delta_node_price = delta_price.Node;
d0 = 10^-2;
d1 = 10^-1;
stop_cond1 = ~isempty(find(delta_price.Link > d0 * link_uc, 1));
stop_cond2 = ~isempty(find(delta_price.Node > d0 * node_uc, 1));
stop_cond3 = this.checkProfitRatio(prices, options);
partial_link_violate = false(NL, 1);
partial_node_violate = false(NN, 1);
b_first = true;
while stop_cond1 && stop_cond2 && stop_cond3
    if ~isempty(DEBUG) && DEBUG
        disp('----link price    delta link price----')
    end
    disp([prices.Link delta_price.Link]);
    b_link = prices.Link > delta_price.Link;
    prices.Link(b_link) = prices.Link(b_link) - delta_price.Link(b_link);
    if ~isempty(DEBUG) && DEBUG
        disp('----node price    delta node price----')
    end
    disp([prices.Node delta_price.Node]);
    b_node = prices.Node > delta_price.Node;
    prices.Node(b_node) = prices.Node(b_node) - delta_price.Node(b_node);

    sliceOptimization;
    
    b_link_violate = (link_capacity - aggr_link_load)<0;
    b_node_violate = (node_capacity - aggr_node_load)<0;   
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
    
    partitionWeightUpdate;
    
    partitionFactorUpdate;

    %     stop_cond1 = norm(delta_price.Link) > norm(10^-4 * link_uc);
    %     stop_cond2 = norm(delta_price.Node) > norm(10^-4 * node_uc);
    stop_cond1 = ~isempty(find(delta_price.Link > d0 * link_uc, 1));
    stop_cond2 = ~isempty(find(delta_price.Node > d0 * node_uc, 1));
    number_iter = number_iter + 1;
end

%% Stage 2
%% Profit ratio aware price adjustment
delta_price.Link = t1 * prices.Link;
delta_price.Node = t1 * prices.Node;
stop_cond3 = ~this.checkProfitRatio(prices, options);
% only adjust price when the network's profit lower than the threshold.
if stop_cond3
    while stop_cond3
        prices.Link = prices.Link + delta_price.Link;
        prices.Node = prices.Node + delta_price.Node;
        
        sliceOptimization;
        
        stop_cond3 = ~this.checkProfitRatio(prices, options);
        if ~stop_cond3
            delta_price.Link = delta_price.Link * 2;
            delta_price.Node = delta_price.Node * 2;
        end
        
        partitionWeightUpdate;
        
        partitionFactorUpdate;
        
    end
    
    l = 0.5;
    h = 1;
    new_opts = options;
    new_opts.Epsilon = 10^-3;
    while this.checkProfitRatio(prices, new_opts)
        alpha = (l+h)/2;
        prices.Link = prices.Link + delta_price.Link * alpha;
        prices.Node = prices.Node + delta_price.Node * alpha;
        
        sliceOptimization;
        
        if this.checkProfitRatio(prices, options)
            h = alpha;
        else
            l = alpha;
        end
        
        partitionWeightUpdate;
        
        partitionFactorUpdate;
    end
end

%% Finalize substrate network
% # The resource allocation variables, virtual node/link load, and flow rate of each
% slice.
% # After the optimization, each network slice has record the final prices.
% # Record the substrate network's node/link load, price.
for s = 1:NS
    sl = this.slices{s};
    sl.Variables.x = sl.temp_vars.x;
    sl.Variables.z = sl.temp_vars.z;
    sl.setPathBandwidth;
    sl.VirtualNodes.Load = sl.getNodeLoad;
    sl.VirtualLinks.Load = sl.getLinkLoad;
    sl.FlowTable.Rate = sl.getFlowRate;
end
this.writeLink('Price', prices.Link);
this.writeDataCenter('Price', prices.Node);
this.writeDataCenter('Load', aggr_node_load);
this.writeLink('Load', aggr_link_load);

%% Calculate the output
% The output variables includes,
%
% # Price and Load of physical nodes and links;
% # The approximation of net social welfare, and the accurate net social welfare;
% # The net profit of each slice and the substrate network, including four results from
% different methods, i.e. |ApproximatePercent|, |ApproximatePrice|, |AccuratePercent|,
% |AccuratePrice|.
% # Flow rate of all flows in the network.
%% TODO: replace with calculateOutput
output.node_price = prices.Node;
output.link_price = prices.Link;
output.node_load = aggr_node_load;
output.link_load = aggr_link_load;
output.welfare_approx = 0;
output.welfare_accurate = 0;
t = zeros(NS+1, 1);
output.profit = table(t, t, t, t, 'VariableNames',...
    {'ApproximatePercent', 'ApproximatePrice', 'AccuratePercent', 'AccuratePrice'});
clear t;
output.flow_rate = [];
options.Model = 'Accurate';
options.bFinal = true;
%% Calculate the net social welfare
% The resource price is also keep fixed since it related to the linear resouce dynamic
% cost, static cost and the prcing factor. Here, we use two model to assess the net
% social welfare,
%
% * *Approximate Model*: the net social welfare is the sum of raw profit of all slices.
% * *Accurate Model*: the net social welfare is the total utility less the total network
% cost.
for s = 1:NS
    sl = this.slices{s};
    output.flow_rate = [output.flow_rate; sl.FlowTable.Rate];
    
    var_x = [sl.Variables.x; sl.Variables.z];
    %% Calculate the profit
    % * *Proportional Net profit with approximate cost*
    %
    % For slices,
    %
    %      slice_profit = (utility - approximate_cost) * pf,
    %
    % where the approximate cost is computed by our approximation formula. The
    % proportional net profit means only one part of the raw net profit is attributed to
    % the slice, another part should be deliver to the substrate network. As a result,
    % for substrate network,
    %
    %      net_profit = sum(slice_profit) * (1-pf).
    %
    % The net profit of substrate network is the sum of proportion profit from each slice
    % computed by the _Approximate Model_. It is the same computing method, when using
    % _Accurate Model_ with proportion profit.
    %
    % *NOTE*: fcnNetProfit here evaluate the raw profit with resource cost.
    p = SimpleSlice.fcnNetProfit(var_x, sl);
    output.welfare_approx = output.welfare_approx + p;
    output.welfare_accurate = output.welfare_accurate...
        + sl.weight*sum(fcnUtility(sl.FlowTable.Rate));
    output.profit.ApproximatePercent(s) = options.PercentFactor * p;
    %%%
    % * *Proportional Net profit with accurate cost*
    %
    % For slices,
    %
    %      net_profit = (utility - accurate_cost) * pf.
    %
    % The accurate cost is based on the actual resource allocation after the optimization.
    % Slices which share a physical resource will also proportionally share the static
    % cost of the resource.
    %
    % *NOTE*: it is not accurate to calculate the static cost of each slice with this
    % method.
    p = SimpleSlice.fcnNetProfit(var_x, sl, options);
    if this.static_factor ~= 0
        idx = sl.VirtualNodes.Load>0;
        nid = sl.VirtualNodes.PhysicalNode;
        p = p - dot(sl.VirtualNodes.Load(idx)./this.readNode('Load', nid(idx)),...
            this.readNode('StaticCost', nid(idx)));
    end
    output.profit.AccuratePercent(s) = options.PercentFactor * p;
    %%%
    % * *Net profit with offered price*
    %
    % For slices,
    %
    %      net_profit = utility - payment(price).
    %
    % The price is calculated accoding to optimization procedure.
    %
    % When the resource price is offered, the slice's profit is not related to how the
    % cost is computed. Therefore, |output.profit.ApproximatePrice| and
    % |output.profit.AccuratePercent| of slices is the same. The cost only affects the
    % substrate network's profit. For substrate network with _Approximate Model_ and
    % _Accurate Model_, the net profit is given by accordingly
    %
    %      net_profit = net_social_welfare - sum(slice_profit)
    %                 = sum(slice_payment(price))) - cost(approximate/accurate)
    %
    % *NOTE*: fcnProfit evalute the profit using offered price.
    output.profit.ApproximatePrice(s) = SimpleSlice.fcnProfit(sl,options);
end
if ~isempty(this.eta)
    embed_profit_approx = this.eta*this.getNetworkCost;
    embed_profit_accurate = this.getNetworkCost([], [], options.Model);
else
    embed_profit_approx = 0;
    embed_profit_accurate = 0;
end
output.welfare_approx = output.welfare_approx + embed_profit_approx;
output.welfare_accurate = output.welfare_accurate - ...
    this.getNetworkCost([],[], options.Model) + embed_profit_accurate;

f = (1-options.PercentFactor)/options.PercentFactor;
output.profit.ApproximatePercent(end) = ...
    sum(output.profit.ApproximatePercent(1:(end-1))*f) + embed_profit_approx;
output.profit.AccuratePercent(end) = ...
    sum(output.profit.AccuratePercent(1:(end-1))*f) + embed_profit_accurate;
output.profit.ApproximatePrice(end) = output.welfare_approx - ...
    sum(output.profit.ApproximatePrice(1:(end-1))) + embed_profit_approx;
output.profit.AccuratePrice = output.profit.ApproximatePrice;
output.profit.AccuratePrice(end) = output.welfare_accurate - ...
    sum(output.profit.AccuratePrice(1:(end-1))) + embed_profit_accurate;
% output the optimization results
if ~isempty(DEBUG) && DEBUG
    fprintf('Optimization results:\n');
    fprintf('\tThe optimization procedure contains %d iterations.\n', number_iter);
    fprintf('\tOptimal objective value: %d.\n', new_net_welfare);
    fprintf('\tNormalized network utilization is %G.\n', this.utilizationRatio);
end

%% Sub-procedure
    function sliceOptimization
        for si = 1:NS
            sl = this.slices{si};
            node_id = sl.VirtualNodes.PhysicalNode;
            link_id = sl.VirtualLinks.PhysicalLink;
            sl.VirtualLinks.Capacity = inf*ones(sl.NumberVirtualLinks,1);
            sl.VirtualNodes.Capacity = inf*ones(sl.NumberVirtualNodes,1);
            if b_violate
                sb_link_violate = b_link_violate(link_id);
                violate_link_id = link_id(sb_link_violate);
                partition_ratio = link_load(violate_link_id,si) ./ aggr_link_load(violate_link_id);
                sl.VirtualLinks.Capacity(sb_link_violate) = ...
                    part_factor*partition_ratio.*this.readLink('Capacity', violate_link_id);
                sb_node_violate = b_node_violate(node_id);
                violate_node_id = node_id(sb_node_violate);
                partition_ratio = node_load(violate_node_id,si) ./ aggr_node_load(violate_node_id);
                sl.VirtualNodes.Capacity(sb_node_violate) = ...
                    part_factor*partition_ratio.*this.readNode('Capacity', violate_node_id);
            end
            sl.prices.Link = prices.Link(link_id);
            sl.prices.Node = prices.Node(node_id);
            %%%
            % optimal each slice with price and resource constraints.
            if nargout == 2
                tic;
            end
            [net_profit(si), node_load(:,si), link_load(:,si)] = ...
                sl.optimalFlowRate(options);  % with capacity constraints
            if nargout == 2
                t = toc;
                slice_runtime = max(slice_runtime, t);
                runtime.Serial = runtime.Serial + t;
            end
            sl.prices.Link = [];
            sl.prices.Node = [];
        end
        if nargout == 2
            runtime.Parallel = runtime.Parallel + slice_runtime;
        end
        aggr_node_load = sum(node_load, 2);
        aggr_link_load = sum(link_load, 2);
    end
end