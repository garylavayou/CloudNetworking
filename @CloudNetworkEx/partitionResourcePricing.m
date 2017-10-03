function [output, runtime] = partitionResourcePricing(this, init_price)
global InfoLevel;
options = getstructfields(this.options, ...
    {'Method', 'ProfitType', 'WelfareType', 'PercentFactor'});

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
node_capacity = this.getNodeField('Capacity');
link_capacity = this.getLinkField('Capacity');
link_uc = this.getLinkCost;  % dynamic and static unit cost
node_uc = this.getNodeCost;

%% Social-welfare aware price adjustment
% Initial Price
t1 = 1;           % {0.1|1}
if nargin >= 2 && ~isempty(init_price)
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
part_factor = 1;

net_profit = zeros(NS,1);
b_violate = false;
while true
    number_iter = number_iter + 1;

    sliceOptimization;
    
    %% Price adjustment    
    % Slices solve P1 with ¦Ñ_k, return the node (link) load v(y);
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
    link_price(b_link_violate)  = link_price(b_link_violate) + delta_link_price(b_link_violate);
    delta_link_price(b_link_violate) = delta_link_price(b_link_violate) * 2;
    node_price(b_node_violate) = node_price(b_node_violate) + delta_node_price(b_node_violate);
    delta_node_price(b_node_violate) = delta_node_price(b_node_violate) * 2;
end

%% Stage 1-2
% delta_link_price = t0 * link_price;  % 0.01 * init_price.link
% delta_node_price = t0 * node_price;
step_weight = 1./link_price / min(1./link_price);
delta_link_price = step_weight.*link_uc;  
step_weight = 1./node_price / min(1./node_price);
delta_node_price = step_weight.*node_uc;
min_delta_link_price = delta_link_price;
min_delta_node_price = delta_node_price;
d0 = 10^-2;
d1 = 10^-1;
stop_cond1 = ~isempty(find(delta_link_price > d0 * link_uc, 1));
stop_cond2 = ~isempty(find(delta_node_price > d0 * node_uc, 1));
stop_cond3 = this.checkProfitRatio(aggr_node_load, aggr_link_load, node_price, link_price, options);
partial_link_violate = false(NL, 1);
partial_node_violate = false(NN, 1);
b_first = true;
while stop_cond1 && stop_cond2 && stop_cond3
    if InfoLevel.UserModelDebug == DisplayLevel.Iteration
        disp('----link price    delta link price----')
    end
    disp([link_price delta_link_price]);
    b_link = link_price > delta_link_price;
    link_price(b_link) = link_price(b_link) - delta_link_price(b_link);
    if InfoLevel.UserModelDebug == DisplayLevel.Iteration
        disp('----node price    delta node price----')
    end
    disp([node_price delta_node_price]);
    b_node = node_price > delta_node_price;
    node_price(b_node) = node_price(b_node) - delta_node_price(b_node);

    sliceOptimization;
    
    b_link_violate = (link_capacity - aggr_link_load)<0;
    b_node_violate = (node_capacity - aggr_node_load)<0;   
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
        stop_cond3 = this.checkProfitRatio(aggr_node_load, aggr_link_load, node_price, link_price, options);
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
    
    partitionWeightUpdate;
    
    partitionFactorUpdate;

    %     stop_cond1 = norm(delta_link_price) > norm(10^-4 * link_uc);
    %     stop_cond2 = norm(delta_node_price) > norm(10^-4 * node_uc);
    stop_cond1 = ~isempty(find(delta_link_price > d0 * link_uc, 1));
    stop_cond2 = ~isempty(find(delta_node_price > d0 * node_uc, 1));
    number_iter = number_iter + 1;
end

%% Stage 2
%% Profit ratio aware price adjustment
delta_link_price = t1 * link_price;
delta_node_price = t1 * node_price;
stop_cond3 = ~this.checkProfitRatio(aggr_node_load, aggr_link_load, node_price, link_price, options);
% only adjust price when the network's profit lower than the threshold.
if stop_cond3
    while stop_cond3
        link_price = link_price + delta_link_price;
        node_price = node_price + delta_node_price;
        
        sliceOptimization;
        
        stop_cond3 = ~this.checkProfitRatio(aggr_node_load, aggr_link_load, node_price, link_price, options);
        if ~stop_cond3
            delta_link_price = delta_link_price * 2;
            delta_node_price = delta_node_price * 2;
        end
        
        partitionWeightUpdate;
        
        partitionFactorUpdate;
        
    end
    
    l = 0.5;
    h = 1;
    new_opts = options;
    new_opts.Epsilon = 10^-3;
    while this.checkProfitRatio(aggr_node_load, aggr_link_load, node_price, link_price, new_opts)
        alpha = (l+h)/2;
        link_price = link_price + delta_link_price * alpha;
        node_price = node_price + delta_node_price * alpha;
        
        sliceOptimization;
        
        if this.checkProfitRatio(aggr_node_load, aggr_link_load, node_price, link_price, options)
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
    sl.Variables.x = sl.x_path;
    sl.Variables.z = sl.z_npf;
    sl.setPathBandwidth;
    sl.VirtualNodes.Load = sl.getNodeLoad;
    sl.VirtualLinks.Load = sl.getLinkLoad;
    sl.FlowTable.Rate = sl.getFlowRate;
end
this.setLinkField('Price', link_price);
this.setNodeField('Price', node_price);
this.setNodeField('Load', aggr_node_load);
this.setLinkField('Load', aggr_link_load);

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
output.node_price = node_price;
output.link_price = link_price;
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
    p = Slice.fcnNetProfit(var_x, sl);
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
    p = Slice.fcnNetProfit(var_x, sl, options);
    if this.static_factor ~= 0
        idx = sl.VirtualNodes.Load>0;
        nid = sl.VirtualNodes.PhysicalNode;
        p = p - dot(sl.VirtualNodes.Load(idx)./this.getNodeField('Load', nid(idx)),...
            this.getNodeField('StaticCost', nid(idx)));
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
    output.profit.ApproximatePrice(s) = Slice.fcnProfit(var_x, sl);
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
if InfoLevel.UserModelDebug >= DisplayLevel.Final
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
                    part_factor*partition_ratio.*this.getLinkField('Capacity', violate_link_id);
                sb_node_violate = b_node_violate(node_id);
                violate_node_id = node_id(sb_node_violate);
                partition_ratio = node_load(violate_node_id,si) ./ aggr_node_load(violate_node_id);
                sl.VirtualNodes.Capacity(sb_node_violate) = ...
                    part_factor*partition_ratio.*this.getNodeField('Capacity', violate_node_id);
            end
            sl.VirtualLinks.Price = link_price(link_id);
            sl.VirtualNodes.Price = node_price(node_id);
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
        end
        if nargout == 2
            runtime.Parallel = runtime.Parallel + slice_runtime;
        end
        aggr_node_load = sum(node_load, 2);
        aggr_link_load = sum(link_load, 2);
    end
end