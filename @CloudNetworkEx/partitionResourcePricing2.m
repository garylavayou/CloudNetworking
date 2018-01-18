function [output, runtime] = partitionResourcePricing2(this, init_price)
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
node_capacity = this.getNodeField('Capacity');
link_capacity = this.getLinkField('Capacity');
link_uc = this.getLinkCost;  % dynamic and static unit cost
node_uc = this.getNodeCost;
%% Count the usage of links and nodes
% Here, the simplest method is adopted. If a link/node might be used by a slice, then its
% count is increased by the slice's weight.
% Another method: if a link/node might be used by a flow in a slice, then its count is
% increased by the slice's weight.
link_usage = zeros(NL, NS);
node_usage = zeros(NN, NS);
aggr_link_usage = zeros(NL,1);
aggr_node_usage = zeros(NS,1);
weighted_link_usage = zeros(NL, NS);
weighted_node_usage = zeros(NN, NS);
for s = 1:NS
    link_usage(this.slices{s}.VirtualLinks.PhysicalLink,s) = 1;
    node_usage(this.slices{s}.VirtualNodes.PhysicalNode,s) = 1;
end
partition_weight = struct('Node', ones(1, NS), 'Link', ones(1, NS));

%% Social-welfare aware price adjustment
% Initial Price
t1 = 0.5;           % {0.1|1}
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
part_factor = struct('Link', 1.1, 'Node', 1.1);
res_occupy_ratio = struct('Link', ones(NS,1), 'Node', ones(NS,1));

partitionWeightUpdate;

net_profit = zeros(NS,1);
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
        break;
    end
    link_price(b_link_violate)  = link_price(b_link_violate) + delta_link_price(b_link_violate);
    delta_link_price(b_link_violate) = delta_link_price(b_link_violate) * 2;
    node_price(b_node_violate) = node_price(b_node_violate) + delta_node_price(b_node_violate);
    delta_node_price(b_node_violate) = delta_node_price(b_node_violate) * 2;
    
    %% Update partition weight and factor
    % If all slices and the substrate network have residual capacity on all resources,
    % then the partition factor can be scale-down. Since the resource requirement on node
    % and link has different ratio among slices. So the partition factor may be different
    % between node and links.
    %
    % When update the partition weight, we use the resource occupation ratio to
    % represent the resource demand of each slice:
    % 1. calculate a single partition weight;
    % 2. calculate a partition weight for node and link separately (current);
    % 3. calculate a partition weight for each resource;
    %%% 
    partitionWeightUpdate;
    
    partitionFactorUpdate;
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
    if ~isempty(DEBUG) && DEBUG
        disp('----link price    delta link price----')
    end
    disp([link_price delta_link_price]);
    b_link = link_price > delta_link_price;
    link_price(b_link) = link_price(b_link) - delta_link_price(b_link);
    if ~isempty(DEBUG) && DEBUG
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
    sl.Variables.x = sl.temp_vars.x;
    sl.Variables.z = sl.temp_vars.z;
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

%% TODO replace with calculateOutput.
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
    output.profit.ApproximatePrice(s) = Slice.fcnProfit(sl, options);
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
            link_id = sl.VirtualLinks.PhysicalLink;
            partition_ratio = weighted_link_usage(link_id,si) ./ aggr_link_usage(link_id);
            sl.VirtualLinks.Capacity = ...
                part_factor.Link*partition_ratio.*this.getLinkField('Capacity', link_id);
            node_id = sl.VirtualNodes.PhysicalNode;
            partition_ratio = weighted_node_usage(node_id,si) ./ aggr_node_usage(node_id);
            sl.VirtualNodes.Capacity = ...
                part_factor.Node*partition_ratio.*this.getNodeField('Capacity', node_id);
            sl.prices.Link = link_price(link_id);
            sl.prices.Node = node_price(node_id);
            %%%
            % optimal each slice with price and resource constraints.
            if nargout == 2
                tic;
            end
            [net_profit(si), node_load(:,si), link_load(:,si)] = ...
                sl.optimalFlowRate(options);
            if nargout == 2
                t = toc;
                slice_runtime = max(slice_runtime, t);
                runtime.Serial = runtime.Serial + t;
            end
            res_occupy_ratio.Link(si) = sum(link_load(:,si))/sum(link_capacity(link_id));
            res_occupy_ratio.Node(si) = sum(node_load(:,si))/sum(node_capacity(node_id));
            sl.prices.Link = [];
            sl.prices.Node = [];
        end
        if nargout == 2
            runtime.Parallel = runtime.Parallel + slice_runtime;
        end
        aggr_node_load = sum(node_load, 2);
        aggr_link_load = sum(link_load, 2);
    end

    function partitionFactorUpdate
        b_slice_link_violate = false(NS,1);
        b_slice_node_violate = false(NS,1);
        for si = 1:NS
            sl = this.slices{si};
            link_id = sl.VirtualLinks.PhysicalLink;
            if find(sl.VirtualLinks.Capacity-link_load(link_id,si)<=eps,1)
                b_slice_link_violate(si) = true;
            end
            node_id = sl.VirtualNodes.PhysicalNode;
            if find(sl.VirtualNodes.Capacity-node_load(node_id,si)<=eps,1)
                b_slice_node_violate(si) = true;
            end
        end
        if isempty(find(b_slice_link_violate,1)) && isempty(find(b_slice_node_violate,1))
            slice_scale_factor = struct('Link', ones(NS,1), 'Node', ones(NS,1));
            for si = 1:NS
                sl = this.slices{si};
                link_id = sl.VirtualLinks.PhysicalLink;
                node_id = sl.VirtualNodes.PhysicalNode;
                slice_scale_factor.Link(si) = ...
                    min(sl.VirtualLinks.Capacity./link_load(link_id,si));
                slice_scale_factor.Node(si) = ...
                    min(sl.VirtualNodes.Capacity./node_load(node_id,si));
            end
            part_factor.Link = part_factor.Link/min(slice_scale_factor.Link);
            part_factor.Node = part_factor.Node/min(slice_scale_factor.Node);
        end
    end

    %% partition the network resource
    % the partition ratio is determined by the usage number
    function partitionWeightUpdate
        % Normalized partition weight
        partition_weight.Link = sqrt(res_occupy_ratio.Link'/min(res_occupy_ratio.Link));
        partition_weight.Node = sqrt(res_occupy_ratio.Node'/min(res_occupy_ratio.Node));
        weighted_link_usage = link_usage .* partition_weight.Link;  % compatable arithmetic operation
        weighted_node_usage = node_usage .* partition_weight.Node;  % compatable arithmetic operation
        aggr_link_usage = sum(weighted_link_usage, 2);
        aggr_node_usage = sum(weighted_node_usage, 2);
    end
end