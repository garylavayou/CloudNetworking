%% Static Network Slicing
% In the static slicing method, once the resource is allocated to a slice, the allocation
% scheme is not changed during its lifetime.
%
% *FIXME*: when link and node resources are exhausted, some slice request might be
% rejected.
%%
function output = staticSlicing(this, slice, options)

if nargin <= 2
    options.Display = 'final';
    options.Method = 'slice';
else
    if ~isfield(options, 'Display')
        options.Display = 'final';
    end
    if ~isfield(options, 'Method')
        options.Method = 'slice';
    end
end

NS = this.NumberSlices;
NN = this.NumberNodes;
link_uc = this.getLinkField('UnitCost');
node_uc = this.getNodeField('UnitCost');
epsilon = this.unitStaticNodeCost;
phis_n = epsilon*this.delta*(NN-1)/this.totalNodeCapacity;
phis_l = epsilon*(1-this.delta)*(NN-1)/this.totalLinkCapacity;
link_price = (link_uc + phis_l) * (1 + options.PricingFactor);
node_price = (node_uc + phis_n) * (1 + options.PricingFactor);
this.setLinkField('Price', link_price);
this.setNodeField('Price', node_price);

if nargin>=2 && ~isempty(slice)
    %% Allocate Resource to the new arrival slice
    % The residual capacity of the substrate network is available to the slice.
    ss = Slice(slice);
    ss.VirtualNodes.Capacity = ...
        this.getNodeField('ResidualCapacity', ss.VirtualNodes.PhysicalNode);
    ss.VirtualLinks.Capacity = ...
        this.getLinkField('ResidualCapacity', ss.VirtualLinks.PhysicalLink);
    ss.optimalFlowRate(options);
    %% Finalize the new slice and the substrate network
    % # After the optimization, the resource allocation variables, flow rate, virtual
    % node/link load of the last slice have been recorded.
    % # Calculate and announce the resource prices to the new slice. The price is fixed in
    % the static slicing method, so the price has been calculated in advance.
    % # Record the substrate network's node/link load, price. When a slice arrive or
    % depart, the network load changes.
    slice.Variables = ss.Variables;
    slice.setPathBandwidth;
    slice.FlowTable.Rate = ss.FlowTable.Rate;
    slice.VirtualLinks.Load = ss.VirtualLinks.Load;
    slice.VirtualNodes.Load = ss.VirtualNodes.Load;
    slice.VirtualLinks.Price = link_price(slice.VirtualLinks.PhysicalLink);
    slice.VirtualNodes.Price = node_price(slice.VirtualNodes.PhysicalNode);
end

node_load = zeros(NN, NS);
link_load = zeros(this.NumberLinks, NS);
for s = 1:NS
    sl = this.slices{s};
    node_load(sl.VirtualNodes.PhysicalNode, s) = sl.VirtualNodes.Load;
    link_load(sl.VirtualLinks.PhysicalLink, s) = sl.VirtualLinks.Load;
end
aggr_node_load = sum(node_load, 2);
aggr_link_load = sum(link_load, 2);
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
    %      net_profit = (utility - approximate_cost) * pf,
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
    p = -Slice.fcnNetProfit(var_x, sl);
    output.welfare_approx = output.welfare_approx + p;
    output.welfare_accurate = output.welfare_accurate...
        + sl.weight*sum(log(sl.FlowTable.Rate));
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
    p = -Slice.fcnNetProfit(var_x, sl, options);
    idx = sl.VirtualNodes.Load>0;
    nid = sl.VirtualNodes.PhysicalNode;
    p = p - dot(sl.VirtualNodes.Load(idx)./this.getNodeField('Load', nid(idx)),...
        this.getNodeField('StaticCost', nid(idx)));
    output.profit.AccuratePercent(s) = options.PercentFactor * p;
    %%%
    % * *Net profit with offered price*
    %
    % For slices,
    %
    %      net_profit = utility - payment(price).
    %
    % The price is calculated accoding to resources cost and a pricing factor. The pricing
    % factor is configured to ensure that the substrate network can make profit.
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
    output.profit.ApproximatePrice(s) = -Slice.fcnProfit(var_x, sl);
end
f = (1-options.PercentFactor)/options.PercentFactor;
output.profit.ApproximatePercent(end) = ...
    sum(output.profit.ApproximatePercent(1:(end-1))*f);
output.profit.AccuratePercent(end) = ...
    sum(output.profit.AccuratePercent(1:(end-1))*f);
output.profit.ApproximatePrice(end) = output.welfare_approx - ...
    sum(output.profit.ApproximatePrice(1:(end-1)));
output.profit.AccuratePrice = output.profit.ApproximatePrice;
output.welfare_accurate = output.welfare_accurate - ...
    this.getNetworkCost([],[], options.Model);
output.profit.AccuratePrice(end) = output.welfare_accurate - ...
    sum(output.profit.AccuratePrice(1:(end-1)));

end

