%% Resource Partition Optimization
% link and node resources are firstly evenly or proportionaly partiotioned into each
% slice, and then each slice independently run optimization.
%
% the process is contrary to VNE, which known the explicit demand of network resources,
% and try to find the embedding solution.
function [ output ] = resourcePartitionOptimization( this, slice_weight, options )
if nargin <= 1 || isempty(slice_weight)
    slice_weight = ones(this.NumberSlices,1);
end
if nargin <= 2
    options.Display = 'final';
end
this.clearStates;

NN = this.NumberNodes;
NS = this.NumberSlices;
%% count the usage of links and nodes
% Here, the simplest method is adopted. If a link/node might be used by a slice, then its
% count is increased by the slice's weight.
% another method: if a link/node might be used by a flow in a slice, then its count is
% increased by the slice's weight.
link_usage_number = zeros(this.NumberLinks, NS);
node_usage_number = zeros(NN, NS);
for s = 1:NS
    link_usage_number(this.slices{s}.VirtualLinks.PhysicalLink,s) = slice_weight(s);
    node_usage_number(this.slices{s}.VirtualNodes.PhysicalNode,s) = slice_weight(s);
end
aggr_link_usage_number = sum(link_usage_number, 2);
aggr_node_usage_number = sum(node_usage_number, 2);

%% partition the network resource
% the partition ratio is determined by the usage number
for s = 1:NS
    link_id = this.slices{s}.VirtualLinks.PhysicalLink;
    partition_ratio = link_usage_number(link_id,s) ./ aggr_link_usage_number(link_id);
    this.slices{s}.VirtualLinks.Capacity = ...
        partition_ratio.*this.getLinkField('Capacity', link_id);
    node_id = this.slices{s}.VirtualNodes.PhysicalNode;
    partition_ratio = node_usage_number(node_id,s) ./ aggr_node_usage_number(node_id);
    this.slices{s}.VirtualNodes.Capacity = ...
        partition_ratio.*this.getNodeField('Capacity', node_id);
end

%% Independently optimize each network slice
utility = zeros(NS,1);
node_load = zeros(NN, NS);
link_load = zeros(this.NumberLinks, NS);
for s = 1:NS
    utility(s) = this.slices{s}.optimalFlowRate(options);
    node_load(this.slices{s}.VirtualNodes.PhysicalNode, s) = ...
        this.slices{s}.VirtualNodes.Load;
    link_load(this.slices{s}.VirtualLinks.PhysicalLink, s) = ...
        this.slices{s}.VirtualLinks.Load;   
end
aggr_node_load = sum(node_load, 2);
aggr_link_load = sum(link_load, 2);

%% Finalize substrate network
link_uc = this.getLinkField('UnitCost');
node_uc = this.getNodeField('UnitCost');
epsilon = this.unitStaticNodeCost;
phis_n = epsilon*this.delta*(NN-1)/this.totalNodeCapacity;
phis_l = epsilon*(1-this.delta)*(NN-1)/this.totalLinkCapacity;
link_price = (link_uc + phis_l) * (1 + options.PricingFactor);
node_price = (node_uc + phis_n) * (1 + options.PricingFactor);
for s = 1:NS
    sl = this.slices{s};
    sl.VirtualLinks.Price = link_price(sl.VirtualLinks.PhysicalLink);
    sl.VirtualNodes.Price = node_price(sl.VirtualNodes.PhysicalNode);
end
this.setNodeField('Load', aggr_node_load);
this.setLinkField('Load', aggr_link_load);
this.setLinkField('Price', link_price);
this.setNodeField('Price', node_price);

%% calculate the output (net social welfare)
output.node_price = node_price;
output.link_price = link_price;
output.node_load = aggr_node_load;
output.link_load = aggr_link_load;
output.welfare_approx = sum(utility);
output.welfare_accurate = 0;
t = zeros(NS+1, 1);
output.profit = table(t, t, t, t, 'VariableNames',...
    {'ApproximatePercent', 'ApproximatePrice', 'AccuratePercent', 'AccuratePrice'});
clear t;
output.flow_rate = [];
options.Model = 'Accurate';
for s = 1:NS
    sl = this.slices{s};
    nid = sl.VirtualNodes.PhysicalNode;
    output.flow_rate = [output.flow_rate; sl.FlowTable.Rate];

    var_x = [sl.Variables.x; sl.Variables.z];
    p = -Slice.fcnNetProfit(var_x, sl);
    output.profit.ApproximatePercent(s) = options.PercentFactor * p;
    
    p = -Slice.fcnNetProfit(var_x, sl, options);
    idx = sl.VirtualNodes.Load>0;
    p = p - dot(sl.VirtualNodes.Load(idx)./this.getNodeField('Load',(nid(idx))),...
        this.getNodeField('StaticCost', nid(idx)));
    output.profit.AccuratePercent(s) = options.PercentFactor * p;
    
    output.profit.ApproximatePrice(s) = -Slice.fcnProfit(var_x, sl);
    
    output.welfare_accurate = output.welfare_accurate...
        + sl.weight*sum(log(sl.FlowTable.Rate));
end
f = (1-options.PercentFactor)/options.PercentFactor;
output.profit.ApproximatePercent(end) = ...
    sum(output.profit.ApproximatePercent(1:(end-1))*f);
output.profit.AccuratePercent(end) = ...
    sum(output.profit.AccuratePercent(1:(end-1))*f);
output.profit.ApproximatePrice(end) = ...
    output.welfare_approx - sum(output.profit.ApproximatePrice(1:(end-1)));
output.profit.AccuratePrice = output.profit.ApproximatePrice;
output.profit.AccuratePrice(end) = ...
    output.welfare_accurate - sum(output.profit.AccuratePrice(1:(end-1)));
output.welfare_accurate = output.welfare_accurate - ...
    this.getNetworkCost([],[], options.Model);
if strncmp(options.Display,'iter', 4) || ...
        strncmp(options.Display,'notify', 6) || strncmp(options.Display,'final', 5)
    fprintf('\tOptimal net social welfare: fx = %G.\n', output.welfare_approx);
end
end