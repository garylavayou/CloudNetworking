%% Static Network Slicing

%%
function output = StaticSlicing(this, slice, options)
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

if nargin>1 && ~isempty(slice)
    %% Allocate Resource to the new arrival slice
    % The residual capacity of the substrate network is available to the slice.
    ss = Slice(slice);
    ss.VirtualNodes.Capacity = ...
        this.getNodeField('ResidualCapacity', ss.VirtualNodes.PhysicalNode);
    ss.VirtualLinks.Capacity = ...
        this.getLinkField('ResidualCapacity', ss.VirtualLinks.PhysicalLink);
    net_profit = ss.optimalFlowRate(options);  
    slice.Variables = ss.Variables;
    slice.setPathBandwidth;
    slice.FlowTable.Rate = ss.FlowTable.Rate;
    slice.VirtualLinks.Load = ss.VirtualLinks.Load;
    slice.VirtualNodes.Load = ss.VirtualNodes.Load;
end

%% Finalize substrate network
% # After the optimization, the flow rate, node/link load have been recorded.
% # Announce the resource prices to the new slice.
% # Record the substrate network's node and link load.
NS = this.NumberSlices;
NN = this.NumberNodes;
link_uc = this.getLinkField('UnitCost');
node_uc = this.getNodeField('UnitCost');
epsilon = this.unitStaticNodeCost;
phis_n = epsilon*this.delta*(NN-1)/this.totalNodeCapacity;
phis_l = epsilon*(1-this.delta)*(NN-1)/this.totalLinkCapacity;
link_price = (link_uc + phis_l) * (1 + options.PricingFactor);
node_price = (node_uc + phis_n) * (1 + options.PricingFactor);
node_load = zeros(NN, NS);
link_load = zeros(this.NumberLinks, NS);
sl = this.slices{NS};
sl.VirtualLinks.Price = link_price(sl.VirtualLinks.PhysicalLink);
sl.VirtualNodes.Price = node_price(sl.VirtualNodes.PhysicalNode);
for s = 1:NS
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
% # Price of physical nodes and links;
% # Load of physical nodes and links;
% # The approximation of net social welfare, and the accurate net social welfare;
% # The net profit of each slice and the substrate network, including four results from
% different methods, i.e. |ApproximatePercent|, |ApproximatePrice|, |AccuratePercent|,
% |AccuratePrice|. 
% # flow rate of all flows in the network.
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

end

