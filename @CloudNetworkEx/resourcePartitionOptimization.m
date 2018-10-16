%% Resource Partition Optimization
% Link and node resources are firstly evenly or proportionaly partiotioned into each
% slice, and then each slice independently run optimization. Each slice is given a fixed
% price, then each slice run optimization independently.
%
% The process is different from VNE, which knows the explicit demand of network resources,
% and try to find the embedding solution.
function [ output, runtime ] = resourcePartitionOptimization( this, slice_weight )
global DEBUG;
options = getstructfields(this.options, ...
    {'SlicingMethod', 'ProfitType', 'WelfareType', 'PercentFactor'});

this.clearStates;

NN = this.NumberNodes;
NS = this.NumberSlices;
%% count the usage of links and nodes
% Here, the simplest method is adopted. If a link/node might be used by a slice, then its
% count is increased by the slice's weight.
% Another method: if a link/node might be used by a flow in a slice, then its count is
% increased by the slice's weight. (NOTE: This two method do not make huge difference)
link_usage = zeros(this.NumberLinks, NS);
node_usage = zeros(NN, NS);
if nargin <= 1 || isempty(slice_weight)
%     slice_weight = ones(this.NumberSlices,1);
%     for s = 1:NS
%         sl = this.slices{s};
%         link_usage(sl.VirtualLinks.PhysicalLink,s) = slice_weight(s);
%         node_usage(sl.VirtualNodes.PhysicalNode,s) = slice_weight(s);
%     end
    for s = 1:NS
        sl = this.slices{s};
        % count the usage of edge with slice weight.
        % I_edge_path*I_flow_path' = I_edge_flow _> 
        I_edge_flow = sl.I_edge_path*sl.I_flow_path';
        I_node_flow = sl.I_dc_path*sl.I_flow_path';
        link_id = sl.VirtualLinks.PhysicalLink;
        node_id = sl.VirtualNodes.PhysicalNode;
        link_usage(link_id,s) = sum(I_edge_flow,2)*sl.weight;
        node_usage(node_id,s) = sum(I_node_flow,2)*sl.weight;
    end
else
    for s = 1:NS
        sl = this.slices{s};
        link_usage(sl.VirtualLinks.PhysicalLink,s) = slice_weight(s);
        node_usage(sl.VirtualNodes.PhysicalNode,s) = slice_weight(s);
    end
end

aggr_link_usage = sum(link_usage, 2);
aggr_node_usage = sum(node_usage, 2);

%% partition the network resource
% the partition ratio is determined by the usage number
for s = 1:NS
    link_id = this.slices{s}.VirtualLinks.PhysicalLink;
    partition_ratio = link_usage(link_id,s) ./ aggr_link_usage(link_id);
    this.slices{s}.VirtualLinks.Capacity = ...
        partition_ratio.*this.getLinkField('Capacity', link_id);
    node_id = this.slices{s}.VirtualNodes.PhysicalNode;
    partition_ratio = node_usage(node_id,s) ./ aggr_node_usage(node_id);
    this.slices{s}.VirtualNodes.Capacity = ...
        partition_ratio.*this.getNodeField('Capacity', node_id);
end

%% Independently optimize each network slice
if nargout == 2
    [prices, runtime] = pricingFactorAdjustment(this);
else
    [prices] = pricingFactorAdjustment(this);
end

% Finalize substrate network
this.finalize(prices);

%% calculate the output (net social welfare)
output = this.calculateOutput();
if ~isempty(DEBUG) && DEBUG
    fprintf('\tOptimal net social welfare: fx = %G.\n', output.welfare_approx);
end
end