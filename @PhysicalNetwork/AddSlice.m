%% Add a Slice to Substrate Network
% Add a slice to physical network. See also <RemoveSlice>

%% Method Prototype
%
%     AddSlice(phy_net, slice_opt)
%
% |slice_opt|: Options on how to creat a slice.
%
% * *Type*:
% * *Weight*:
% * *FlowPattern*: flow pattern of the slice, see alse _FlowPattern_.
% * *NumberFlows*: this field is specified when |Pattern = 'RandomMultiFlow'|.
% * *NumberPaths*: number of candidate paths for each service flow.
% * *NumberVNFs*:
% * *RandomSeed*:
% * *Identifier*: the identifer of this slice. This is different from the index of the
% slice, and this field is not required to set.
% * *DelayConstraint*:
% * *ConstantProfit*:
function sl = AddSlice(this, slice_opt, varargin)
%% when arguments is not given, provide default value for slice data.
% if nargin < 2 || isempty(slice_opt)
%     slice_opt.Weight = 1;
%     slice_opt.FlowPattern = FlowPattern.RandomSingleFlow;
%     slice_opt.NumberPaths = 3;
%     slice_opt.VNFList = unique_randi(this.NumberVNFs, 3);
%     slice_opt.DelayConstraint = inf;
% end
slice_opt = this.preAddingSlice(slice_opt);

%% create flow table
graph = this.residualgraph(slice_opt);
[slice_opt.FlowTable, phy_adjacent, flag] = this.generateFlowTable(graph, slice_opt);
sl = Slice.empty(0,1);
if flag == false
    return;
end
number_flow = height(slice_opt.FlowTable);
if number_flow == 0
    cprintf('comments', 'Info: all flow are rejected.\n');
    return;
end
slice_opt.FlowTable.Properties.VariableNames = ...
    {'Source', 'Target', 'Rate', 'Delay', 'Paths'};

phy_node_id = find(transpose(sum(phy_adjacent,1)~=0) | sum(phy_adjacent,2)~=0);
slice_opt.Adjacent = phy_adjacent(phy_node_id, phy_node_id);
slice_opt.NumberNodes = length(slice_opt.Adjacent);
slice_opt.NumberEdges = nnz(slice_opt.Adjacent);
slice_opt.NodeMapS2P = phy_node_id;
% slice_opt.node_price = ones(length(phy_node_id),1);
slice_opt.NodeMapP2S = zeros(this.NumberNodes,1);
slice_opt.NodeMapP2S(phy_node_id) = (1:slice_opt.NumberNodes)';
slice_opt.LinkMapS2P = zeros(slice_opt.NumberEdges,1);
slice_opt.LinkMapP2S = zeros(this.NumberLinks,1);
[head, tail] = find(slice_opt.Adjacent~=0);
%%
% |head| and |tail| is the index of node in slice, <this.graph.IndexEdge> require the physical
% ID of the node.
slice_opt.LinkMapS2P = this.graph.IndexEdge(phy_node_id(head),phy_node_id(tail));
slice_opt.LinkMapP2S(slice_opt.LinkMapS2P) = 1:length(head);
slice_opt.FlowTable.Source = slice_opt.NodeMapP2S(slice_opt.FlowTable.Source);
slice_opt.FlowTable.Target = slice_opt.NodeMapP2S(slice_opt.FlowTable.Target);
slice_opt.Parent = this;

% To perform the slice admitting control, we should first add the new slice into the
% network, and perform resource allocation.
%% TODO
% change slice's storage to ListArray
sl = this.createslice(slice_opt, varargin{:});
%%%
% Assign identifier to flow/path of the slice.
global DEBUG; %#ok<NUSED>
sl.FlowTable.Identifier = this.flow_identifier_generator.next(sl.NumberFlows);
this.AllocatePathId(sl);

%% resource usage
link_id = sl.Links.PhysicalLink;
this.AggregateLinkUsage(link_id) = this.AggregateLinkUsage(link_id) + 1;
node_id = sl.Nodes.PhysicalNode;
this.AggregateNodeUsage(node_id) = this.AggregateNodeUsage(node_id) + 1;
end