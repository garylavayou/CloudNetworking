%% Add a Slice to Substrate Network
% Add a slice to physical network.

%% Method Prototype
%
%     AddSlice(phy_net, slice_opt) 
%
% |slice_opt|: Options on how to creat a slice.
%
% * *Type*: 
% * *Weight*:
% * *Pattern*: flow pattern of the slice, see alse _FlowPattern_.
% * *NumberFlows*: this field is specified when |Pattern = 'RandomMultiFlow'|.
% * *NumberPaths*: number of candidate paths for each service flow.
% * *NumberVNFs*:
% * *RandomSeed*:
% * *Identifier*: the identifer of this slice. This is different from the index of the
% slice, and this field is not required to set.
% * *DelayConstraint*:
function sl = AddSlice(this, slice_opt)
% when arguments is not given, provide default value for slice data.
if nargin < 2 || isempty(slice_opt)
    slice_opt.Weight = 1;
    slice_opt.Pattern = FlowPattern.RandomSingleFlow;
    slice_opt.NumberPaths = 3;
    slice_opt.VNFList = unique_randi(this.NumberVNFs, 3);
    slice_opt.DelayConstraint = inf;
end

if ~isfield(slice_opt,'Weight') || isempty(slice_opt.Weight) || slice_opt.Weight == 0
    slice_opt.Weight = 1;
end
if ~isfield(slice_opt,'Pattern') || isempty(slice_opt.Pattern)
    slice_opt.Pattern = FlowPattern.RandomSingleFlow;
end
if ~isfield(slice_opt,'NumberPaths') || isempty(slice_opt.NumberPaths ) || ...
        slice_opt.NumberPaths == 0
    slice_opt.NumberPaths = 3;
end
if isfield(slice_opt, 'RandomSeed') && ~isempty(slice_opt.RandomSeed)
    rng(slice_opt.RandomSeed);
end
if ~isfield(slice_opt, 'VNFList') || isempty(slice_opt.VNFList)
    if ~isfield(slice_opt, 'NumberVNFs') || isempty(slice_opt.NumberVNFs)
        slice_opt.VNFList = unique_randi(this.NumberVNFs, randi([1 4]));
    else
        slice_opt.VNFList = unique_randi(this.NumberVNFs, slice_opt.NumberVNFs);
    end
else
    if isfield(slice_opt, 'NumberVNFs') && ~isempty(slice_opt.NumberVNFs)
        slice_opt.VNFList = slice_opt.VNFList(1:slice_opt.NumberVNFs);
    end
end
if ~isfield(slice_opt,'DelayConstraint') || isempty(slice_opt.DelayConstraint) || ...
        slice_opt.DelayConstraint == 0
    slice_opt.DelayConstraint = inf;
end
if slice_opt.Pattern ~= FlowPattern.RandomSingleFlow && ...
        ~isfield(slice_opt, 'NumberFlows')
    error('FlowPattern.RandomSingleFlow must be specified with NumberFlows.');
end
% create flow table
slice_opt.flow_table = table;
phy_adjacent = spalloc(this.NumberNodes, this.NumberNodes, this.NumberLinks);
switch slice_opt.Pattern
    case FlowPattern.RandomSingleFlow
        number_flow = 1;
    case FlowPattern.RandomMultiFlow
        number_flow = min(slice_opt.NumberFlows, this.NumberNodes*(this.NumberNodes-1));
    otherwise
        error('slice type cannot be loaded');
end
switch slice_opt.Pattern
    case {FlowPattern.RandomSingleFlow, FlowPattern.RandomMultiFlow}
        k = 0;
        while k < number_flow
            end_points = unique_randi(this.Topology.numnodes,2, 'stable');
            if height(slice_opt.flow_table)==0 ||...
                    isempty(find(slice_opt.flow_table{:,1}==end_points(1) & ...
                    slice_opt.flow_table{:,2}==end_points(2),1))
                k = k+1;
                path_list = this.graph.CandidatePaths(slice_opt.NumberPaths, ...
                    end_points(1), end_points(2), slice_opt.DelayConstraint, ...
                    this.LinkOption.delay);
                slice_opt.flow_table(end+1,:) = {end_points(1),end_points(2),0,0,path_list};
                for p = 1:path_list.Width
                    path = path_list.paths{p};
                    for l = 1:(path.Length-1)
                        phy_adjacent(path.Node(l), path.Node(l+1)) = ...
                            this.graph.Adjacent(path.Node(l), path.Node(l+1));  %#ok<SPRIX>
                    end
                end
            end
        end
    otherwise
        error('slice type cannot be loaded');
end
slice_opt.flow_table.Properties.VariableNames = ...
    {'Source', 'Target', 'Rate', 'Delay', 'Paths'};

phy_node_id = find(transpose(sum(phy_adjacent,1)~=0) | sum(phy_adjacent,2)~=0);
slice_opt.adjacent = phy_adjacent(phy_node_id, phy_node_id);
slice_opt.num_nodes = length(slice_opt.adjacent);
slice_opt.num_edges = nnz(slice_opt.adjacent);
slice_opt.node_map_S2P = phy_node_id;
slice_opt.node_price = ones(length(phy_node_id),1);
slice_opt.node_map_P2S = zeros(this.NumberNodes,1);
for i = 1:slice_opt.num_nodes
    slice_opt.node_map_P2S(slice_opt.node_map_S2P(i)) = i;
end

slice_opt.link_map_S2P = zeros(slice_opt.num_edges,1);
slice_opt.link_map_P2S = zeros(this.NumberLinks,1);
[head, tail] = find(slice_opt.adjacent~=0);
%%
% head and tail is the index of node in slice, this.graph.IndexEdge require the physical
% ID of the node
slice_opt.link_map_S2P = this.graph.IndexEdge(phy_node_id(head),phy_node_id(tail));
slice_opt.link_price = ones(length(head),1);
slice_opt.link_map_P2S(slice_opt.link_map_S2P) = 1:length(head);
slice_opt.flow_table.Source = slice_opt.node_map_P2S(slice_opt.flow_table.Source);
slice_opt.flow_table.Target = slice_opt.node_map_P2S(slice_opt.flow_table.Target);
for k = 1:number_flow
    path_list = slice_opt.flow_table{k,{'Paths'}};
    for p = 1:path_list.Width
        path = path_list.paths{p};
        path.node_list = slice_opt.node_map_P2S(path.node_list);
    end
end
slice_opt.parent = this;

this.NumberSlices = this.NumberSlices+1;
start_slice_id = this.NumberSlices;
this.slices{start_slice_id} = Slice(slice_opt);
sl = this.slices{start_slice_id};
this.AllocateFlowId(start_slice_id);
this.AllocatePathId(start_slice_id);

%% Update flow and path identifers
% After adding a slice, update the flow id and path id.
%
% TODO: assume the existing network slices has continuous id for flow and path
% thus this slice's flow and path id follow's the last existing network's flow and path.
% after remove a slice, all network slices' path and flow id should be reallocated.
% there is one exception, when the removed slice is the last one slice.
end
