%% Add a Slice to Substrate Network
% Generate Slice Data
%     AddSlice(phy_net, slice_opt) add a slice to phy_net use options in slice_opt.
function AddSlice(this, slice_opt)
if nargin < 2 || isempty(slice_opt)
    slice_opt.DelayConstraint = inf;
    slice_opt.NumberPaths = 3;
    slice_opt.Weight = 1;
    slice_data.VNFList = unique_randi(this.NumberVNFs, 3);
else % slice_opt is not empty
    if ~isfield(slice_opt,'DelayConstraint') || isempty(slice_opt.DelayConstraint) ||...
            slice_opt.DelayConstraint == 0
        slice_opt.DelayConstraint = inf;
    end
    if ~isfield(slice_opt,'NumberPaths') || isempty(slice_opt.NumberPaths ) || ...
            slice_opt.NumberPaths == 0
        slice_opt.NumberPaths = 3;
    end    
    if ~isfield(slice_opt,'Weight') || isempty(slice_opt.Weight) || slice_opt.Weight == 0
        slice_data.Weight = 1;
    else
        slice_data.Weight = slice_opt.Weight;
    end
    if isfield(slice_opt, 'RandomSeed') && ~isempty(slice_opt.RandomSeed)
        rng(slice_opt.RandomSeed);
    end
    if ~isfield(slice_opt, 'VNFList') || isempty(slice_opt.VNFList)
        % default service chain has 3 VNF.
        slice_data.VNFList = unique_randi(this.NumberVNFs, 3);
    else
        slice_data.VNFList = slice_opt.VNFList;
    end
end

% create flow table
slice_data.flow_table = table;
phy_adjacent = spalloc(this.NumberNodes, this.NumberNodes, this.NumberLinks);
switch slice_opt.Type
    case SliceType.RandomSingleFlow
        number_flow = 1;
    otherwise
        error('slice type cannot be loaded');
end
switch slice_opt.Type
    case {SliceType.RandomSingleFlow}
        k = 0;
        while k < number_flow
            end_points = unique_randi(this.Topology.numnodes,2, 'stable');
            if height(slice_data.flow_table)==0 ||...
                    isempty(slice_data.flow_table.Source==end_points(1) && ...
                    slice_data.flow_table.Target==end_points(2))
                k = k+1;
                path_list = this.graph.CandidatePaths(slice_opt.NumberPaths, ...
                    end_points(1), end_points(2), slice_opt.DelayConstraint, ...
                    this.LinkOption.delay);
                slice_data.flow_table(end+1,:) = {end_points(1),end_points(2),0,0,path_list};
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
slice_data.flow_table.Properties.VariableNames = ...
    {'Source', 'Target', 'Rate', 'Delay', 'Paths'};

phy_node_id = find(transpose(sum(phy_adjacent,1)~=0) | sum(phy_adjacent,2)~=0);
slice_data.adjacent = phy_adjacent(phy_node_id, phy_node_id);
slice_data.num_nodes = length(slice_data.adjacent);
slice_data.num_edges = nnz(slice_data.adjacent);
slice_data.node_map_S2P = phy_node_id;
slice_data.node_price = ones(length(phy_node_id),1);
slice_data.node_map_P2S = zeros(this.NumberNodes,1);
for i = 1:slice_data.num_nodes
    slice_data.node_map_P2S(slice_data.node_map_S2P(i)) = i;
end

slice_data.link_map_S2P = zeros(slice_data.num_edges,1);
slice_data.link_map_P2S = zeros(this.NumberLinks,1);
[head, tail] = find(slice_data.adjacent~=0);
%%
% head and tail is the index of node in slice, this.graph.IndexEdge require the physical
% ID of the node
slice_data.link_map_S2P = this.graph.IndexEdge(phy_node_id(head),phy_node_id(tail));
slice_data.link_price = ones(length(head),1);
slice_data.link_map_P2S(slice_data.link_map_S2P) = 1:length(head);
slice_data.flow_table.Source = slice_data.node_map_P2S(slice_data.flow_table.Source);
slice_data.flow_table.Target = slice_data.node_map_P2S(slice_data.flow_table.Target);
for k = 1:number_flow
    path_list = slice_data.flow_table{k,{'Paths'}};
    for p = 1:path_list.Width
        path = path_list.paths{p};
        path.node_list = slice_data.node_map_P2S(path.node_list);
    end
end


this.slices{this.NumberSlices+1} = Slice(slice_data);
this.NumberSlices = this.NumberSlices+1;

%% TODO after adding a slice, update the flow id and path id
% assume the existing network slices has continuous id for flow and path
% thus this slice's flow and path id follow's the last existing network's flow and path.
% after remove a slice, all network slices' path and flow id should be reallocated.
% there is one exception, when the removed slice is the last one slice.
end
