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
% * *Pattern*: flow pattern of the slice, see alse _FlowPattern_.
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
%     slice_opt.Pattern = FlowPattern.RandomSingleFlow;
%     slice_opt.NumberPaths = 3;
%     slice_opt.VNFList = unique_randi(this.NumberVNFs, 3);
%     slice_opt.DelayConstraint = inf;
% end
if ~isfield(slice_opt,'Pattern') || isempty(slice_opt.Pattern)
    error('error: invalid flow pattern.'); %     slice_opt.Pattern = FlowPattern.RandomSingleFlow;
end
if ~isfield(slice_opt,'NumberPaths') || isempty(slice_opt.NumberPaths ) || ...
        slice_opt.NumberPaths == 0
    error('error: invalid number of paths.'); %     slice_opt.NumberPaths = 3;
end
if slice_opt.Pattern ~= FlowPattern.RandomSingleFlow && ...
        ~isfield(slice_opt, 'NumberFlows')
    error('FlowPattern.RandomSingleFlow must be specified with NumberFlows.');
end
if isfield(slice_opt, 'RandomSeed') && ~isempty(slice_opt.RandomSeed)
    rng(slice_opt.RandomSeed);
else
    warning('random number seed is not specified (set as %d).', floor(now));
    rng(floor(now));        % this value should be the same on the same day.
end
if ~isfield(slice_opt, 'VNFList') || isempty(slice_opt.VNFList)
    if ~isfield(slice_opt, 'NumberVNFs') || isempty(slice_opt.NumberVNFs) ||...
            slice_opt.NumberVNFs > this.NumberVNFs
        error('error: invalid VNF list options.'); % slice_opt.VNFList = unique_randi(this.NumberVNFs, randi([1 4]));
    else
        slice_opt.VNFList = unique_randi(this.NumberVNFs, slice_opt.NumberVNFs);
    end
else
    if max(slice_opt.VNFList) > this.NumberVNFs
        error('error: invalid VNF list options.');
    end
    % TODO check VNF list
    if isfield(slice_opt, 'NumberVNFs') && ~isempty(slice_opt.NumberVNFs)
        slice_opt.VNFList = slice_opt.VNFList(1:slice_opt.NumberVNFs);
%     else
%         slice_opt.VNFList = slice_opt.VNFList;
    end
end
if ~isfield(slice_opt,'DelayConstraint') || isempty(slice_opt.DelayConstraint) || ...
        slice_opt.DelayConstraint == 0
    slice_opt.DelayConstraint = inf;
end
if this.NumberDataCenters < this.NumberNodes
    % if only part of the forwarding nodes is VNF-capable, we should make sure that the
    % path at least transit one VNF-capable node.
    slice_opt.MiddleNodes = this.DataCenters.NodeIndex;
    % no matter when, the DataCenters is the middle nodes. However, if the MiddleNodes
    % option is not provided, the route calculation will be performed in a different way.
end

%% create flow table
graph = this.residualgraph(slice_opt);
[slice_opt.flow_table, phy_adjacent, flag] = this.generateFlowTable(graph, slice_opt);
sl = [];
if flag == false
    return;
end
number_flow = height(slice_opt.flow_table);
if number_flow == 0
    cprintf('comments', 'Info: all flow are rejected.\n');
    return;
end
slice_opt.flow_table.Properties.VariableNames = ...
    {'Source', 'Target', 'Rate', 'Delay', 'Paths'};

phy_node_id = find(transpose(sum(phy_adjacent,1)~=0) | sum(phy_adjacent,2)~=0);
slice_opt.adjacent = phy_adjacent(phy_node_id, phy_node_id);
slice_opt.num_nodes = length(slice_opt.adjacent);
slice_opt.num_edges = nnz(slice_opt.adjacent);
slice_opt.node_map_S2P = phy_node_id;
% slice_opt.node_price = ones(length(phy_node_id),1);
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
% slice_opt.link_price = ones(length(head),1);
slice_opt.link_map_P2S(slice_opt.link_map_S2P) = 1:length(head);
slice_opt.flow_table.Source = slice_opt.node_map_P2S(slice_opt.flow_table.Source);
slice_opt.flow_table.Target = slice_opt.node_map_P2S(slice_opt.flow_table.Target);
slice_opt.parent = this;
slice_opt = rmfield(slice_opt, 'method');

%% TODO
% To perform the slice admitting control, we should first add the new slice into the
% network, and perform resource allocation.
sl = this.createslice(slice_opt, varargin{:});     % Todo, change slice's storage to ListArray 
%%%
% Assign identifier to flow/path of the slice.
global DEBUG; %#ok<NUSED>
sl.FlowTable.Identifier = this.flow_identifier_generator.next(sl.NumberFlows);
this.allocatepathid(sl);

%% TODO
% used for resource partitioning and pricing.
slice_link_usage = slice_opt.link_map_P2S~=0;
this.link_usage = [this.link_usage slice_link_usage];
slice_node_usage = slice_opt.node_map_P2S~=0;
this.node_usage = [this.node_usage slice_node_usage];
end

