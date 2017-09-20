function [flow_table, phy_adjacent, flag] = generateFlowTable( this, graph, slice_opt )
%generateFlowTable Generate with specified flow pattern.
%      Flow pattern are listed in <FlowPattern>.
switch slice_opt.FlowPattern
    case {FlowPattern.RandomSingleFlow, FlowPattern.RandomMultiFlow}
        if nargout == 1
            flow_table = ...
                generateFlowTable@CloudNetwork( this, graph, slice_opt );
        else
            [flow_table, phy_adjacent, flag] = ...
                generateFlowTable@CloudNetwork( this, graph, slice_opt );
        end
        return;
    case FlowPattern.RandomInterDataCenter
        number_flow = this.NumberDataCenters*(this.NumberDataCenters-1);
        if isfield(slice_opt, 'NodeSet')
            node_set = slice_opt.NodeSet;
        else
            node_set = this.DataCenters.NodeIndex;
        end
        numnodes = length(node_set);
    case FlowPattern.RandomInterBaseStation
        number_flow = this.NumberBaseStations*(this.NumberBaseStations-1);
        if isfield(slice_opt, 'NodeSet')
            node_set = slice_opt.NodeSet;
        else
            node_set = this.base_stations;
        end
        numnodes = length(node_set);
    case FlowPattern.RandomDataCenter2BaseStation
        number_flow = this.NumberBaseStations*this.NumberDataCenters;
        if isfield(slice_opt, 'BSNodeSet')
            bs_set = slice_opt.BSNodeSet;
        else
            bs_set = this.base_stations;
        end
        if isfield(slice_opt, 'DCNodeSet')
            dc_set = slice_opt.DCNodeSet;
        else
            dc_set = this.DataCenters.NodeIndex;
        end
        numdc = length(dc_set);
    otherwise
        error('error: slice type cannot be loaded.');
end
if isfield(slice_opt, 'DuplicateFlow') && slice_opt.DuplicateFlow==true
    b_duplicate = true;
    number_flow = slice_opt.NumberFlows;
else
    b_duplicate = false;
    number_flow = min(slice_opt.NumberFlows, number_flow);
end
flow_table = table;
if nargout >= 3
    flag = true;
end
if nargout >=2
    phy_adjacent = spalloc(this.NumberNodes, this.NumberNodes, this.NumberLinks);
end
k = 0;
options = this.updatePathConstraints(slice_opt);

while k < number_flow
    %% TODO the source and targets can be the same node
    % In this case, the flow only consumes computing resources on data centers.
    switch slice_opt.FlowPattern
        case {FlowPattern.RandomInterDataCenter,FlowPattern.RandomInterBaseStation}
            id = unique_randi(numnodes, 2, 'stable');
            end_points = node_set(id);
        case FlowPattern.RandomDataCenter2BaseStation
            id = randi(numdc, 1);
            end_points(1) = dc_set(id);
            bs = bs_set;
            bs(bs==end_points(1)) = [];     % avoid the same node as BS and DC.
            end_points(2) = bs(randi(length(bs),1));
    end
    if height(flow_table)>0 && ~b_duplicate && ...
            ~isempty(find(flow_table{:,1}==end_points(1) & flow_table{:,2}==end_points(2),1)) 
        continue;
    end
    k = k+1;
    path_list = graph.CandidatePaths(slice_opt.NumberPaths, ...
        end_points(1), end_points(2), options);
    % assert path list
    switch this.assert_path_list(end_points, path_list, slice_opt)
        case -1  % failed with error;
            error('error: cannot find feasible path between %d and %d.',...
                end_points(1), end_points(2));
        case 0  % succeed
        case 1  % failed but continue
            continue;
        case 2  % failed and return;
            flag = false;
            return;
    end
    % Path list records the physical nodes on paths. Latterly, it will be transformed
    % into virtual nodes. See also <AddSlice> and <createflow>.
    flow_table(end+1,:) = ...
        {end_points(1),end_points(2),0,path_list.Latency,path_list}; %#ok<AGROW>
    if nargout >= 2
        for p = 1:path_list.Width
            path = path_list.paths{p};
            for l = 1:(path.Length-1)
                phy_adjacent(path.Node(l), path.Node(l+1)) = ...
                    graph.Adjacent(path.Node(l), path.Node(l+1));  %#ok<SPRIX>
            end
        end
    end
end

end

