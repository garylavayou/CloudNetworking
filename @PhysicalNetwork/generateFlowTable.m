function [flow_table, phy_adjacent, flag] = generateFlowTable( this, graph, slice_opt )
% generateFlowTable Generate random flows with random selected node pairs.
% called by _AddSlice_ and _creatflow_.
flow_table = table;
if nargout >= 3
    flag = true;
end
if nargout >= 2
    phy_adjacent = spalloc(this.NumberNodes, this.NumberNodes, this.NumberLinks);
end
options = this.updatePathConstraints(slice_opt);
info = this.updateDemandInfo(slice_opt);
if isfield(slice_opt, 'DuplicateFlow') && slice_opt.DuplicateFlow==true
    b_duplicate = true;
    number_flow = slice_opt.NumberFlows;
else
    b_duplicate = false;
    number_flow = min(slice_opt.NumberFlows, info.NumberFlows);
end

k = 0;
while k < number_flow
    %% TODO the source and targets can be the same node
    % In this case, the flow only consumes computing resources on data centers.
    end_points = this.generateEndPoints(info, slice_opt);
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
            error('PhysicalNetwork:Disconnected', ...
                'error: cannot find feasible path between %d and %d.', ...
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

