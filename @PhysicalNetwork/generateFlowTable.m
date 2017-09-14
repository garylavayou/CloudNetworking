function [flow_table, phy_adjacent, flag] = generateFlowTable( this, graph, slice_opt )
% generateFlowTable Generate random flows with random selected node pairs.
flow_table = table;
if nargout >= 2
    phy_adjacent = spalloc(this.NumberNodes, this.NumberNodes, this.NumberLinks);
end
switch slice_opt.Pattern
    case FlowPattern.RandomSingleFlow
        number_flow = 1;
    case FlowPattern.RandomMultiFlow
        number_flow = min(slice_opt.NumberFlows, this.NumberNodes*(this.NumberNodes-1));
    otherwise
        error('error: cannot handle the flow pattern <%s>.', ...
                        slice_opt.Pattern.char);
end

k = 0;
options = getstructfields(slice_opt, {'DelayConstraint', 'MiddleNodes'});
options.delay_opt = this.LinkOptions.delay;
if isfield(slice_opt, 'NodeSet')
    numnodes = length(slice_opt.NodeSet);
    nodeset = slice_opt.NodeSet;
else
    numnodes = this.Topology.numnodes;
    nodeset = 1:numnodes;
end

flag = true;
while k < number_flow
    end_points = nodeset(unique_randi(numnodes, 2, 'stable'));
    if height(flow_table)==0 || isempty(find(flow_table{:,1}==end_points(1) & ...
            flow_table{:,2}==end_points(2),1))
        k = k+1;
        path_list = graph.CandidatePaths(slice_opt.NumberPaths, ...
            end_points(1), end_points(2), options);
        % assert path list
        switch assert_path_list(end_points, path_list, slice_opt)
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
        flow_table(end+1,:) = {end_points(1),end_points(2),0,0,path_list}; %#ok<AGROW>
        if nargin >= 2
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

end

