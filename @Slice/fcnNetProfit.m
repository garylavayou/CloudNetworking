%% Net profit of the slice
% The net profit is the total utility of the slice less the resources cost.
% The constant part of the objective fucntion can be ignored.
% 
% When calculate the cost of the network as a single slice, _Slice.getSliceCost_ is
% equivalent to _PhysicalNetwork.getNetworkCost_. Otherwise, use _Slice.getSliceCost_ to
% calculate each slice's cost.
%
% Called by _singleSliceOptimization_ as objective function;
% Called by _optimalFlowRate_ as objective function;
function [profit, grad] = fcnNetProfit(var_x, S, options)
var_path = var_x(1:S.NumberPaths);
var_node = var_x((S.NumberPaths+1):end);
link_load = S.getLinkLoad(var_path);
node_load = S.getNodeLoad(var_node);
if isempty(S.weight)        % for network as a single slice.
    weight = S.FlowTable.Weight;
    profit = -sum(weight.*log(S.getFlowRate(var_path))) + ...
        S.Parent.getNetworkCost(node_load, link_load);
else                        % for network as multiple slice.
    weight = S.weight*ones(S.NumberFlows, 1);
    if nargin>=3 && isfield(options, 'Model')
        profit = -sum(weight.*log(S.getFlowRate(var_path))) + ...
            S.getSliceCost(node_load, link_load, options.Model);
    else
        profit = -sum(weight.*log(S.getFlowRate(var_path))) + ...
            S.getSliceCost(node_load, link_load);
    end
end

if nargout >=2
    link_uc = S.Parent.getLinkField('UnitCost', S.VirtualLinks.PhysicalLink);
    node_uc = S.Parent.getNodeField('UnitCost', S.VirtualNodes.PhysicalNode);
    delta = S.Parent.delta;
    epsilon = S.Parent.unitStaticNodeCost;
    N = S.Parent.NumberNodes;
    phis_n = epsilon*delta*(N-1)/S.Parent.totalNodeCapacity;
    phis_l = epsilon*(1-delta)*(N-1)/S.Parent.totalLinkCapacity;
    grad = spalloc(length(var_x),1, S.NumberPaths+nnz(S.I_node_path)*S.NumberVNFs);
    for p = 1:S.NumberPaths
        i = S.path_owner(p);
        f = S.I_flow_path(:,p)~=0;  % find the weight of the flow on the path.
        grad(p) = -weight(f)/(S.I_flow_path(i,:)*var_path) + ...
            (link_uc + phis_l)'*S.I_edge_path(:,p); %#ok<SPRIX>
    end
    
    nz = S.NumberVirtualNodes*S.NumberPaths;
    z_index = S.NumberPaths+(1:nz);
    for f = 1:S.NumberVNFs
        % compatiable arithmetic operation
        grad(z_index) = (node_uc + phis_n).*S.I_node_path; %#ok<SPRIX>
        z_index = z_index + nz;
    end
end
end