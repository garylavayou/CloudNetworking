% single slice's objective function.
function [profit, grad] = fcnNetWelfare(var_x, S)
var_path = var_x(1:S.NumberPaths);
var_node = var_x((S.NumberPaths+1):end);
link_load = S.getLinkLoad(var_path);
node_load = S.getNodeLoad(var_node);
link_uc = S.Parent.getLinkField('UnitCost');
node_uc = S.Parent.getNodeField('UnitCost');
delta = S.Parent.Delta;
epsilon = S.Parent.staticNodeCost;
N = S.Parent.NumberNodes;
phis_n = epsilon*delta*(N-1)/S.Parent.totalNodeCapacity;
phis_l = epsilon*(1-delta)*(N-1)/S.Parent.totalLinkCapacity;
theta = S.Parent.networkUtilization(node_load, link_load);
weight = S.FlowTable.Weight;

%% Net social welfare of the network
% the constant part of the objective fucntion can be ignored.
profit = -sum(weight.*log(S.getFlowRate(var_path))) ...
    + dot(link_uc, link_load) + dot(node_uc, node_load) +...
    epsilon*((N-1)*theta+1);

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