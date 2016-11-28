function [fval,  grad] = subproblemObjective(var_x, lambda, S)
%SUBPROBLEMOBJECTIVE the subproblem's objective function, used in the dual decomposition
%   method. This subproblem is associate with one of the slice.
var_path = var_x(1:S.NumberPaths);
var_node = var_x((S.NumberPaths+1):end);
link_load = S.getLinkLoad(var_path);    % the load of virtual links in the slice
node_load = S.getNodeLoad(var_node);    % the load of virtual nodes in the slice 
link_uc = S.Parent.getLinkField('UnitCost', S.VirtualLinks.PhysicalLink); % the virtual links's unit cost
node_uc = S.Parent.getNodeField('UnitCost', S.VirtualNodes.PhysicalNode); % the virtual nodes's unit cost
delta = S.Parent.Delta;
epsilon = S.Parent.staticNodeCost;
N = S.Parent.NumberNodes;
phis_n = epsilon*delta*(N-1)/S.Parent.totalNodeCapacity;
phis_l = epsilon*(1-delta)*(N-1)/S.Parent.totalLinkCapacity;
%% net social welfare in the subproblem
% when compute the static cost, the Capacity of all physical nodes and links is included,
% 
% the constant part of the objective fucntion can be ignored.
% 
profit = S.weight*sum(log(S.getFlowRate(var_path))) ...
    - dot(link_uc, link_load) - dot(node_uc, node_load) ...
    - (phis_n*sum(node_load)+phis_l*sum(link_load));

%% dual variable cost
% |S.As_res*lambda.pf(:)| is a compatible arithmetic operation.
dual_cost = dot(lambda.n, node_load) + dot(lambda.e, link_load);
% if isfield(lambda, 'p')
%     dual_cost = dual_cost - dot(lambda.p, var_path);
% end
% if isfield(lambda, 'npf')
%     dual_cost = dual_cost - dot(lambda.npf(:),var_node);
% end

% objective value
fval = -profit + dual_cost;

%% gradient of objective function
% 
grad = zeros(length(var_x),1);
for p = 1:S.NumberPaths
    i = S.path_owner(p);
    grad(p) = -S.weight/(S.I_flow_path(i,:)*var_path) + ...
        dot((link_uc+phis_l+lambda.e),S.I_edge_path(:,p));
end

nz = S.NumberVirtualNodes*S.NumberPaths;
z_index = S.NumberPaths+(1:nz);
for f = 1:S.NumberVNFs
    %% Each iteration: find z(:,:,f)'s derivatives.
    % compatible arithmetic operation: node_unit_cost, s_n and lambda.n is column vectors
    % and (lambda.pf(:,f))' is a row vectors, the result is a matrix.
    grad(z_index) = (node_uc + phis_n + lambda.n).*S.I_node_path;
    z_index = z_index + nz;
end

end