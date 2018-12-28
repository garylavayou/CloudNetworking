function [fval,  grad] = subproblemObjective(var_x, lambda, slice)
%SUBPROBLEMOBJECTIVE the subproblem's objective function, used in the dual decomposition
%   method. This subproblem is associate with one of the slice.
var_path = var_x(1:slice.NumberPaths);
var_node = var_x((slice.NumberPaths+1):end);
link_load = slice.getLinkLoad(false, var_path);    % the load of virtual links in the slice, equal to <getLinkCapacity>
node_load = slice.getNodeLoad(false, var_node);    % the load of virtual nodes in the slice 

%% net profit of a slice
% The net profit of a slice (objective value) is equal to the utility less the slice cost
% and the cost introduced by dual variables.
%
% *Note*: cannot use _CloudNetworkEx.getNetworkCost_ to calculate the slice cost.
profit = slice.Weight*sum(fcnUtility(slice.getFlowRate(var_path))) - ...
    slice.getResourceCost(node_load, link_load);

%% dual variable cost
dual_cost = dot(lambda.n, node_load) + dot(lambda.e, link_load);
% if isfield(lambda, 'p')
%     dual_cost = dual_cost - dot(lambda.p, var_path);
% end
% if isfield(lambda, 'npf')
%     dual_cost = dual_cost - dot(lambda.npf(:),var_node);
% end
% |S.As_res*lambda.pf(:)| is a compatible arithmetic operation.

fval = -profit + dual_cost;

if nargout <= 1
    fval = -fval;
else
    %% gradient of objective function
    % The upper bound number of non-zero elements in the gradient vector:
    %  the gradient on path variable is nonzeros, so there is |P| components;
    %  whether the gradient on node variable is zeros is depend on the node-path
    %  incidence matrix, so the number of non-zero elements is less than i.e.
    %  |nnz(I_dc_path)*F|.
    link_uc = slice.Parent.readLink('UnitCost', slice.VirtualLinks.PhysicalLink);
    node_uc = slice.Parent.readNode('UnitCost', slice.VirtualNodes.PhysicalNode);
    grad = spalloc(length(var_x),1, slice.NumberPaths+nnz(slice.I_dc_path)*slice.NumberVNFs);
    for p = 1:slice.NumberPaths
        i = slice.path_owner(p);
        grad(p) = -slice.Weight/(slice.I_flow_path(i,:)*var_path) + ...
            dot((link_uc+slice.Parent.phis_l+lambda.e),slice.I_edge_path(:,p)); %#ok<SPRIX>
    end
    
    nz = slice.NumberVirtualNodes*slice.NumberPaths;
    z_index = slice.NumberPaths+(1:nz);
    for f = 1:slice.NumberVNFs
        %% Each iteration: find z(:,:,f)'s derivatives.
        % compatible arithmetic operation: node_unit_cost, s_n and lambda.n is column
        % vectors and (lambda.pf(:,f))' is a row vectors, the result is a matrix.
        grad(z_index) = (node_uc + slice.Parent.phis_n + lambda.n).*slice.I_dc_path; %#ok<SPRIX>
        z_index = z_index + nz;
    end
end
end