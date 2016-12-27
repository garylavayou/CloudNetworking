%% Hessian matrix of the Largrangian
% since the problem only contains linear constraint, the hessian matrix of the
% Largrangian is equal to the seconderivatives of the objective function, and the
% Largrangian multipliers $\lambda$ takes no effect.
%
% the Hessian matrix contains only $P^2$ nonzeros elements on the diagonal,
% which is the second derviatives on path variables.
% Called by _singleSliceOptimization_;

function hess = fcnHessian(var_x, ~, S)
if isempty(S.weight)
    weight = S.FlowTable.Weight;    % for single slice;
else
    weight = S.weight*ones(S.NumberFlows, 1);  % for multiple slices
end
hess = spalloc(length(var_x),length(var_x), S.NumberPaths^2);
var_path = var_x(1:S.NumberPaths);
for p = 1:S.NumberPaths
    i = S.path_owner(p);
    f = S.I_flow_path(:,p)~=0;
    hess(p,1:S.NumberPaths) = weight(f)*...
        S.I_flow_path(i,:)/((S.I_flow_path(i,:)*var_path))^2; %#ok<SPRIX>
end
end
