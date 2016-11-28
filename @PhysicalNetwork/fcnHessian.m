%% Hessian matrix of the Largrangian
% See also Slice.fcnHessian
function hess = fcnHessian(var_x, ~, S)
weight = S.FlowTable.Weight;
hess = spalloc(length(var_x),length(var_x), S.NumberPaths^2);
var_path = var_x(1:S.NumberPaths);
for p = 1:S.NumberPaths
    i = S.path_owner(p);
    f = S.I_flow_path(:,p)~=0;
    hess(p,1:S.NumberPaths) = weight(f)*...
        S.I_flow_path(i,:)/((S.I_flow_path(i,:)*var_path))^2; %#ok<SPRIX>
end
end

