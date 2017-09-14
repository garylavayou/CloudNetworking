%% Compact form of Hessian matrix of the Largrangian
function hess = fcnHessianCompact(act_var_x, ~, S)
var_x = zeros(S.num_vars,1);
var_x(S.I_active_variable) = act_var_x;
hess = Slice.fcnHessian(var_x, [], S);
hess = hess(S.I_active_variable, S.I_active_variable);
end