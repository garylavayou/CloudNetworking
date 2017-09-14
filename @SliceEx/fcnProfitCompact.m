%% Evaluate the objective function and gradient
% only active independent variables are passed into the objective function.
% Considering the constraint's coefficient matrix, if the corresponding column of the
% coefficient matrix for a variable is all zero, then this variable is inactive and can be
% directly set as 0. So we can remove it from the optimization problem.
%
% NOTE: we can also remove the all-zero rows of the coefficient matrix, which do not
% influence the number of variables. See also <optimalFlowRate>.
function [profit, grad]= fcnProfitCompact(act_var_x, S)
var_x = zeros(S.num_vars,1);
var_x(S.I_active_variable) = act_var_x;

% we extend the active variables by adding zeros to the inactive ones.
[profit, grad] = Slice.fcnProfit(var_x, S);

% eliminate the inactive variable's derivatives.
grad = grad(S.I_active_variable);
end