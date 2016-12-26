%% Evalute the objective function and gradient
% only active independent variables are passed into the objective function.
function [profit, grad]= fcnProfitCompact(act_var_x, S)
var_x = zeros(S.num_vars,1);
var_x(S.I_active_variable) = act_var_x;

[profit, grad] = Slice.fcnProfit(var_x, S);

grad = grad(S.I_active_variable);
end