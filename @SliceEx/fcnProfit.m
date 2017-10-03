%% Profit of the slice
% Evalute the function value and gradient of the objective.
%
% When prices are specified (explictly or implicitly), the profit is the total utility
% of the slice minus the total payment to the slice provider;
%%
function [profit, grad] = fcnProfit(var_x, S, options)
if nargout <= 1
    if nargin == 2
        profit = fcnProfit@Slice(var_x, S);
    else
        profit = fcnProfit@Slice(var_x, S, options);
    end
    profit = profit + S.constant_profit;
else
    if nargin == 2
        [profit, grad] = fcnProfit@Slice(var_x, S);
    else
        [profit, grad] = fcnProfit@Slice(var_x, S, options);
    end
    profit = profit - S.constant_profit;
end
end