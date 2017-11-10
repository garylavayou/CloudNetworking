%% Profit of the slice
% Evalute the function value and gradient of the objective.
%
% When prices are specified (explictly or implicitly), the profit is the total utility
% of the slice minus the total payment to the slice provider;
%%
function [profit, grad] = fcnProfit(vars, slice, options)
if nargout <= 1
    profit = fcnProfit@Slice(vars, slice, options);
    profit = profit + slice.constant_profit;
else
    [profit, grad] = fcnProfit@Slice(vars, slice, options);
    profit = profit - slice.constant_profit;
end
end