%%
plot(1:NUM_EVENT, results.Reconfig.Profit, ...
    1:NUM_EVENT, results.Fastconfig.Profit);
%%
plot(1:NUM_EVENT, results.Reconfig.Cost, ...
    1:NUM_EVENT, results.Fastconfig.Cost);
legend('Reconfig', 'Fastconfig');
ylabel('Reconfiguration Cost (Const)');
%%
plot(1:NUM_EVENT, results.Reconfig.LinearCost, ...
    1:NUM_EVENT, results.Fastconfig.LinearCost);
legend('Reconfig', 'Fastconfig');
ylabel('Reconfiguration Cost (Linear)');
%%
plot(1:NUM_EVENT, results.Reconfig.NumberReconfigVariables, ...
    1:NUM_EVENT, results.Fastconfig.NumberReconfigVariables);
legend('Reconfig', 'Fastconfig');
ylabel('Number of Reconfigured Variables');
%%
plot(1:NUM_EVENT, results.Reconfig.RatioReconfigVariables, ...
    1:NUM_EVENT, results.Fastconfig.RatioReconfigVariables);
legend('Reconfig', 'Fastconfig');
ylabel('Ratio of Reconfigured Variables');
%%
% plot(1:NUM_EVENT, results.Reconfig.NumberReconfigFlows, ...
%     1:NUM_EVENT, results.Fastconfig.NumberReconfigFlows);
% legend('Reconfig', 'Fastconfig');
% ylabel('Numebr of Reconfigured Variables');
