%% Utilization ratio (temp)
disp('Node utilization:');
disp(sum(this.getNodeLoad(this.temp_vars.z))/...
    sum(sum(reshape(this.temp_vars.v, this.NumberDataCenters, this.NumberVNFs),2)));

disp('Link utilization:');
disp(sum(this.getLinkLoad(this.temp_vars.x))/sum(this.temp_vars.c));
%% Utilization ratio (final)
disp('Node utilization:');
disp(sum(this.VirtualDataCenters.Load)/sum(this.VirtualDataCenters.Capacity));
disp('Link utilization:');
disp(sum(this.VirtualLinks.Load)/sum(this.VirtualLinks.Capacity));
disp('Overal utilization:');
disp(this.utilizationRatio());
%%
idx = find(results.Dimconfig{1}.ReconfigType==2);
idx = idx(idx>1);
r1 = mean(results.Dimconfig{1}.NumberReconfigVariables(idx));
%
idx = find(results.DimconfigReserve{1}.ReconfigType==2);
idx = idx(idx>1);
r2 = mean(results.DimconfigReserve{1}.NumberReconfigVariables(idx));
disp([r1,r2]);
%%
idx = results.Dimconfig{1}.ReconfigType==1;
r1 = mean(results.Dimconfig{1}.NumberReconfigVariables(idx));
%
idx = results.DimconfigReserve{1}.ReconfigType==1;
r2 = mean(results.DimconfigReserve{1}.NumberReconfigVariables(idx));
disp([r1,r2]);
%%
idx = results.Dimconfig{1}.ReconfigType==1;
r1 = mean(results.Dimconfig{1}.NumberReconfigFlows(idx));
%
idx = results.DimconfigReserve{1}.ReconfigType==1;
r2 = mean(results.DimconfigReserve{1}.NumberReconfigFlows(idx));
disp([r1,r2]);
%%
disp([mean(results.Dimconfig{1}.NumberReconfigVariables(200:end)), ...
mean(results.DimconfigReserve{1}.NumberReconfigVariables(200:end))]);

%%
disp([mean(results.Dimconfig{1}.NumberReconfigFlows(100:end)), ...
mean(results.DimconfigReserve{1}.NumberReconfigFlows(100:end))]);
%%
disp([sum(results.DimconfigReserve{1}.NumberReconfigFlows) ...
    sum(results.DimconfigReserve{1}.NumberReconfigVariables)]);
%   2775       12507    update reconfig prices each dim iteration
%   3061       12779    shorten dim reconfig cost distribute interval (1/3)
%   3028       12659    shorten dim reconfig cost distribute interval (1/4)
%   2884       11364    hsr + fsr (1/4)
disp([sum(results.Dimconfig{1}.NumberReconfigFlows) ...
    sum(results.Dimconfig{1}.NumberReconfigVariables)]);
%   3329       14926    shorten dim reconfig cost distribute interval (1/3)
%   3221       14593    shorten dim reconfig cost distribute interval (1/4)
%   3366       14289
%%
m = zeros(length(results.DimBaseline),1);
for i = 1:length(results.DimBaseline)
   m(i) = mean(results.DimBaseline{i}{:,'NumberReconfigVariables'}./results.DimBaseline{i}{:,'NumberVariables'});
end
disp(m)
%% 
idx = 3;
1-sum(results.Dimconfig{idx}{:, 'NumberReconfigVariables'})/sum(results.DimBaseline{idx}{:, 'NumberReconfigVariables'})
1-sum(results.DimconfigReserve{idx}{:, 'NumberReconfigVariables'})/sum(results.Dimconfig{idx}{:, 'NumberReconfigVariables'})
1-sum(results.Dimconfig{idx}{:, 'NumberReconfigFlows'})/sum(results.DimBaseline{idx}{:, 'NumberReconfigFlows'})
1-sum(results.DimconfigReserve{idx}{:, 'NumberReconfigFlows'})/sum(results.DimBaseline{idx}{:, 'NumberReconfigFlows'})
