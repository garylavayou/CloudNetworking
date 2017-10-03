function [stat, slice_stat] = createStatTable(num_point, num_type, type)
zero_stat0 = zeros(num_point, num_type);
zero_stat1 = zeros(num_point, 1);
zero_stat2 = zeros(num_point, 2);
zero_stat4 = zeros(num_point, 4);

stat = table;
if strcmpi(type, 'optimal-spp')
    stat.Welfare = zero_stat2;
%     stat.AccurateWelfare = zero_stat2;
    stat.Properties.VariableDescriptions{1} = 'Upperbound|ByPrice';
%     stat.ApproximateWelfare = zero_stat2;
%     stat.Properties.VariableDescriptions{2} = 'Upperbound|ByPrice';
else
    stat.Welfare = zero_stat1;
    %     stat.AccurateWelfare = zero_stat1;
    %     stat.ApproximateWelfare = zero_stat1;
end
stat.Profit = zero_stat1;
stat.Cost = zero_stat1;
stat.Utilization = zero_stat1;
stat.NodeUtilization = zero_stat4;
stat.Properties.VariableDescriptions{6} = 'Average|Max|Min|StdDev';
stat.LinkUtilization = zero_stat4;
stat.Properties.VariableDescriptions{7} = 'Average|Max|Min|StdDev';
slice_stat = cell(num_type,1);
for st = 1:num_type
    slice_stat{st} = table;
    slice_stat{st}.Rate = zero_stat4;
    slice_stat{st}.Properties.VariableDescriptions{1} = 'Average|Max|Min|StdDev';
    slice_stat{st}.Profit = zero_stat4;
    slice_stat{st}.Properties.VariableDescriptions{2} = 'Average|Max|Min|StdDev';
end

if strcmpi(type, 'static')
    stat.NumberFlowsStatic = zero_stat1;
    stat.NumberReject = zero_stat0;
    stat.NumberPartialReject = zero_stat0;
    stat.NumberSlicesStatic = zero_stat0;
else
    stat.NumberSlices = zero_stat0;
    stat.NumberFlows = zero_stat1;
end

if cellstrfind({'optimal-spp', 'dynamic-price'}, type)    
    stat.Runtime = zero_stat2;      % {Parallel|Serial}
else
    stat.Runtime = zero_stat1;
end