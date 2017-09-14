function [tb, stbs] = saveStatTable(PN, output, rt, slice_types, method)
tb = table;
% tb{1, 'Runtime'} = max(rt(rt~=0));
switch method 
    case 'static'
        tb{1, 'Runtime'} = max(rt(rt~=0));
        tb{1, 'NumberFlowsStatic'} = PN.NumberFlows;  
        tb{1, 'NumberSlicesStatic'} = PN.CountSlices';
    case 'optimal-spp'
        tb{1, 'ApproximateWelfare'} = ...
            [output.welfare_approx_optimal, output.welfare_approx];
        tb{1, 'AccurateWelfare'} = ...
            [output.welfare_accurate_optimal, output.welfare_accurate];
    case 'dynamic-price'
    otherwise
        error('error: invalid type (%s)', method);
end
if cellstrfind({'optimal-spp', 'dynamic-price'}, method)
    tb{1, 'Runtime'} = [rt.Parallel, rt.Serial];
    tb{1, 'NumberFlows'} = PN.NumberFlows;
    tb{1, 'NumberSlices'} = PN.CountSlices';
end
if cellstrfind({'static', 'dynamic-price'}, method)
    tb{1, 'ApproximateWelfare'} = output.welfare_approx;
    tb{1, 'AccurateWelfare'} = output.welfare_accurate;
end
tb{1, 'Utilization'} = PN.utilizationRatio;
[r1, r2, r3, r4] = PN.nodeUtilization;
tb{1, 'NodeUtilization'} = [r1, r2, r3, r4];
[r1, r2, r3, r4] = PN.linkUtilization;
tb{1, 'LinkUtilization'} = [r1, r2, r3, r4];
tb{1, 'Cost'} = PN.getNetworkCost;
tb{1, 'Profit'} = output.profit.ApproximatePrice(end);

stbs = table;
for j = 1:length(slice_types)
    [p,r] = PN.statSlice(slice_types(j), output.profit.ApproximatePrice);
    stbs(j, {'Profit','Rate'}) = {p,r};
end
end

