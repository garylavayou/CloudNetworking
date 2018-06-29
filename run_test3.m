%% Physical Network Speification for SD_RAN
% This program used to test if the parameters are properly set.
% NOTE: link utilization ratio is lower, since not all links will be used for DC-2-DC
% traffic.

%% Specification of VNFs and Network Slices
% |StaticCostOption| is not set, the default value is |None|;
% |RandomSeed|: the first seed is for random static cost, the second is for process
% efficiency.
clear variables global;
link_opt.RandomSeed = 20170421;
link_opt.DelayModel = LinkDelayOption.Random;
link_opt.CostUnit = 0.2;
node_opt.CostUnit = 0.1;
net_opt.AdmitPolicy = 'reject-flow';
net_opt.PricingFactor = 2;      % 1 | 2 | 3 ,
net_opt.Threshold = 'average';

VNF_opt.Number = 6;            % number of VNF type
VNF_opt.Model = VNFIntegrateModel.AllInOne;
VNF_opt.RandomSeed = [20161101 0]; 
net_opt.AdmitPolicy = 'reject-flow';
net_opt.Form = 'compact';

%% control variables
b_static = true;
b_optimal = true;
b_repeat = true;

%% Construct Network
% 13 23 33
slice_type = 33;

%% 
if b_optimal
    net_opt.SlicingMethod = SlicingMethod.SingleNormal;
    PN = CloudAccessNetwork(node_opt, link_opt, VNF_opt, net_opt);
    PN.slice_template = Slice.loadSliceTemplate(slice_type);
    link_capacity = PN.getLinkField('Capacity');
    node_capacity = PN.getDataCenterField('Capacity');
    seed = floor(now);
    slice_opt = PN.slice_template(1);
    fprintf('\nSingle Slice Optimization:\n');
    fprintf('\tAverage unit link cost: %.2G, average unit node cost: %.2G.\n', ...
        mean(PN.getLinkField('UnitCost')), ...
        mean(PN.getDataCenterField('UnitCost')));
    fprintf('\t\t(Ratio of unit node cost to unit link cost: %.2G.)\n\n',...
        mean(PN.getDataCenterField('UnitCost'))/mean(PN.getLinkField('UnitCost')));
    N = 5;
    while true && N > 0
        slice_opt.RandomSeed = seed;
        seed = seed + 1;
        PN.AddSlice(slice_opt);
        %     output = PN.singleSliceOptimization();
        output = PN.singleSliceOptimization(struct('bCompact', false));
    
        fprintf('\tNumber of slices: %d.\n', PN.NumberSlices);
        fprintf('\tOptimal net social welfare (without pricing) is %.4e.\n', ...
            output.WelfareOptimal);
        fprintf('\tOptimal net social welfare (with pricing) is %.4e.\n', output.Welfare);
        fprintf('\tnet profit of each slice:\n');
        fprintf('\t\t%f\n',output.Profit(1:(end-1),:));
        fprintf('\tnet profit of substrate network:\n');
        fprintf('\t\t%f\n',output.Profit(end,:));
        fprintf('\tNetwork utilization ratio %f.\n',PN.utilizationRatio);
        fprintf('\t\t(Node utilization: %.2G)\n', sum(PN.getDataCenterField('Load')/sum(node_capacity)));
        fprintf('\t\t(Link utilization: %.2G)\n\n', sum(PN.getLinkField('Load')/sum(link_capacity)));
        if ~b_repeat
            break;
        else
            N = N - 1;
        end
    end
end
%%
if b_static
    net_opt.SlicingMethod = SlicingMethod.StaticPricing;
    PN_static = CloudAccessNetwork(node_opt, link_opt, VNF_opt, net_opt);
    PN_static.slice_template = Slice.loadSliceTemplate(slice_type);
    link_capacity = PN_static.getLinkField('Capacity');
    node_capacity = PN_static.getDataCenterField('Capacity');
    link_price = PN_static.getLinkCost * (1 + net_opt.PricingFactor);
    node_price = PN_static.getNodeCost * (1 + net_opt.PricingFactor);
    PN_static.setLinkField('Price', link_price);
    PN_static.setDataCenterField('Price', node_price);
    slice_opt = PN_static.slice_template(1);
    if b_repeat
        slice_opt.AdmitPolicy = 'reject-slice';
        fprintf('\nStatic Slicing Repeat:\n');
    else
        fprintf('\nStatic Slicing:\n');
    end
    seed = floor(now);
    while true
        slice_opt.RandomSeed = seed;
        seed = seed + 1;
        sl = PN_static.AddSlice(slice_opt);
        if isempty(sl)
            break;
        end
        output = PN_static.staticSlicing(sl);
        fprintf('\tNumber of slices: %d.\n', PN_static.NumberSlices);
        fprintf('\tOptimal net social welfare (with pricing) is %.4e.\n', output.Welfare);
        fprintf('\tnet profit of each slice:\n');
        fprintf('\t\t%f\n',output.Profit(1:(end-1),:));
        fprintf('\tnet profit of substrate network:\n');
        fprintf('\t\t%f\n',output.Profit(end,:));
        fprintf('\tNetwork utilization ratio %f.\n',PN_static.utilizationRatio);
        fprintf('\t\t(Node utilization: %.2G)\n', ...
            sum(PN_static.getDataCenterField('Load')/sum(node_capacity)));
        fprintf('\t\t(Link utilization: %.2G)\n\n', ...
            sum(PN_static.getLinkField('Load')/sum(link_capacity)));
        if ~b_repeat
            break;
        end
    end
    fprintf('\tAverage unit link cost: %.2G, average unit node cost: %.2G.\n', ...
        mean(PN_static.getLinkField('UnitCost')), ...
        mean(PN_static.getDataCenterField('UnitCost')));
    fprintf('\t\t(Ratio of unit node cost to unit link cost: %.2G.)\n\n',...
        mean(PN_static.getDataCenterField('UnitCost'))/...
        mean(PN_static.getLinkField('UnitCost')));
end
