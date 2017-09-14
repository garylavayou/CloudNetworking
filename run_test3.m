%% Physical Network Speification for SD_RAN
% This program used to test if the parameters are properly set.
%% Configurable parameters in this experiment
link_opt.CostUnit = 0.15;
node_opt.CostUnit = 0.1;
options.PricingFactor = 2;      % 1 | 2 | 3 ,
options.Threshold = 'average';

%% Specification of VNFs and Network Slices
% |StaticCostOption| is not set, the default value is |None|;
% |RandomSeed|: the first seed is for random static cost, the second is for process
% efficiency.
link_opt.RandomSeed = 20170421;
link_opt.delay = LinkDelayOption.Random;
net_opt.AdmitPolicy = 'reject-flow';

VNF_opt.Number = 6;            % number of VNF type
VNF_opt.Model = VNFIntegrateModel.AllInOne;
VNF_opt.RandomSeed = [20161101 0]; 

%% Algorithm options
options.Display = 'off';
options.ProfitType = {'ApproximatePrice','AccuratePrice'};
options.WelfareType = {'Accurate', 'Approximate'};

%% control variables
b_static_slice = false;
b_optimal = false;
b_static_slice_repeat = true;

%% Construct Network
% 7 8 9
slice_type = 9;

if b_static_slice
    CNs = CloudNetwork(node_opt, link_opt, VNF_opt, net_opt);
    CNs.slice_template = Slice.loadSliceTemplate(slice_type);
    link_capacity = CNs.getLinkField('Capacity');
    node_capacity = CNs.getDataCenterField('Capacity');
    link_price = CNs.getLinkCost * (1 + options.PricingFactor);
    node_price = CNs.getNodeCost * (1 + options.PricingFactor);
    CNs.setLinkField('Price', link_price);
    CNs.setDataCenterField('Price', node_price);
    slice_opt = CNs.slice_template(1);
    slice_opt.method = 'static-slicing';
    slice_opt.admit_ploicy = 'reject-flow';
    sl = CNs.AddSlice(slice_opt);
    options.Method = 'slice-price';
    output = CNs.staticSlicing(sl, options);

    fprintf('\nStatic Slicing\n');
    fprintf('\tOptimal net social welfare (with pricing) is %.4e(%.4e).\n', ...
        output.welfare_accurate, output.welfare_approx);
    fprintf('\tnet profit of each slice:\n');
    display(output.profit(1:(end-1),:));
    fprintf('\tnet profit of substrate network:\n');
    display(output.profit(end,:));
    fprintf('\tNetwork utilization ratio %f.\n',CNs.utilizationRatio);
    node_load = CNs.getDataCenterField('Load');
    node_index = CNs.DataCenters.NodeIndex(node_load>0);
    fprintf('\t\t(Node utilization: %.2G)\n', ...
        sum(node_load)/sum(node_capacity(node_index)));
    link_load = CNs.getLinkField('Load');
    link_index = link_load>0;
    fprintf('\t\t(Link utilization: %.2G)\n\n', ...
        sum(link_load)/sum(link_capacity(link_index)));
end
%% single slice optimize
if b_optimal
    CN = CloudNetwork(node_opt, link_opt, VNF_opt, net_opt);
    CN.slice_template = Slice.loadSliceTemplate(slice_type);
    link_capacity = CN.getLinkField('Capacity');
    node_capacity = CN.getDataCenterField('Capacity');
    slice_opt = CN.slice_template(1);
    CN.AddSlice(slice_opt);
    options.Method = 'normal';
    output = CN.singleSliceOptimization(options);

    fprintf('\nSingle Slice Optimization:\n');
    fprintf('\tAverage unit link cost: %.2G, average unit node cost: %.2G.\n', ...
        mean(CN.getLinkField('UnitCost')), mean(CN.getDataCenterField('UnitCost')));
    fprintf('\t\t(Ratio of unit node cost to unit link cost: %.2G.)\n',...
        mean(CN.getDataCenterField('UnitCost'))/mean(CN.getLinkField('UnitCost')));
    fprintf('\tOptimal net social welfare (without pricing) is %.4e(%.4e).\n', ...
        output.welfare_accurate_optimal, output.welfare_approx_optimal);
    fprintf('\tOptimal net social welfare (with pricing) is %.4e(%.4e).\n', ...
        output.welfare_accurate, output.welfare_approx);
    fprintf('\tnet profit of each slice:\n');
    display(output.profit(1:(end-1),:));
    fprintf('\tnet profit of substrate network:\n');
    display(output.profit(end,:));
    fprintf('\tNetwork utilization ratio %f.\n',CN.utilizationRatio);
    node_load = CN_static.getDataCenterField('Load');
    node_index = CN_static.data_centers(node_load>0);
    fprintf('\t\t(Node utilization: %.2G)\n', ...
        sum(node_load)/sum(node_capacity(node_index)));
    link_load = CN_static.getLinkField('Load');
    link_index = link_load>0;
    fprintf('\t\t(Link utilization: %.2G)\n\n', ...
        sum(link_load)/sum(link_capacity(link_index)));end
%%
if b_static_slice_repeat
    CN_static = CloudNetwork(node_opt, link_opt, VNF_opt, net_opt);
    CN_static.slice_template = Slice.loadSliceTemplate(slice_type);
    link_capacity = CN_static.getLinkField('Capacity');
    node_capacity = CN_static.getDataCenterField('Capacity');
    link_price = CN_static.getLinkCost * (1 + options.PricingFactor);
    node_price = CN_static.getNodeCost * (1 + options.PricingFactor);
    CN_static.setLinkField('Price', link_price);
    CN_static.setDataCenterField('Price', node_price);
    slice_opt = CN_static.slice_template(1);
    slice_opt.method = 'static-slicing';
    slice_opt.admit_ploicy = 'reject-slice';
    seed = floor(now);
    options.Method = 'slice-price';
    fprintf('\nStatic Slicing Repeat:\n');
    while true
        slice_opt.RandomSeed = seed;
        slice_opt.admit_ploicy = net_opt.AdmitPolicy;
        seed = seed + 1;
        sl = CN_static.AddSlice(slice_opt);
        if isempty(sl)
            break;
        end
        output = CN_static.staticSlicing(sl, options);
        fprintf('\tNumber of slices: %d.\n', CN_static.NumberSlices);
        fprintf('\tOptimal net social welfare (with pricing) is %.4e(%.4e).\n', ...
            output.welfare_accurate, output.welfare_approx);
        fprintf('\tnet profit of each slice:\n');
        display(output.profit(1:(end-1),:));
        fprintf('\tnet profit of substrate network:\n');
        display(output.profit(end,:));
        fprintf('\tNetwork utilization ratio %f.\n',CN_static.utilizationRatio);
        node_load = CN_static.getDataCenterField('Load');
        node_index = CN_static.DataCenters.NodeIndex(node_load>0);
        fprintf('\t\t(Node utilization: %.2G)\n', ...
            sum(node_load)/sum(node_capacity));
        link_load = CN_static.getLinkField('Load');
        link_index = link_load>0;
        fprintf('\t\t(Link utilization: %.2G)\n\n', ...
            sum(link_load)/sum(link_capacity));
    end
    fprintf('\tAverage unit link cost: %.2G, average unit node cost: %.2G.\n', ...
        mean(CN_static.getLinkField('UnitCost')), mean(CN_static.getDataCenterField('UnitCost')));
    fprintf('\t\t(Ratio of unit node cost to unit link cost: %.2G.)\n\n',...
        mean(CN_static.getDataCenterField('UnitCost'))/mean(CN_static.getLinkField('UnitCost')));
    node_load = CN_static.getDataCenterField('Load');
    link_load = CN_static.getLinkField('Load');
    residual_link_capacity = link_capacity - link_load;
end