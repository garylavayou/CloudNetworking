%% Physical Network Speification
%% Specification of Substrate Network

node_opt.model = NetworkModel.Sample1;
node_opt.capacity = NodeCapacityOption.BandwidthProportion;
node_opt.CapacityFactor = 1.5;     % [0.3; 0.5; 0.8; 1.5]
node_opt.cost = NodeCostOption.CapacityInverse;
node_opt.alpha = 1; % [1; 3]
                    % the ratio of unit node cost to unit link cost.
link_opt.delay = LinkDelayOption.BandwidthPropotion;
link_opt.cost = LinkCostOption.LengthDependent;
link_opt.delay2cost = 0.1;        % 0.1
net_opt.delta = 0.7;
%% Specification of VNFs and Network Slices 

VNF_opt.Number = 4;            % number of VNF type
VNF_opt.Model = VNFIntegrateModel.AllInOne;
VNF_opt.StaticCostOption = NodeStaticCostOption.Random;     % 
VNF_opt.static_cost_range = [0.1 0.3];
% VNF_opt.ProcessEfficiency is not set, using random value.
VNF_opt.RandomSeed = [20161101 20161031];       % the first seed is for random static cost, 
                                                % the second is for process efficiency
if node_opt.model == NetworkModel.Sample1
    type.index = [11 21 31];
else
    type.index = [12 22 32];
end
type.count = [1 7 9];           % [1 7 9]
type.partition_weight = [1 1 1];    % [1 2 4]
seed = 20161231;

%% options
options.PricingFactor = 3;
options.PercentFactor = 0.8;
                                                
%% Construct Network
% Initialize substrate network and add network slices.
PN = CloudNetwork(node_opt, link_opt, VNF_opt, net_opt);
PN.slice_template = Slice.loadSliceTemplate(type.index);
link_capacity = PN.getLinkField('Capacity');
node_capacity = PN.getNodeField('Capacity');
node_static_cost = PN.getNodeField('StaticCost');
partition_weight = zeros(dot(type.index,type.count),1);
s = 1;
for i = 1:length(type.index)
%     PN.slice_template(i).Weight = PN.slice_template(i).Weight * 4;
    if PN.slice_template(i).NumberVNFs > VNF_opt.Number
        error_info = strcat('error: confict VNF configuration where ',...
            'the number of requested VNF type is more than the number', ...
            ' of VNF type the network provides');
        error(error_info);
    end
    for j = 1:type.count(i)
        slice_opt = PN.slice_template(i);
        slice_opt.RandomSeed = seed;
        seed = seed + 1;
        PN.AddSlice(slice_opt);
        partition_weight(s) = type.partition_weight(i);
        s = s + 1;
    end
end
% clearvars -except PN node_opt link_opt VNF_opt slice_opt net_opt
% PN.plot;

%% Optimize
% * *Global optimization*
disp('--------- Single Slice Optimization ----------')
declare_info_level('Global', DisplayLevel.Off);
% options.Method = 'single-function';
tic;
[output_optimal, ss] = PN.singleSliceOptimization(options);
toc;
fprintf('\tOptimal net social welfare is %.4e(%.4e).\n', ...
    output_optimal.welfare_accurate, output_optimal.welfare_approx);
fprintf('\tnet profit of each slice:\n');
display(output_optimal.profit(1:(end-1),:));
fprintf('\tnet profit of substrate network:\n');
display(output_optimal.profit(end,:));
fprintf('\tNetwork utilization ratio %f.\n\n',PN.utilizationRatio);
%% 
% * *Dual Decomposition*
% 
% options.Display = 'iter';
% tic
% [output_dual] = PN.optimizeNetSocialWelfare1(options);
% toc
% fprintf('\nOptimization with Dual Decomposition\n')
% fprintf('\tNet social welfare is %.4e(%.4e).\n', ...
%     output_dual.welfare_accurate, output_dual.welfare_approx);
% fprintf('\tnet profit of each slice:\n');
% display(output_dual.profit(1:(end-1),:));
% fprintf('\tnet profit of substrate network:\n');
% display(output_dual.profit(end,:));
% fprintf('\tNetwork utilization ratio %f.\n\n',PN.utilizationRatio);
%%
% * *Dynamic Pricing*
disp('--------- Optimization with Adjusting Resource Price ----------')
tic
[output_price] = PN.optimizeResourcePrice([], options);
toc
fprintf('\tNet social welfare is %.4e(%.4e).\n', ...
    output_price.welfare_accurate, output_price.welfare_approx);
fprintf('\tnet profit of each slice:\n');
display(output_price.profit(1:(end-1),:));
fprintf('\tnet profit of substrate network:\n');
display(output_price.profit(end,:));
fprintf('\tNetwork utilization ratio %f.\n\n',PN.utilizationRatio);
%% 
disp('--------- Optimization with Adjusting Resource Price (2) ----------')
declare_info_level('Global', DisplayLevel.Off);
tic
[output_price2] = PN.optimizeResourcePrice2([], options);
toc
fprintf('\tNet social welfare is %.4e(%.4e).\n', ...
    output_price2.welfare_accurate, output_price2.welfare_approx);
fprintf('\tnet profit of each slice:\n');
display(output_price2.profit(1:(end-1),:));
fprintf('\tnet profit of substrate network:\n');
display(output_price2.profit(end,:));
fprintf('\tNetwork utilization ratio %f.\n\n',PN.utilizationRatio);
%% 
disp('--------- Optimization with Adjusting Resource Price (3) ----------')
declare_info_level('Global', DisplayLevel.Off);
tic
[output_price3] = PN.optimizeResourcePrice3([], options);
toc
fprintf('\tNet social welfare is %.4e(%.4e).\n', ...
    output_price3.welfare_accurate, output_price3.welfare_approx);
fprintf('\tnet profit of each slice:\n');
display(output_price3.profit(1:(end-1),:));
fprintf('\tnet profit of substrate network:\n');
display(output_price3.profit(end,:));
fprintf('\tNetwork utilization ratio %f.\n\n',PN.utilizationRatio);

%% 
% * *Resource Partition*
disp('--------- Optimization with Resource Partition ----------')
declare_info_level('Global', DisplayLevel.Off);
tic
[output_part] = PN.resourcePartitionOptimization(partition_weight, options);
% [output_part] = PN.resourcePartitionOptimization([], options);
toc
fprintf('\tNet social welfare is %.4e(%.4e).\n', ...
    output_part.welfare_accurate, output_part.welfare_approx);
fprintf('\tnet profit of each slice:\n');
display(output_part.profit(1:(end-1),:));
fprintf('\tnet profit of substrate network:\n');
display(output_part.profit(end,:));
fprintf('\tNetwork utilization ratio %f.\n\n',PN.utilizationRatio);

%%
% * *Price Adjustment Based on Resource Partition*
disp('--------- Price Adjustment Based on Resource Partition ----------')
declare_info_level('Global', DisplayLevel.Off);
tic
[output_partprice] = PN.partitionResourcePricing([], options);
toc
fprintf('\tNet social welfare is %.4e(%.4e).\n', ...
    output_partprice.welfare_accurate, output_partprice.welfare_approx);
fprintf('\tnet profit of each slice:\n');
display(output_partprice.profit(1:(end-1),:));
fprintf('\tnet profit of substrate network:\n');
display(output_partprice.profit(end,:));
fprintf('\tNetwork utilization ratio %f.\n\n',PN.utilizationRatio);
%%
bar([output_optimal.welfare_approx, output_price.welfare_approx, ...
    output_price2.welfare_approx, output_part.welfare_approx]);