%% Physical Network Speification
%% Specification of Substrate Network
node_opt.model = NetworkModel.Sample1;
node_opt.capacity = NodeCapacityOption.BandwidthProportion;
node_opt.capacity_factor = 1;     % [0.3; 0.5; 0.8; 1]
node_opt.cost = NodeCostOption.Uniform;
node_opt.alpha = 1; % [1; 3]
                    % the ratio of unit node cost to unit link cost.
link_opt.delay = LinkDelayOption.BandwidthPropotion;
link_opt.cost = LinkCostOption.LengthDependent;
link_opt.delay2cost = 0.1;        % 0.1
net_opt.delta = 0.5;

%% Specification of VNFs and Network Slices 
VNF_opt.Number = 6;            % number of VNF type
VNF_opt.Model = VNFIntegrateModel.AllInOne;
VNF_opt.StaticCostOption = NodeStaticCostOption.Random;     % 
VNF_opt.static_cost_range = [0.1 0.3];
% VNF_opt.ProcessEfficiency is not set, using random value.
VNF_opt.RandomSeed = [20161031 20161101];       % the first seed is for random static cost, 
                                                % the second is for process efficiency
options.Display = 'off';
options.PricingFactor = 1;
options.PercentFactor = 0.8;

%% Construct Network
% Initialize substrate network

PN = PhysicalNetwork(node_opt, link_opt, VNF_opt, net_opt);
PN.slice_template = Slice.loadSliceTemplate(1);
link_capacity = PN.getLinkField('Capacity');
node_capacity = PN.getNodeField('Capacity');
node_staticcost = PN.getNodeField('StaticCost');

%% add network slices and run optimize
slice_opt = PN.slice_template(1);
PN.AddSlice(slice_opt);
tic;
[output_optimal, ss] = PN.singleSliceOptimization(options);
toc
fprintf('\nSingle Slice Optimization\n');
fprintf('\tOptimal net social welfare is %.4e(%.4e).\n', ...
    output_optimal.welfare_accurate, output_optimal.welfare_approx);
fprintf('\tnet profit of each slice:\n');
display(output_optimal.profit(1:(end-1),:));
fprintf('\tnet profit of substrate network:\n');
display(output_optimal.profit(end,:));
fprintf('\tNetwork utilization ratio %f.\n\n',PN.utilizationRatio);
