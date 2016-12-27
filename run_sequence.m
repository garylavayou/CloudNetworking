%% Physical Network Speification
%% Specification of Substrate Network
node_opt.model = NetworkModel.Sample2;  % NetworkModel.Sample1;
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
VNF_opt.Number = 8;            % number of VNF type
VNF_opt.Model = VNFIntegrateModel.AllInOne;
VNF_opt.StaticCostOption = NodeStaticCostOption.Random;     % 
VNF_opt.static_cost_range = [0.1 0.3];
% VNF_opt.ProcessEfficiency is not set, using random value.
VNF_opt.RandomSeed = [20161031 20161101];       % the first seed is for random static cost, 
                                                % the second is for process efficiency
                                                
%% Construct Network
% Initialize substrate network
PN = PhysicalNetwork(node_opt, link_opt, VNF_opt, net_opt);
PN.slice_template = Slice.loadSliceTemplate([4 5 6]);

%% Network Slice as a Poisson Arrvial Process
% Three types of network slice as a compound Poisson process.
arrival.Number = 30;
num_events = arrival.Number*2;
arrival.Interval = 1;
event_set = struct;
for i = 1:length(PN.slice_template)
    event_set(i).ServiceInterval = PN.slice_template(i).ServiceInterval;
    event_set(i).Probability = PN.slice_template(i).Probability;
end
seed=20161231;
RE = RequestEvent(event_set, arrival, seed);

%% Running Options
options.Display = 'off';
options.PricingFactor = 1;
options.PercentFactor = 0.8;
b_single_optimal = false;
b_dual_decomp = false;
b_price_adjust = false;
b_resource_part = false;
if b_single_optimal
    profit_approx.optimal = zeros(num_events,1); 
    profit_accurate.optimal = zeros(num_events,1);
    utilization.optimal = zeros(num_events,1);
    runtime.optimal = zeros(num_events,1);
end
if b_dual_decomp
    profit_approx.dual = zeros(num_events,1); %#ok<UNRCH>
    profit_accurate.dual = zeros(num_events,1);
    utilization.dual = zeros(num_events,1);
    runtime.dual = zeros(num_events,1);
end
if b_price_adjust
    profit_approx.price = zeros(num_events,1);
    profit_accurate.price = zeros(num_events,1);
    utilization.price = zeros(num_events,1);
    runtime.price = zeros(num_events,1);
end
if b_resource_part
    profit_approx.part = zeros(num_events,1); 
    profit_accurate.part = zeros(num_events,1);
    utilization.part = zeros(num_events,1);
    runtime.part = zeros(num_events,1);
end
stat.number_slices = zeros(num_events,3);
stat.number_flows = zeros(num_events,1);
RE.reset;
for i = 1:num_events
    e = RE.nextEvent;
    T = RE.countCurrentType;
    display(T);
    stat.number_slices(i,T.Value) = T.Count;
    if strcmp(e.Description, 'arrival')
        fprintf('\tSlice type %d.\n', e.Type);
        slice_opt = PN.slice_template(e.Type);
        slice_opt.RandomSeed = e.RandomSeed;
        slice_opt.Identifier = e.Identifier;
        PN.AddSlice(slice_opt);
        stat.number_flows(i) = PN.NumberFlows;
        if b_single_optimal
            tic;
            [output_optimal, ss] = PN.singleSliceOptimization(options); 
            runtime.optimal(i) = toc;
            profit_approx.optimal(i) = output_optimal.welfare_approx;
            profit_accurate.optimal(i) = output_optimal.welfare_accurate;
            utilization.optimal(i) = PN.utilizationRatio;
        end
        if b_dual_decomp
            tic;
            [output_dual] = PN.optimizeNetSocialWelfare1(options); 
            runtime.dual(i) = toc;
            profit_approx.dual(i) = output_dual.welfare_approx;
            profit_accurate.dual(i) = output_dual.welfare_accurate;
            utilization.dual(i) = PN.utilizationRatio;
        end
        if b_price_adjust
            tic;
            [output_price] = PN.optimizeResourcePrice([], options); 
            runtime.price(i) = toc;
            profit_approx.price(i) = output_price.welfare_approx;
            profit_accurate.price(i) = output_price.welfare_accurate;
            utilization.price(i) = PN.utilizationRatio;
        end
        if b_resource_part
            tic;
            [output_part] = PN.resourcePartitionOptimization([], options); 
            runtime.part(i) = toc;
            profit_approx.part(i) = output_part.welfare_approx;
            profit_accurate.part(i) = output_part.welfare_accurate;
            utilization.part(i) = PN.utilizationRatio;
        end
    else
        PN.RemoveSlice(e.Id);
        stat.number_flows(i) = PN.NumberFlows;
        if PN.NumberSlices > 0
            if b_single_optimal
                [output_optimal, ss] = PN.singleSliceOptimization(options); 
                profit_approx.optimal(i) = output_optimal.welfare_approx;
                profit_accurate.optimal(i) = output_optimal.welfare_accurate;
                utilization.optimal(i) = PN.utilizationRatio;
            end
            if b_dual_decomp
                [output_dual] = PN.optimizeNetSocialWelfare1(options); %#ok<UNRCH>
                profit_approx.dual(i) = output_dual.welfare_approx;
                profit_accurate.dual(i) = output_dual.welfare_accurate;
                utilization.dual(i) = PN.utilizationRatio;
            end
            if b_price_adjust
                [output_price] = PN.optimizeResourcePrice([], options); 
                profit_approx.price(i) = output_price.welfare_approx;
                profit_accurate.price(i) = output_price.welfare_accurate;
                utilization.price(i) = PN.utilizationRatio;
            end
            if b_resource_part
                [output_part] = PN.resourcePartitionOptimization([], options); 
                profit_approx.part(i) = output_part.welfare_approx;
                profit_accurate.part(i) = output_part.welfare_accurate;
                utilization.part(i) = PN.utilizationRatio;
            end
        end
    end
end
display(RE.countArriveType);
%%
x = 10:(num_events-10);
if b_single_optimal
    %     plot(x, profit_approx.optimal(x), x, profit_accurate.optimal(x));
end
if b_price_adjust
%     plot(x, profit_approx.price(x), x, profit_accurate.price(x)); 
end
%%
subplot(2,3,1);
plot(x, profit_accurate.optimal(x),'-r', ...
    x, profit_accurate.price(x), '-b', ...
    x, profit_accurate.part(x), '-g'); 
legend({'optimal', 'price', 'partition'}, 'Location', 'northwest');
xlabel('Network Events (slice arrival/departure)');
ylabel('Net Social Welfare');
title('Network social welfare from different methods');
h = subplot(2,3,2);
devi1 = profit_accurate.optimal(x) - profit_accurate.price(x);
devi2 = profit_accurate.optimal(x) - profit_accurate.part(x);
plot(x, devi1,'-r', x, devi2, '-b'); 
lim = h.YLim;
lim(1) = 0;
h.YLim = lim;
legend({'diff\_price', 'diff\_part'}, 'Location', 'northwest')
xlabel('Network Events (slice arrival/departure)');
ylabel('Difference');
title('Distance to optimal net social welfare');
h = subplot(2,3,3);
plot(x, utilization.optimal(x), '-r', x, utilization.price(x), '-b', ...
    x, utilization.part(x), '-g');
lim = h.YLim;
lim(2) = 1;
h.YLim = lim;
legend({'optimal', 'price', 'partition'},'Location', 'northeast')
xlabel('Network Events (slice arrival/departure)');
ylabel('Network Utilization');
title('Network utilization from different methods');
subplot(2,3,4)
plot(x, stat.number_slices(x,1), '-g', x, stat.number_slices(x,2), '-b', ...
    x, stat.number_slices(x,3), '-r', x, sum(stat.number_slices(x,:),2), '-k');
legend({'Type-1', 'Type-2', 'Type-3', 'Total'},'Location', 'northwest')
xlabel('Network Events (slice arrival/departure)');
ylabel('Number of Slices');
title('Number of Slices');
subplot(2,3,5)
plot(x, stat.number_flows(x,1), '-b');
xlabel('Network Events (slice arrival/departure)');
ylabel('Number of Flows');
title('Total number of flows in the network');