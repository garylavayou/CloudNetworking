%% Fast slice reconfiguration
% In this experiment, we only use permanent slices, and thus test only the fast slice
% reconfiguration. See <run_exp6> for testing dynamic slice dimensioning with fast slice
% reconfiguration.

%% Specification of Substrate Network
% The parameters are determined in <run_test1.m>.
link_opt.delay = LinkDelayOption.Random;
link_opt.cost = LinkCostOption.CapacityInverse;
link_opt.CostUnit = 150;
link_opt.CapacityFactor = 30;
net_opt.delta = 0.7;

node_opt.model = NetworkModel.Sample1;
node_opt.capacity = NodeCapacityOption.BandwidthProportion;
node_opt.cost = NodeCostOption.CapacityInverse;
node_opt.CostUnit = 500;
node_opt.CapacityFactor = 1.5;     % [0.3; 0.5; 0.8; 1]
net_opt.ClassName = 'CloudNetwork';

%% Specification of VNFs and Network Slices
% |StaticCostOption| is not set, the default value is |None|;
% |RandomSeed|: the first seed is for random static cost, the second is for process
% efficiency.
VNF_opt.Number = 4;            % number of VNF type
VNF_opt.Model = VNFIntegrateModel.AllInOne;
VNF_opt.RandomSeed = [20161101 0];       

%% Arrival sequence confiuration
% Assume that the arrival interval is fixed, *only adjust the service lifetime and the
% arrival probability* of each type of slices to control the number of slices in the
% network. See also <RequestEvent>.
type.Index = [11; 21; 31];
type.Fixed = [1; 2; 3];
type.FixedCount = [3; 4; 8];      % Number of persistent slices: {1|2|3...}
seed_dynamic = 20161231;
arrival.Number = 100;
arrival.Interval = 0.1;     % Average of arrival interval.

%% Algorithm options
options.PricingFactor = 1;
options.ProfitType = {'ApproximatePrice','AccuratePrice'};
options.WelfareType = {'Accurate', 'Approximate'};

%% Control Variables
declare_info_level('Global', DisplayLevel.Off);
b_single_optimal = true;
b_price_adjust = true;
b_static_slice = true;
% b_dual_decomp = false;
% b_price_adjust2 = false;
% b_resource_part = false;
% b_part_price = false;

%%
PN = instantiateclass(net_opt.ClassName, ...
    node_opt, link_opt, VNF_opt, net_opt);
PN.slice_template = Slice.loadSliceTemplate(type.Index);
for t = 1:num_fix_type
    for s = 1:type.FixedCount(t)
        slice_opt = PN.slice_template(type.Fixed(t));
        slice_opt.RandomSeed = seed_dynamic;
        seed_dynamic = seed_dynamic + 1;
        PN.AddSlice(slice_opt);
    end
end
%%
flow_interval = [0.01, 0.02, 0.03];
flow_arrive_rate = [3000, 1500, 1000];
sum_arrive_rate = sum(flow_arrive_rate);
event_set = struct;
for i = 1:PN.NumberSlices
    event_set(i).ServiceInterval = flow_interval(i);
    event_set(i).Probability = flow_arrive_rate(i)/sum_arrive_rate;
end
RE = RequestEvent(event_set, arrival, seed_dynamic);
%%
RE.reset;
if isfield(arrival, 'StopTime')
    num_events = arrival.StopTime;
else
    num_events = arrival.Number*2;
end
%%
for i = 1:num_events
    e = RE.nextEvent;
    T = RE.countCurrentType;
    disp(T);
    stat.times(i) = e.Time;
    if strcmp(e.Description, 'arrival')
        slice_index = e.Type;
        %% dynamic allocate flow index
        flow_id = 
    else
    end
end