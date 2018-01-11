%% Dynamic Fast Slice Reconfiguration
% Preconfigurtion for experiment 42xx: evaluate the performance of fast slice
% reconfiguration ('fastconfig' and 'fastconfig2') and hybrid slicing scheme ('dimconfig'
% and 'dimconfig2'). 
% # This experiment uses network topology Sample-(2) see also <run_test21>.
clear link_opt node_opt options VNF_opt type ;
if ~exist('debug_info', 'var') && exist('breakpoint.mat', 'file')
    load('breakpoint.mat');
    dbstop(debug_info);
else
    debug_info = dbstatus;
    save('breakpoint', 'debug_info');
end
EXPNAME = 'EXP42XX';

%% Physical Network Specification for Sample-2
options.NetworkType = 'DynamicCloudNetwork';

link_opt.DelayModel = LinkDelayOption.Random;
link_opt.CostModel = LinkCostOption.CapacityInverse;
link_opt.CapacityFactor = 1000;
link_opt.CostUnit = 50;
link_opt.RandomSeed = 20170423;

node_opt.Model = NetworkModel.Sample2;
node_opt.CapacityModel = NodeCapacityOption.NetworkSpecified;
node_opt.Capacity = [0, 4000, 0, 4000, 0, 4000, 0, ...
    0, 0, 6000, 0, 3000, 0, 0, 2000];  % user defined capacity;
node_opt.CapacityFactor = 3;
node_opt.CostModel = NodeCostOption.CapacityInverse;
node_opt.CostUnit = 500;

%% Specification of VNFs and Network Slices
VNF_opt.Number = 6;            % number of VNF type
VNF_opt.Model = VNFIntegrateModel.AllInOne;
VNF_opt.RandomSeed = [20161101 0];
type.Index = [44; 54; 64];
type.Permanent = 1;
type.PermanentCount = 1;      % Number of persistent slices: {1|2|3...}

%% Algorithm options
options.PricingFactor = 1;      % 1 | 2 | 3 ,
options.Threshold = 'min';
options.Form = 'compact';       % {'compact'|'normal'}
options.NonzeroTolerance = 10^-4;
options.ConstraintTolerance = 10^-3;
options.DiffNonzeroTolerance = 10^-3;
options.PostProcessing = 'round';

%% Experiment Control
warning('off', 'backtrace');
warning('off', 'verbose');
b_reconfig = true;
b_fastconfig = true;
b_fastconfig2 = true;
b_dimconfig = true;
b_dimconfig2 = true;
b_dimconfig0 = false;
SEED = 20170410;            % 20161231
