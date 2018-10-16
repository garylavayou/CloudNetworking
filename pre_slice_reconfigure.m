%% Dynamic Fast Slice Reconfiguration
% Experiment 4: evaluate the performance of fast slice reconfiguration
% ('fastconfig' and 'fastconfig2') and hybrid slicing scheme ('dimconfig'
% and 'dimconfig2'). 

%% Experiment Control
global DEBUG; 
DEBUG = true;
warning('off', 'backtrace');
warning('off', 'verbose');
b_baseline = false;
b_dimbaseline = false;			
b_fastconfig0 = false;
b_fastconfig = false;
b_dimconfig = false;				
b_dimconfig0 = false;			
b_dimconfig1 = false;
b_fastconfig2 = false;
b_dimconfig2 = false;
b_save = false;
b_plot = false;
b_plotsave = false;

%% Debug Information
%{
if ~exist('debug_info', 'var') && exist('breakpoint.mat', 'file')
    load('breakpoint.mat');
    dbstop(debug_info);
else
    debug_info = dbstatus;
    save('breakpoint', 'debug_info');
end
%}

%% Physical Network Specification
clear link_opt node_opt options VNF_opt type ;
options.NetworkType = 'DynamicCloudNetwork';
node_opt.Model = NetworkModel.Sample2;
node_opt.CapacityModel = NodeCapacityOption.NetworkSpecified;
node_opt.Capacity = [0, 4000, 0, 4000, 0, 4000, 0, ...
	0, 0, 6000, 0, 3000, 0, 0, 2000];  % user defined capacity;
node_opt.CapacityFactor = 3;
node_opt.CostModel = NodeCostOption.CapacityInverse;
node_opt.CostUnit = 500;
link_opt.DelayModel = LinkDelayOption.Random;
link_opt.CostModel = LinkCostOption.CapacityInverse;
link_opt.CapacityFactor = 1000;
link_opt.CostUnit = 50;
link_opt.RandomSeed = 20170423;

%% Specification of VNFs and Network Slices
VNF_opt.Number = 6;            % number of VNF type
VNF_opt.Model = VNFIntegrateModel.AllInOne;
VNF_opt.RandomSeed = [20161101 0];
type.Index = [144; 154; 164; 174; 184];
type.Permanent = 4;
type.PermanentCount = 1;  % Number of persistent slices: {1|2|3...}
type.Static = [1; 2; 3];
type.StaticCount = [1; 2; 2];
type.StaticClass = {'SimpleSlice'};

%% Algorithm options
options.PricingFactor = 1;      % 1 | 2 | 3 ,
options.Threshold = 'min';
options.Form = 'compact';       % {'compact'|'normal'}
options.NonzeroTolerance = 10^-3;
options.ConstraintTolerance = 10^-3;
options.DiffNonzeroTolerance = 10^-3;
options.PostProcessing = 'round';
SEED = 20170410;            % 20161231
