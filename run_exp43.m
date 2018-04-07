%% Dynamic Fast Slice Reconfiguration
% using a fixed reconfiguration cost settings ($\theta=0.5$), the other settings are the
% same as see also <run_test42>. 
clear link_opt node_opt options VNF_opt type ; 
if ~exist('debug_info', 'var')
    if exist('breakpoint.mat', 'file')
        load('breakpoint.mat');
    end    
    dbstop(debug_info);
else
  debug_info = dbstatus;
  save('breakpoint', 'debug_info');  
end

%% Physical Network Specification for Sample-2
node_opt.model = NetworkModel.Sample2;
node_opt.capacity = NodeCapacityOption.NetworkSpecified;
node_opt.node_capacity = [0, 4000, 0, 4000, 0, 4000, 0, ...
    0, 0, 6000, 0, 3000, 0, 0, 2000];  % user defined capacity;
node_opt.CapacityFactor = 3;
node_opt.cost = NodeCostOption.CapacityInverse;
node_opt.CostUnit = 500;

link_opt.RandomSeed = 20170423;
link_opt.CapacityFactor = 1000;
link_opt.CostUnit = 300;
link_opt.delay = LinkDelayOption.Random;
link_opt.cost = LinkCostOption.CapacityInverse;
options.NetworkType = 'DynamicCloudNetwork';

%% Specification of VNFs and Network Slices
VNF_opt.Number = 6;            % number of VNF type
VNF_opt.Model = VNFIntegrateModel.AllInOne;
VNF_opt.RandomSeed = [20161101 0];

%% Algorithm options
options.Threshold = 'average';

%% Slice arrival sequence confiuration
% Assume that the arrival interval is fixed, *only adjust the service lifetime and the
% arrival probability* of each type of slices to control the number of slices in the
% network. See also <RequestEvent>.
type.Index = [44; 54; 64]; 
type.Permanent = 1;
type.PermanentCount = 3;      % Number of persistent slices: {1|2|3...}
SEED = 20170410;   % 20161231
NUM_EVENT = 150;             %600;

%% Experiment Control
options.theta = 0.5;
b_baseline = true;
b_fastconfig = true;

%% Run script
% single slice reconfiguration
TOTAL_NUM= NUM_EVENT*(b_fastconfig+b_baseline);
global total_iter_num;
total_iter_num = 0;
if exist('progress_bar', 'var') && isvalid(progress_bar)
    close(progress_bar);
end
progress_bar = waitbar(total_iter_num/TOTAL_NUM, ...
    sprintf('Simulation Progress: %d/%d', total_iter_num, TOTAL_NUM));

if b_fastconfig
    options.ReconfigMethod = 'fastconfig';    % {'reconfig', 'fastconfig', 'dimension', 'fastconfig2'}
    progress_bar.Name = 'Fast Reconfiguration';
    pause(0.01);
    clear functions; %#ok<CLFUNC>
    seed_dynamic = SEED; %#ok<NASGU>
    options.UnitReconfigureCost = etas;
    SingleSliceReconfiguration;
    results.fastconfig = g_results;
end
if b_baseline
    % Only need to compute once for different reconfiguration cost efficient, since the
    % optimization procedure is independent on the cosefficient, and the reconfiguration
    % cost with other coefficient can be derieved from the one results by using the
    % coefficient.
    options.ReconfigMethod = 'reconfig';   
    progress_bar.Name = 'Reconfiguration';
    pause(0.01);
    clear functions; %#ok<CLFUNC>
    seed_dynamic = SEED;
    options.UnitReconfigureCost = etas;
    SingleSliceReconfiguration;
    results.reconfig = g_results;
end
close(progress_bar);

%% 
%data_plot3;
%data_plot31;