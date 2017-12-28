%% Dynamic Fast Slice Reconfiguration
% using configuration on Sample-2, see also <run_test21>
% TO BE REMOVED.
clear link_opt node_opt options VNF_opt type ;
clearvars -global InfoLevel; 
if ~exist('debug_info', 'var') && exist('breakpoint.mat', 'file')
        load('breakpoint.mat');
        dbstop(debug_info);
else
  debug_info = dbstatus;
  save('breakpoint', 'debug_info');  
end

%% Physical Network Specification for Sample-2
% See also <run_test21>
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

options.NetworkType = 'DynamicCloudNetwork';

%% Specification of VNFs and Network Slices
VNF_opt.Number = 6;            % number of VNF type
VNF_opt.Model = VNFIntegrateModel.AllInOne;
VNF_opt.RandomSeed = [20161101 0];

%% Slice arrival sequence confiuration
% Assume that the arrival interval is fixed, *only adjust the service lifetime and the
% arrival probability* of each type of slices to control the number of slices in the
% network. See also <RequestEvent>.
% These types do not support Adhoc Mode, and they do not support Dimensioning-Trigger mechanisms.
% Use <run_exp421> to evaluate Adhoc Mode and Dimensioning-Trigger.
type.Index = [44; 54; 64]; 
type.Permanent = 3;
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
declare_info_level({'Global', 'ClassDebug'}, [DisplayLevel.Notify, DisplayLevel.Notify]);
warning('off', 'backtrace');
warning('off', 'verbose');
etas = linspace(1, 4, 30);  %etas = etas([11 12]);
b_reconfig = true;
b_fastconfig = true;
b_fastconfig2 = true;
SEED = 20170410;            % 20161231
NUM_EVENT = 100;            % {40|100|600};

%% Run script
% single slice reconfiguration
NUM_TEST = length(etas);
TOTAL_NUM= NUM_EVENT*(NUM_TEST*(b_fastconfig+b_fastconfig2)+b_reconfig);
global total_iter_num;
total_iter_num = 0;
if exist('progress_bar', 'var') && isvalid(progress_bar)
    close(progress_bar);
end
progress_bar = waitbar(total_iter_num/TOTAL_NUM, ...
    sprintf('Simulation Progress: %d/%d', total_iter_num, TOTAL_NUM));
jframe=getJFrame(progress_bar); 
jframe.setAlwaysOnTop(1);
% WindowAPI(progress_bar, 'topmost');
%%
if b_fastconfig
    options.Method = 'fastconfig';    % {'reconfig', 'fastconfig', 'dimension', 'fastconfig2'}
    progress_bar.Name = 'Fast Reconfiguration';
    pause(0.01);
    for i = 1:length(etas)
        clear functions; %#ok<CLFUNC>
        seed_dynamic = SEED;  %#ok<NASGU>
        DynamicSlice.ETA(etas(i));
        SingleSliceReconfiguration;
        if i == 1
            results.Fastconfig = {g_results};
        else
            results.Fastconfig{i} = g_results;
        end
    end
end

%%
if b_fastconfig2
    options.Method = 'fastconfig2';    % {'reconfig', 'fastconfig', 'dimension', 'fastconfig2'}
    progress_bar.Name = 'Fast Reconfiguration 2';
    pause(0.01);
    for i = 1:length(etas)
        clear functions; %#ok<CLFUNC>
        seed_dynamic = SEED;  %#ok<NASGU>
        DynamicSlice.ETA(etas(i));
        SingleSliceReconfiguration;
        if i == 1
            results.Fastconfig2 = {g_results};
        else
            results.Fastconfig2{i} = g_results;
        end
    end
end

%%
if b_reconfig
    % Only need to compute once for different reconfiguration cost efficient, since the
    % optimization procedure is independent on the cosefficient, and the reconfiguration
    % cost with other coefficient can be derieved from the one results by using the
    % coefficient.
    options.Method = 'reconfig';   
    progress_bar.Name = 'Reconfiguration';
    pause(0.01);
    clear functions; %#ok<CLFUNC>
    seed_dynamic = SEED;
    DynamicSlice.ETA(1);
    SingleSliceReconfiguration;
    results.Reconfig = g_results;
end

close(progress_bar);

%% 
%data_plot3;
%data_plot31;
% description = 'Experiment 4-2: topology = Sample-2.';
% save('Results\EXP4_OUTPUT02.mat', 'results', 'NUM_EVENT', 'etas', 'description');