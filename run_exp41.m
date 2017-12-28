%% Dynamic Fast Slice Reconfiguration
% using configuration on SD-RAN, see also <run_test3>
clear link_opt node_opt options VNF_opt type ; 
clearvars -global InfoLevel; 
if ~exist('debug_info', 'var')
    if exist('breakpoint.mat', 'file')
        load('breakpoint.mat');
        dbstop(debug_info);
    end    
else
  debug_info = dbstatus;
  save('breakpoint', 'debug_info');  
end

%% Specification of VNFs and Network Slices
link_opt.DelayModel = LinkDelayOption.Random;
link_opt.CostUnit = 0.1;        % 75 | 100 | 150
link_opt.RandomSeed = 20170423;
node_opt.CostUnit = 0.2;        % 250 | 300 | 500

options.NetworkType = 'DynamicCloudAccessNetwork';
options.AdmitPolicy = 'reject-flow';
% options.VNFReconfigCoefficient = 3;

VNF_opt.Number = 6;            % number of VNF type
VNF_opt.Model = VNFIntegrateModel.AllInOne;
VNF_opt.RandomSeed = [20161101 0];

%% Slice arrival sequence confiuration
% Assume that the arrival interval is fixed, *only adjust the service lifetime and the
% arrival probability* of each type of slices to control the number of slices in the
% network. See also <RequestEvent>.
type.Index = [14; 24; 34]; 
type.Permanent = 1;
type.PermanentCount = 3;      % Number of persistent slices: {1|2|3...}

%% Algorithm options
options.PricingFactor = 1;      % 1 | 2 | 3 ,
options.Threshold = 'min';
options.Form = 'compact';
options.NonzeroTolerance = 10^-4;
options.ConstraintTolerance = 10^-3;

%% Experiment Control
declare_info_level({'Global', 'ClassDebug'}, [DisplayLevel.Notify, DisplayLevel.Notify]);
etas = linspace(0.5, 3, 10);
b_reconfig = true;
b_fastconfig = true;
b_fastconfig2 = true;
SEED = 20170410;   % 20161231
NUM_EVENT = 100;             %600;

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
            results.Fastconfig = g_results;
        else
            results.Fastconfig(i) = g_results;
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
            results.Fastconfig2 = g_results;
        else
            results.Fastconfig2(i) = g_results;
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
description = 'Experiment 4-1: topology = SD-RAN.';
save('Results\EXP4_OUTPUT01.mat', 'results', 'NUM_EVENT', 'etas', 'description');