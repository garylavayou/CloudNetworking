%% Dynamic Fast Slice Reconfiguration
% using configuration on SD-RAN, see also <run_test3>
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

%% Configurable parameters in this experiment
link_opt.CostUnit = 0.15;        % 75 | 100 | 150
node_opt.CostUnit = 0.1;        % 250 | 300 | 500
options.PricingFactor = 2;      % 1 | 2 | 3 ,
options.Threshold = 'average';

%% Specification of VNFs and Network Slices
link_opt.RandomSeed = 20170423;
link_opt.delay = LinkDelayOption.Random;
options.AdmitPolicy = 'reject-flow';

VNF_opt.Number = 6;            % number of VNF type
VNF_opt.Model = VNFIntegrateModel.AllInOne;
VNF_opt.RandomSeed = [20161101 0];

%% Algorithm options
options.Display = 'off';
options.ProfitType = {'ApproximatePrice','AccuratePrice'};
options.WelfareType = {'Accurate', 'Approximate'};

%% Arrival sequence confiuration
% Assume that the arrival interval is fixed, *only adjust the service lifetime and the
% arrival probability* of each type of slices to control the number of slices in the
% network. See also <RequestEvent>.
type.Index = [41; 42; 43]; 
type.Permanent = 1;
type.PermanentCount = 3;      % Number of persistent slices: {1|2|3...}
SEED = 20170410;   % 20161231
NUM_EVENT = 150;             %600;

%% Experiment Control
options.Display = 'off';
options.theta = 1/4:1/4:4;
b_reconfig = true;
b_fastconfig = true;

%% Run script
% single slice reconfiguration
NUM_TEST = length(options.theta);
TOTAL_NUM= NUM_EVENT*(NUM_TEST*b_fastconfig+b_reconfig);
total_iter_num = 0;
if exist('progress_bar', 'var') && isvalid(progress_bar)
    close(progress_bar);
end
progress_bar = waitbar(total_iter_num/TOTAL_NUM, ...
    sprintf('Simulation Progress: %d/%d', total_iter_num, TOTAL_NUM));
% pos = progress_bar.OuterPosition;
% pos(3) = 400;
% progress_bar.OuterPosition = pos;

clear results;
if b_fastconfig
    options.Method = 'fastconfig';    % {'reconfig', 'fastconfig', 'dimension', 'fastconfig2'}
    progress_bar.Name = 'Fast Reconfiguration';
    pause(0.01);
    for i = 1:length(options.theta)
        clear functions; %#ok<CLFUNC>
        seed_dynamic = SEED; %#ok<NASGU>
        DynamicSlice.THETA(options.theta(i));
        SingleSliceReconfiguration;
        results.fastconfig(i) = g_results;
    end
end
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
    DynamicSlice.THETA(options.theta(1));
    SingleSliceReconfiguration;
    results.reconfig = g_results;
end
close(progress_bar);

%% 
%data_plot3;
%data_plot31;