%% Dynamic Fast Slice Reconfiguration
% Experiment 4: evaluate the performance of fast slice reconfiguration ('fastconfig' and
% 'fastconfig2') and hybrid slicing scheme ('dimconfig' and 'dimconfig2').
%
% <runexp_4232>:
% # This experiment uses network topology Sample-(2) and slice type-(3), see also
%   <run_test21>.
% # 'fastconfig' and 'fastconfig2' do not support topology change (Adhoc flow is
%   disabled).
% # This experiment has a warm-up phase (2). Two possible procedure: 
%     a. reuse the warm-up phase so as to reduce the computation (<runexp_4xxx2>). 
%     b. warm-up phase is separately executed for each experiment (slice-type should specify the
%        'test' digit, <runexp_4xxx>). 
preconfig_42xx;
EXPNAME = 'EXP4232';
type.Index = [00144; 00154; 00164];
type.Permanent = 3;
options.NumberEventWarmUp = 50;
thetas = linspace(0.1, 5, 5);  %thetas = thetas([1 13]); thetas = logspace(log10(0.01),log10(10), 50); 
NUM_EVENT = 50;            % {40|100|600};
% idx = (options.NumberEventWarmUp+1):NUM_EVENT;
% runexp_4xxx;
idx = 1:NUM_EVENT;
runexp_4xxx2;

%% Output
% # plot figure
%{
data_plot3;
data_plot31;
dataplot3s;
%}
% # Save Results
%{
description = sprintf('%s\n%s\n%s\n%s',...
    'Experiment 4-2-3-2: Fast slice reconfiguration scheme and Hybrid slicing schemes.',...
    'Topology=Sample-2.',...
    'Slice Type 10164 (disable ad-hoc mode, enable dimension-trigger).',...
    'Start with warm-up phase.');
save('Results\EXP4_OUTPUT232.mat', 'description', 'results', 'NUM_EVENT', 'thetas', ...
    'options', 'node_opt', 'link_opt', 'VNF_opt', 'slice_opt', 'type', 'idx');
%}