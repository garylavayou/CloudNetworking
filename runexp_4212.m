%% Dynamic Fast Slice Reconfiguration
% Experiment 4: evaluate the performance of fast slice reconfiguration ('fastconfig' and
% 'fastconfig2') and hybrid slicing scheme ('dimconfig' and 'dimconfig2').
%
% <runexp_4212>:
% # This experiment uses network topology Sample-(2) and slice type-(1), see also
%   <run_test21>.
% # 'fastconfig' and 'fastconfig2' do not support topology change (Adhoc flow is
%   disabled).
% # This experiment has a warm-up phase (2).
preconfig_42xx;
EXPNAME = 'EXP4212';
type.Index = [144; 154; 164];
type.Permanent = 1;
etas = linspace(0.1, 5, 50);  %etas = etas([1 13]); etas = logspace(log10(0.1),log10(10),50); 
options.NumberEventWarmUp = 50;
% b_reconfig = false;
% b_fastconfig = false;
% b_fastconfig2 = true;
% b_dimconfig = false;
% b_dimconfig2 = true;
NUM_EVENT = 200;            % the trigger-interval is set to 50. {200}
idx = 1:NUM_EVENT;
% runexp_4xxx;
runexp_4xxx2;

%% Output
% # plot figure
%{
data_plot3;
data_plot31;
%}
% # Save Results
%{
description = sprintf('%s\n%s\n%s\n%s',...
    'Experiment 4-2-1-2: Fast slice reconfiguration scheme and Hybrid slicing schemes.',...
    'Topology=Sample-2.',...
    'Experiment has a warm-up phase.',...
    'Slice Type 0144 (disable ad-hoc mode, enable dimension-trigger).');
save('Results\EXP4_OUTPUT212.mat', 'description', 'results', 'NUM_EVENT', 'etas', ...
    'options', 'node_opt', 'link_opt', 'VNF_opt', 'slice_opt', 'type', 'idx', 'EXPNAME');
%}