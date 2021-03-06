%% Dynamic Fast Slice Reconfiguration
% Experiment 4: evaluate the performance of fast slice reconfiguration ('fastconfig' and
% 'fastconfig2') and hybrid slicing scheme ('dimconfig' and 'dimconfig2').
%
% <run_exp4222>:
% # This experiment (4) uses network topology Sample-(2) and slice type-(2), see also
%   <run_test21>.
% # 'fastconfig' and 'fastconfig2' do not support topology change (Adhoc flow is
%   disabled).
% # This experiment has warm-up phase (2).
preconfig_42xx;
EXPNAME = 'EXP4222';
type.Index = [10144; 10154; 10164];
type.Permanent = 2;
etas = linspace(0.1, 5, 50);  %etas = etas([11 12]);
% etas = 1;
NUM_EVENT = 100;            % {40|100|600};
idx = 51:NUM_EVENT;
% b_reconfig = true;
% b_fastconfig = false;
% b_fastconfig2 = false;
% b_dimconfig = false;
% b_dimconfig2 = false;
runexp_4xxx;

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
    'Experiment 4-2-2-2: Fast slice reconfiguration scheme and Hybrid slicing schemes.',...
    'Topology=Sample-2.',...
    'Start with warm-up phase.',...
    'Slice Type 10154 (disable ad-hoc mode, enable dimension-trigger).');
save('Results\EXP4_OUTPUT222.mat', 'description', 'results', 'NUM_EVENT', 'etas', ...
    'options', 'node_opt', 'link_opt', 'VNF_opt', 'slice_opt', 'type', 'idx', 'EXPNAME');
%}
