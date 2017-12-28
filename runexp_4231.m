%% Dynamic Fast Slice Reconfiguration
% Experiment 4: evaluate the performance of fast slice reconfiguration ('fastconfig' and
% 'fastconfig2') and hybrid slicing scheme ('dimconfig' and 'dimconfig2').
%
% <runexp_4231>:
% # This experiment uses network topology Sample-(2) and slice type-(3), see also
%   <run_test21>.
% # 'fastconfig' and 'fastconfig2' do not support topology change (Adhoc flow is
%   disabled).
% # This experiment has no warm-up phase (1).
preconfig_42xx;
EXPNAME = 'EXP4231';
type.Index = [144; 154; 164];
type.Permanent = 3;
etas = linspace(0.1, 5, 50);  %etas = etas([1 13]); logspace(log10(0.1),log10(10),50)
NUM_EVENT = 100;            % {40|100|600};
idx = 1:NUM_EVENT;
% b_reconfig = false;
% b_fastconfig = false;
% b_fastconfig2 = false;
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
    'Experiment 4-2-3-1: Fast slice reconfiguration scheme and Hybrid slicing schemes.',...
    'Topology=Sample-2.',...
    'No warm-up phase.',...
    'Slice Type 0164 (disable ad-hoc mode, enable dimension-trigger).');
save('Results\EXP4_OUTPUT231.mat', 'description', 'results', 'NUM_EVENT', 'etas', ...
    'options', 'node_opt', 'link_opt', 'VNF_opt', 'slice_opt', 'type', 'idx',...
    'EXPNAME');
%}
