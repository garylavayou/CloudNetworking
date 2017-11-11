%% Dynamic Fast Slice Reconfiguration
% Experiment 4: evaluate the performance of fast slice reconfiguration ('fastconfig' and
% 'fastconfig2') and hybrid slicing scheme ('dimconfig' and 'dimconfig2').
%
% <runexp_4211>:
% # This experiment uses network topology Sample-(2) and slice type-(1), see also
%   <run_test21>.
% # 'fastconfig' and 'fastconfig2' do not support topology change (Adhoc flow is
%   disabled).
% # This experiment has no warm-up phase (1).
preconfig_42xx;
EXPNAME = 'EXP4211';
type.Index = [144; 154; 164];
type.Permanent = 1;
thetas = 1;
% thetas = linspace(0.1, 5, 50);  %thetas = thetas([1 13]); thetas = logspace(log10(0.1),log10(10),50); 
% b_reconfig = false;
% b_fastconfig = false;
% b_fastconfig2 = true;
% b_dimconfig = false;
% b_dimconfig2 = true;
NUM_EVENT = 200;            % the trigger-interval is set to 50. {200}
idx = 1:NUM_EVENT;
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
    'Experiment 4-2-1-1: Fast slice reconfiguration scheme and Hybrid slicing schemes.',...
    'Topology=Sample-2.',...
    'No warm-up phase.',...
    'Slice Type 0144 (disable ad-hoc mode, enable dimension-trigger).');
save('Results\EXP4_OUTPUT211.mat', 'description', 'results', 'NUM_EVENT', 'thetas', ...
    'options', 'node_opt', 'link_opt', 'VNF_opt', 'slice_opt', 'type', 'idx', 'EXPNAME');
%}
% Single experiment.
%{
description = sprintf('%s\n%s\n%s\n%s\n%s',...
    'Experiment 4-2-1-1: Fast slice reconfiguration scheme and Hybrid slicing schemes.',...
    'No warm-up phase.',...
    'Topology=Sample-2.',...
    'eta = 1, period = 20.',...
    'Slice Type 0144 (disable ad-hoc mode, enable dimension-trigger).');
save('Results\singles\EXP4_OUTPUT211s010.mat', 'description', 'results', 'NUM_EVENT', 'thetas', ...
    'options', 'node_opt', 'link_opt', 'VNF_opt', 'slice_opt', 'type', 'idx', 'EXPNAME');
%}