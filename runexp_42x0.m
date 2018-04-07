%% Dynamic Fast Slice Reconfiguration
% Experiment 4: evaluate the performance of fast slice reconfiguration ('fastconfig' and
% 'fastconfig2') and hybrid slicing scheme ('dimconfig' and 'dimconfig2').
%
% <runexp_42x0>: 
% # This experiment uses network topology Sample-(2) and slice type-(x), see also
%   <run_test21>.
% # 'reconfig' do not support topology change (Adhoc flow is disabled).
% # This experiment is used to caribrate the nomalizer of the L1-approximation (0). We
%   only test the 'reconfig' method, which is independent of the reconfiguration cost
%   coefficient.
preconfig_42xx;
EXPNAME = 'EXP42x0';
type.Index = [144; 154; 164];
type.Permanent = 1;         %  => change to test different type of slice {1,2,3}
etas = 1;
b_baseline = true;          
b_fastconfig = false;
b_fastconfig2 = false;
b_dimconfig = false;
b_dimconfig2 = false;
NUM_EVENT = 200;             % the trigger-interval is set to 50. {200}
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
    sprintf('Experiment 4-2-%d-0: Find the L1-approximation regulizer (beta=2).', type.Permanent),...
    'Topology=Sample-2.',...
    'No warm-up phase.',...
    sprintf('Slice Type %d (disable ad-hoc mode, enable dimension-trigger).', type.Index(type.Permanent))...
    );
output_name = sprintf('Results/calibrate/EXP4_OUTPUT2%d0.mat', type.Permanent);
save(output_name, 'description', 'results', 'NUM_EVENT', 'etas', ...
    'options', 'node_opt', 'link_opt', 'VNF_opt', 'slice_opt', 'type', 'idx', 'EXPNAME');
%}
