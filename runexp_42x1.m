%% Dynamic Fast Slice Reconfiguration
% Experiment 4: evaluate the performance of fast slice reconfiguration ('fastconfig' and
% 'fastconfig2') and hybrid slicing scheme ('dimconfig' and 'dimconfig2').
%
% <runexp_42x1>:
% # This experiment uses network topology Sample-(2) and slice type-(x), see also
%   <run_test21>.
% # 'fastconfig' and 'fastconfig2' do not support topology change (Adhoc flow is
%   disabled).
% # This experiment has no warm-up phase (1).
preconfig_42xx;
type.Index = [144; 154; 164];
type.Permanent = 1;
type.Static = [1; 2; 3];
type.StaticCount = [1; 2; 2];
type.StaticClass = {'Slice'};
% thetas = [1/16 1/8 1/4 1/2 1 2 4 8];  
thetas = 1/2;
% thetas = linspace(0.1, 5, 50);  thetas = logspace(log10(0.1),log10(10),50); 
b_reconfig = false;
b_fastconfig = false;
b_fastconfig2 = false;
b_dimconfig = false;
b_dimconfig2 = false;
b_dimconfig0 = true;
NUM_EVENT = 100;            % the trigger-interval is set to 50. {200}
idx = 1:NUM_EVENT;
EXPNAME = sprintf('EXP42%d1', type.Permanent);
runexp_4xxx;

%% Output
% # plot figure
%{
data_plot3;
data_plot31;
dataplot3s;
dataplot3sa;
%}
% # Save Results
%{
description = sprintf('%s\n%s\n%s\n%s',...
    sprintf('Experiment 42%d1: Fast slice reconfiguration scheme and Hybrid slicing schemes.', type.Permanent),...
    'Topology=Sample-2.',...
    'Experiment without warm-up phase.',...
    sprint('Slice Type %d (disable ad-hoc mode, enable dimension-trigger).', type.Index(type.Permanent))...
    );
output_name = sprintf('Results\\EXP4_OUTPUT2%d1.mat', type.Permanent);
save(output_name, 'description', 'results', 'NUM_EVENT', 'thetas', 'options', 'node_opt', ...
    'link_opt', 'VNF_opt', 'slice_opt', 'type', 'idx', 'EXPNAME');
%}
% Single experiment.
%{
description = sprintf('%s\n%s\n%s\n%s\n%s',...
    sprintf('Experiment 42%d1: Fast slice reconfiguration scheme and Hybrid slicing schemes.', type.Permanent),...
    'Topology=Sample-2.',...
    'Experiment without warm-up phase.',...
    'Slice Type 0144 (disable ad-hoc mode, enable dimension-trigger).');
output_name = sprintf('Results\\singles\\EXP4_OUTPUT2%d1s%04d.mat', type.Permanent, round(thetas*100));
save(output_name, 'description', 'results', 'NUM_EVENT', 'thetas', ...
    'options', 'node_opt', 'link_opt', 'VNF_opt', 'slice_opt', 'type', 'idx', 'EXPNAME');
%}