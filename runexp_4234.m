%% Dynamic Fast Slice Reconfiguration
% Experiment 4: evaluate the performance of fast slice reconfiguration ('fastconfig' and
% 'fastconfig2') and hybrid slicing scheme ('dimconfig' and 'dimconfig2').
%
% <runexp_4233>:
% # This experiment uses network topology Sample-(2) and slice type-(3), see also
%   <run_test21>.
% # 'fastconfig' and 'fastconfig2' do not support topology change (Adhoc flow is
%   disabled).
% # This experiment has a warm-up phase and a number of static slices (4). We should try to reuse
%   the warm-up phase so as to reduce the computation.
preconfig_42xx;
EXPNAME = 'EXP4234';
type.Index = [10144; 10154; 10164]; 
type.Permanent = 3;
type.Static = [1; 2; 3];
type.StaticCount = [1; 2; 2];
type.StaticClass = {'CloudNetwork','CloudNetwork','CloudNetwork'};
thetas = linspace(0.1, 5, 50);  %thetas = thetas([11 12]);
options.NumberEventWarmUp = 50;
NUM_EVENT = 100;            % {40|100|600};

%% Output
% # plot figure
%{
data_plot3;
data_plot31;
%}
% # Save Results
%{
description = sprintf('%s\n%s\n%s',...
    'Experiment 4-2: Fast slice reconfiguration scheme and Hybrid slicing schemes.',...
    'Topology=Sample-2.',...
    'Slice Type 0154 (disable ad-hoc mode, enable trigger).');
save('Results\EXP4_OUTPUT22.mat', 'description', 'results', 'NUM_EVENT', 'thetas', ...
    'options', 'node_opt', 'link_opt', 'VNF_opt', 'slice_opt', 'type');
%}