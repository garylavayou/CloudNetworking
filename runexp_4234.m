%% Dynamic Fast Slice Reconfiguration
% Experiment 4: evaluate the performance of fast slice reconfiguration ('fastconfig' and
% 'fastconfig2') and hybrid slicing scheme ('dimconfig' and 'dimconfig2').
%
% <runexp_4233>:
% # This experiment uses network topology Sample-(2) and slice type-(3), see also
%   <run_test21>.
% # 'fastconfig' and 'fastconfig2' do not support topology change (Adhoc flow is
%   disabled).
% # This experiment has a warm-up phase and a number of static slices (4). We reuse the
%   warm-up phase so as to reduce the computation (<runexp_4xxx2>). 
preconfig_42xx;
EXPNAME = 'EXP4234';
type.Index = [144; 154; 164]; 
type.Permanent = 3;
type.Static = [1; 2; 3];
type.StaticCount = [1; 2; 2];
type.StaticClass = {'Slice'};
etas = [1/16 1/8 1/4 1/2 1 2 4 8];  %linspace(0.1, 5, 50);  %etas = etas([11 12]);
options.NumberEventWarmUp = 50;
NUM_EVENT = 100;            % {40|100|600};
runexp_4xxx2;

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
save('Results\EXP4_OUTPUT22.mat', 'description', 'results', 'NUM_EVENT', 'etas', ...
    'options', 'node_opt', 'link_opt', 'VNF_opt', 'slice_opt', 'type');
%}