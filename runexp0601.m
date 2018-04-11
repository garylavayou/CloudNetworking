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
global DEBUG;
global computime ITER_LIMIT;
DEBUG = true;
type.Index = [144; 154; 164; 174; 184];
type.Permanent = 4;
type.Static = [1; 2; 3];
type.StaticCount = [1; 2; 2];
type.StaticClass = {'Slice'};
mode = 'var-penalty'; etas = 1; numberflow = 1000; weight = 30; penalty = [2;4;8;12;16;20];   % [1;2;4;8;12;16]
b_dimconfig0 = true;        % HSR
NUM_EVENT = 11;           
idx = 1:NUM_EVENT;
computime = zeros(max(length(penalty),1),1);
ITER_LIMIT = inf;    % 20
EXPNAME = sprintf('EXP6');
if exist('penalty', 'var') && ~isempty(penalty) && license('test', 'Distrib_Computing_Toolbox') 
    p = gcp;
    if isempty(p)
        warning('Parallel computing is disabled.');
    else
        fprintf('Starting up parallel computing pool, number of workers: %d.\n', p.NumWorkers);
    end
end
runexp04xxx;


%% Output
% # plot figure
%{
data_plot3;
data_plot31;
dataplot3s;
dataplot3sa;
%}
% # Save Results
% A group of experiments with varying parameters
%{
varname = split(mode, '-'); varname = varname{2};
description = sprintf('%s\n%s\n%s',...
    'Experiment 5001: verify the influence of network settings to HSR (without warm-up phase).',...
    sprintf('Topology=Sample-2, SliceType = %d (disable ad-hoc mode, enable dimensioning).', type.Index(type.Permanent)),...
    sprintf('variables = %s', varname)...
    );
output_name = sprintf('Results/EXP5001_var%s.mat', varname);
save(output_name, 'description', 'results', 'EXPNAME', 'mode', 'etas', 'numberflow', 'weight', ...
    'options', 'node_opt', 'link_opt', 'VNF_opt', 'slice_opt', 'type', 'NUM_EVENT');
%}
% Single experiment.
%{
description = sprintf('%s\n%s\n%s\n%s\n%s',...
    sprintf('Experiment 42%d1: Fast slice reconfiguration scheme and Hybrid slicing schemes.', type.Permanent),...
    'Topology=Sample-2.',...
    'Experiment without warm-up phase.',...
    sprintf('Slice Type %d (disable ad-hoc mode, enable dimension-trigger).', type.Index(type.Permanent))...
    );
output_name = sprintf('Results/singles/EXP4_OUTPUT2%d1s%04d.mat', type.Permanent, round(etas(1)*100));
save(output_name, 'description', 'results', 'NUM_EVENT', 'etas', ...
    'options', 'node_opt', 'link_opt', 'VNF_opt', 'slice_opt', 'type', 'idx', 'EXPNAME');
%}

%% Test
%{
i = 5;
x = (1:height(results.Dimconfig{i}));
plot(x,ema(results.Dimconfig2{i}.Profit,0.3), x,results.Dimconfig2{i}.Profit)
%}
%{
for i = length(runtime.admm):-1:1
    mean_time(i) = mean(runtime.admm{i})/;
end
plot(penalty, mean_time);
%}