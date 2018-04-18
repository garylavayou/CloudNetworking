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
% mode = 'var-penalty'; etas = 1; numberflow = 1000; weight = 30; penalty = [1;2;4;8;12;16;20];   % (24,32) is too large to obtain accurate results.
% mode = 'var-number'; etas = 1; numberflow = 100:100:1000; weight = 30; penalty = 2;   % penalty = [] for normal method.
b_fastconfig0 = true;        % FSR
NUM_EVENT = 11;           
idx = 1:NUM_EVENT;
if strcmpi(mode, 'var-penalty') || strcmpi(mode, 'var-number')
    computime = zeros(NUM_EVENT-1,1);
else
    clear computime;
end
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

%% Test
%{
i = 5;
x = (1:height(results.Dimconfig{i}));
plot(x,ema(results.Dimconfig2{i}.Profit,0.3), x,results.Dimconfig2{i}.Profit)
%}
%{
for i = length(runtime.admm):-1:1
    mean_time(i) = mean(runtime.admm{i});
end
plot(penalty, mean_time);
%}
%% var penalty
% plot
%{
for i = length(runtime.varpenalty):-1:1
    mean_time(i) = mean(runtime.varpenalty(i).admm);
end
plot(penalty,mean_time);
%}
% save results
%{
description = ['Experiment 6001: running time of dual ADMM and the normal method',...
    'varying penalty factor'];
output_name = 'Results/singles/EXP6001.mat';
save(output_name, 'description', 'results', 'runtime', 'penalty');
%}
% NOTE: penalty factor should not be too large, _i.e._, $r<=2$ will be fine; otherwise the
%   results will be inaccurate. On the other hand the penalty should not be too small,
%   otherwise the number of iteration will be large and the convergence rate is slow.

%% var size
% plot
%{
for i = length(runtime.varsize):-1:1
    mean_time(i) = mean(runtime.varsize(i).admm);
end
plot(numberflow,mean_time);
%}
% save results
%{
description = ['Experiment 6002: running time of dual ADMM and the normal method',...
    'varying slice size (number of flows)'];
output_name = 'Results/singles/EXP6002.mat';
save(output_name, 'description', 'results', 'runtime', 'penalty', 'numberflow');
%}