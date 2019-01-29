%% Fast Slice Reconfiguration and the Parallel Computing Algorithm
% Task 3: evaluating the efficacy of the parallel computing algorithm
%		(dual-ADMM) for FSR, in terms of running time and the approximity.
%
%	NOTE: the iteration limit of the dual-ADMM-EB method is derived from the
%	results of running dual-ADMM with varying number of flows.
%
%	ISNTRUCTION:
%		To perform the experiment in <On Fast Slice Reconfiguration>, we need
%		to run the method three times with different algorithms.
% SEE ALSO:
%		Use <runexp04> to evaluating the influcence of parameters to the
%		reconfiguration ratio and the performance of FSR.
%% Tasks
% task = {'varnumber', 'varpenalty', 'test'}
% algorithm = {'normal', 'admm', 'admm_limit'}
function [results, runtime] = runexp05(task, algorithm, test_params)
%%
pre_slice_reconfigure;
global computime NUMBER_ITERS;
%%
if nargin == 0
	i = input(strcat('Please give the task:', newline,...
		'1. ''varnumber'', 2. ''varpenalty'', 3.''test'':{1|2|3}',...
		newline), 's');
	switch i
		case '1', task = 'varnumber';
		case '2', task = 'varpenalty';
		case '3', task = 'test';
		otherwise
			task = 'varnumber';
			warning('invalid value, task is set to ''%s''.', task);
	end
end
if nargin <= 1
	if strcmpi(task, 'varnumber')
		i = input(strcat('Please specify the algorithm:', newline,...
			'1. ''normal'', 2. ''dual-ADMM'', 3.''dual-ADMM-EB'':{1|2|3}',...
			newline), 's');
		switch i
			case '1', algorithm = 'normal';
			case '2', algorithm = 'admm';
			case '3'
				algorithm = 'admm_limit';
				predefined_iter_limit = [93 77 58 35 28 22 28 27 27 23];   % average
			otherwise
				algorithm = 'normal';
				warning('invalid value, algorithm is set to ''%s''.', algorithm);
		end
	elseif strcmpi(task, 'varpenalty')
		algorithm = 'admm';
	end
end

switch task
	case 'varnumber'
		%% Parameter Specification
		% The parameter specification conforms to the paper <On Fast Slice
		% Reconfiguration>.
		etas = 1; numberflow = 100:100:1000; weight = 10; penalty = 4;
		num_vars = length(numberflow);
	case 'varpenalty'
		etas = 1; numberflow = 1000; weight = 10;
		penalty = [1;2;4;8;12;16;20];   % (24,32) is too large to obtain accurate results.
		num_vars = length(penalty);
	case 'test'
		etas = 2; numberflow = 1000; weight = 10; penalty = 0.5;  % [] for normal.
		num_vars  = 1;
		if nargin >= 2
			if isfield(test_params, 'Weight')
				weight = test_params.Weight;
			end
			if isfield(test_params, 'Penalty')
				penalty = test_params.Penalty;
			end
			if isfield(test_params, 'Etas')
				etas = test_params.Etas;
			end
			if isfield(test_params, 'NumberFlows')
				numberflow = test_params.NumberFlows;
			end
		end
end
if strcmpi(algorithm, 'normal')
	penalty = [];
end
test_methods = {'Fastconfig'};   % FSR
NUM_EVENT = 11;
EXPNAME = sprintf('EXP05');
b_save = true;

%% Run script
if ~strcmpi(algorithm, 'normal')
	computime = zeros(NUM_EVENT-1,1);
	NUMBER_ITERS = zeros(NUM_EVENT-1,1);
else
	clear computime;
	clear NUMBER_ITERS;
end
if exist('penalty', 'var') && ...
		~isempty(penalty) && license('test', 'Distrib_Computing_Toolbox')
	p = gcp;
	if isempty(p)
		warning('Parallel computing is disabled.');
	else
		fprintf('Starting up parallel computing pool, number of workers: %d.\n', p.NumWorkers);
	end
end
options.ReconfigMethod = ReconfigMethod.Enumerate(test_methods{1});
title_name = options.ReconfigMethod.FullName;
TOTAL_NUM = NUM_EVENT*num_vars; %#ok<NASGU> used to record progress
setwaitbar;
progress_bar.Name = horzcat(EXPNAME, ' - ', title_name);
pause(0.01);
%%
global ITER_LIMIT;
clear user_opt;
if ~strcmpi(task, 'varnumber') && numel(numberflow) ~= 1
	warning('%s: ''numberflow'' is not a scalar.', calledby);
end
if strcmpi(task, 'varnumber') && numel(penalty) ~= 1
	warning('[%s] ''penalty'' is not a scalar.', calledby);
end

for i = 1:num_vars
	GlobalState.Initialize();
	seed_dynamic = SEED;			%#ok<NASGU> used in each single experiment
	switch task
		case 'varnumber'
			user_opt.Weight = weight(1);
			user_opt.NumberFlows = numberflow(i);
			options.UnitReconfigureCost = etas(1);
			if ~isempty(penalty)
				user_opt.IntraSlicePenalty = penalty(1);
			end
		case 'varpenalty'
			user_opt.IntraSlicePenalty = penalty(i);
			user_opt.Weight = weight(1);
			user_opt.NumberFlows = numberflow(1);
			options.UnitReconfigureCost = etas(1);
	end
	if exist('predefined_iter_limit', 'var')
		if isscalar(predefined_iter_limit)
			ITER_LIMIT = predefined_iter_limit;
		else
			ITER_LIMIT = predefined_iter_limit(i);
		end
	else
		ITER_LIMIT = inf;
	end
	SingleSliceReconfiguration;			% return g_results
	if i == 1
		results= {g_results};
	else
		results{i,1} = g_results;
	end
	if exist('computime', 'var') && ~isempty(computime)
		if ~exist('runtime', 'var') || ~isfield(runtime, task)
			runtime.(task).(algorithm) = computime;
		else
			runtime.(task)(i).(algorithm) = computime;
		end
		if contains(algorithm, 'admm')
			if ~exist('numiters', 'var') || ~isfield(numiters, task)
				numiters.(task).(algorithm) = NUMBER_ITERS;
			else
				numiters.(task)(i).(algorithm) = NUMBER_ITERS;
			end
		end
	end
end
close(progress_bar);

%% Output

%% Test
%{
i = 5;
x = (1:height(results.Dimconfig{i}));
plot(x,ema(results.Dimconfig2{i}.Profit,0.3), x,results.Dimconfig2{i}.Profit)
%}

%% var size
% plot admm runtime
%{
for i = length(runtime.varsize):-1:1
    mean_time(i) = mean(runtime.varsize(i).admm)/numberflow(i)*6;
end
plot(numberflow,mean_time);
%}
% plot normal runtime
%{
for i = length(runtime.varsize):-1:1
    mean_time(i) = mean(runtime.varsize(i).normal);
end
plot(numberflow,mean_time);
%}
% statistics on number of iterations
%{
ll = length(numiters.varsize);
max_iters = zeros(1,ll);
min_iters = zeros(1,ll);
avg_iters = zeros(1,ll);
for i = 1:length(numiters.varsize)
    max_iters(i) = max(numiters.varsize(i).admm);
    min_iters(i) = min(numiters.varsize(i).admm);
    avg_iters(i) = round(mean(numiters.varsize(i).admm));
end
disp(max_iters)
disp(min_iters)
disp(avg_iters)
% 139   130   148   106    67    65    53    45    38    34
%  62    44    27    15    16    13    17    14    17    14
%  93    77    58    35    28    22    28    27    27    23
%}
% save results
if b_save && strcmpi(task, 'varnumber')
	description = [sprintf('Experiment 5-1: running time of %s method, ', algorithm),...
		'varying slice size (number of flows)']; %#ok<NASGU>
	if strcmpi(algorithm, 'normal')
		output_name = sprintf('Results/singles/%s1e%dw%dp_.mat', EXPNAME, etas, weight);
	else
		output_name = sprintf('Results/singles/%s1e%dw%dp%d.mat', EXPNAME, etas, weight, penalty);
	end
	save(output_name, 'description', 'results', 'etas', 'weight', 'penalty', ...
		'numberflow', 'runtime', 'numiters');
end

if b_plot
	dataplot3rtpf;
end

%% var penalty
% plot
%{
for i = length(runtime.varpenalty):-1:1
    mean_time(i) = mean(runtime.varpenalty(i).admm);
end
plot(penalty,mean_time);
%}
% save results
if b_save && strcmpi(task, 'varpenalty')
	description = 'Experiment 5-2: running time of dual ADMM, varying penalty factor'; %#ok<NASGU>
	output_name = sprintf('Results/singles/EXP0502e%dn%dw%d.mat', etas, numberflow, weight);
	save(output_name, 'description', 'results', 'etas', 'weight', 'numberflow', 'penalty', ...
		'runtime', 'numiters');
end
% NOTE: penalty factor should not be too large, _i.e._, $r<=2$ will be fine; otherwise the
%   results will be inaccurate. On the other hand the penalty should not be too small,
%   otherwise the number of iteration will be large and the convergence rate is slow.
