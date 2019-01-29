%% Hybrid Slice Reconfiguration
% Task 11: evaluating the influcence of serveral parameters to the
%		reconfiguration ratio of HSR and HSR-RSV (FSR), compared with FSR and a
%		baseline scheme that does not consider reconfiguration cost.
%	Task 12: evaluting the performance of HSR and HSR-RSV.
% Task 21: evaluating the influcence of parameters to the reconfiguration
%		ratio of FSR. 
%	Task 22: evaluting the performance of FSR.
%
%% Tasks
% task = {'perfeval.hsr', 'testparam.hsr', 'perfeval.fsr', 'testparam.fsr'}
% mode = {'vareta', 'varweight', 'varnumber'}
function results = runexp04(task, mode, ioopts)
global DEBUG;
DEBUG = true;
EXPNAME = sprintf('EXP041');
%% Predefined experiment settings
%		> Network Type:
%			- DynamicCloudNetwork(1)
%		> Topology:
%			- Sample2(2)
%		> Warm-up phase:
%			- No(0)
pre_slice_reconfigure;
if nargin >= 3
	b_save = bitand(ioopts, hex2dec('0001')); % initialized in pre_slice_reconfigure
	b_plot = bitand(ioopts, hex2dec('0002'));
end
%% Handle User Input for Configuration
if nargin == 0
	i = input(['Please give the task: ', newline, ...
		'11. ''perfeval.HSR'', 12. ''testparam.HSR'':', newline, ...
		'21. ''perfeval.FSR'', 22. ''testparam.FSR'':', newline, ...
		'{11|12|21|22}'], 's');
	switch i
		case {'11', ''}, task = 'perfeval.hsr';
		case '12', task = 'testparam.hsr';
		case '21', task = 'perfeval.fsr';
		case '22', task = 'testparam.fsr';
		otherwise
			task = 'perfeval.hsr';
			warning('invalid value, task is set to ''%s''.', task);
	end
end
if nargin <= 1 && startsWith(task,'testparam')
	i = input('Please give the mode: 1. ''vareta'', 2. ''varweight'', 3.''varnumber'':{1|2|3}', 's');
	switch i
		case '1', mode = 'vareta';
		case '2', mode = 'varweight';
		case '3', mode = 'varnumber';
		otherwise
			mode = 'vareta';
			warning('invalid value, mode is set to ''%s''.', mode);
	end
end
%% Set Simulation Parameters
if endsWith(task, 'hsr', 'IgnoreCase', true)
	NUM_EVENT = 401;            % the trigger-interval is set to 60.
	test_methods = {'DimBaseline', 'Dimconfig', 'DimconfigReserve'}; %#ok<*NASGU>
elseif endsWith(task, 'fsr', 'IgnoreCase', true)
	NUM_EVENT = 51;
	test_methods = {'Baseline', 'Fastconfig'};
end
if startsWith(task,'perfeval', 'IgnoreCase', true)
	mode = 'vareta'; etas = 1; numberflow = 100; weight = 10; num_vars = 1;
elseif startsWith(task,'testparam', 'IgnoreCase', true)
	switch mode
		case 'vareta'
			etas = [1/32 1/16 1/8 1/4 1/2 1 2 4 8]; 
			numberflow = 100; weight = 10; num_vars = length(etas);
		case 'varweight'
			weight = 10:10:80; 
			etas = 1;  numberflow = 100; num_vars = length(weight);
		case 'varnumber'
			numberflow = 30:30:240; 
			etas = 1; weight = 10; num_vars = length(numberflow);
		otherwise
			error('error: invalid mode [%s].', mode);
	end
end

%% Run Simulation
runexp04xxx;			% see the script.

%% Output
%% plot figure
if b_plot
	%%
	if ~exist('mode', 'var')
		i = input('Please specify mode: 1. ''vareta'', 2. ''var-weight'', 3. ''varnumber'' (1|2|3)?', 's');
		switch i
			case {'1', 'vareta', ''}, mode = 'vareta';
			case {'2', 'varweight'}, mode = 'varweight';
			case {'3', 'varnumber'}, mode = 'varnumber';
			otherwise
				mode = 'vareta';
				warning('invalid value, mode is set to ''%s''.', mode);
		end
	end
	%%
	options.bSavePlot = b_plotsave;
	if endsWith(task, 'hsr', 'IgnoreCase',true)
		dataplot3l1approx
		%%
		load('Results/singles/EXP41_OUTPUT241s0100.mat', 'results');
		dataplot3sa2(results); %#ok<*NODEF>
		%%
		dataplot3params(mode, [], [], [], options)
	else
		%%
		lines = struct('Sources', {{'Fastconfig', 'DimBaseline'}},...
			'Labels', {{'FSR', 'Baseline'}});
		options.Suffix = '-fastcomp';
		%%
		dataplot3params(mode, lines, [], [], options)
		%%
		if ~exist('results', 'var')
			load('Results/singles/EXP41_OUTPUT241s0100.mat', 'results');
		end
		for i = 1:length(lines.Sources)
			results.(lines.Sources{i}) = results.(lines.Sources{i})(ex_id);
		end
		options.Suffix = '-fast';
		dataplot3sa2(results, lines, [], 1, options);
	end
end
%% Save Results
if b_save
	if endsWith(task, 'hsr', 'IgnoreCase', true)
		tag = 'HSR and HSR-RSV';
		suffix = '';
	elseif endsWith(task, 'fsr', 'IgnoreCase', true)
		tag = 'FSR';
		suffix = '_fast';
	end
	if startsWith(task, 'testparam', 'IgnoreCase', true)
		description = 'Experiment 4-11: verify the influence of network settings to ';
		output_name = sprintf('Results/%s1_%s%s.mat', EXPNAME, mode, suffix);
	elseif startsWith(task, 'perfeval', 'IgnoreCase', true)
		description = 'Experiment 4-12: Performance evaluation of ';
		output_name = sprintf('Results/%s2_%se%04d.mat', EXPNAME, suffix, round(etas(1)*100));
	end
	description = [description, tag, ...
			'(un-ordered SFC, without warm-up phase, disable ad-hoc mode).', newline, ...
			'Topology = ', node_opt.Model.char, newline, ...
			sprintf('SliceType = %d\n', type.Index(type.Permanent))];
	if startsWith(task, 'testparam', 'IgnoreCase', true)
		description = [description, 'variables = ' mode(4:end), '.']; 
	elseif  startsWith(task, 'perfeval', 'IgnoreCase', true)
		description = [description, '.']; 
	end
	save(output_name, 'description', 'results', ...
		'mode', 'etas', 'numberflow', 'weight', 'options', 'node_opt', ...
		'link_opt', 'VNF_opt', 'slice_opt', 'type', 'NUM_EVENT');
end
