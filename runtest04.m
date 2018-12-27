%% Hybrid Slice Reconfiguration
% Task 1: evaluating the influcence of serveral parameters to the
%		reconfiguration ratio of HSR and HSR-RSV, compared with FSR and a
%		baseline scheme that does not consider reconfiguration cost. (testparams)
%	Task 2: evaluting the performance of HSR and HSR-RSV. (perfeval)
%
% Task 3: evaluating the efficacy of the parallel computing algorithm
%		(dual-ADMM) for FSR. 
global DEBUG;
DEBUG = true;

pre_slice_reconfigure;
VNF_opt.Ordered = true;

task = 'perfeval';			% {'perfeval'|'testparam'|'paracomp'}
mode = 'vareta';
% b_save = false; % initialized in pre_slice_reconfigure
% b_plot = false;

% type.Index = [144; 154; 164; 174; 184];
% type.Permanent = 4;
% type.Static = [1; 2; 3];
type.StaticCount = [1; 2; 2];
% type.StaticClass = {'SimpleSlice'};
if startsWith(task, 'perfeval')
	mode = 'vareta'; etas = 1; numberflow = 100; weight = 10; num_vars = 1;
elseif startsWith(task, 'testparams')
	switch mode
		case 'vareta'
			etas = [1/32 1/16 1/8 1/4 1/2 1 2 4 8]; 
			numberflow = 40; weight = 10; num_vars = length(etas);
		case 'varweight'
			weight = 10:10:80; 
			etas = 1;  numberflow = 100; num_vars = length(weight);
		case 'varnumber'
			numberflow = 30:30:240; 
			etas = 1; weight = 10; num_vars = length(numberflow);
		otherwise
			error('error: invalid mode [%s].', task);
	end
elseif strcmpi(task, 'paracomp')
	error('error: to implement.');
end
if VNF_opt.Ordered
	options.NetworkType = 'NormalDynamicNetwork';
	type.StaticClass = {'NormalSlice'};
end

% test_methods = {'Baseline', 'Dimconfig', 'DimconfigReserve', 'Fastconfig', 'DimBaseline'};
test_methods = {'DimconfigReserve'};
NUM_EVENT = 401;            % {200,400} the trigger-interval is set to 50. 
EXPNAME = sprintf('TEST04');
warning('off', 'backtrace');
warning('on', 'verbose');
runexp04xxx;

%% Output
% # plot figure
if b_plot
	data_plot3;
	data_plot31;
	dataplot3s;
	dataplot3sa;
end
% # Save Results
if b_save
	if exist('test_mode', 'var')
		description = sprintf('Experiment 4-1: verify the influence of %s',...
			config_num, ...
			'network settings to HSR (without warm-up phase, disable ad-hoc mode).');
		varname = split(mode, '-'); varname = varname{2};
		description = sprintf('%s\nTopology=%s\nSliceType = %d\nvariables = %s.', ...
			description, node_opt.Model.char, type.Index(type.Permanent), ...
			varname);
		output_name = sprintf('Results/EXP0401_var%s.mat', varname);
	else
		description = sprintf('Experiment 4-2: %s%s',...
			'Performance evaluation of FSR and HSR (%s).', ...
			'without warm-up phase, disable ad-hoc mode');
		description = sprintf('%s\nTopology=%s\nSliceType = %d.', ...
			description, node_opt.Model.char, type.Index(type.Permanent));
		output_name = sprintf('Results/EXP0402e%04d.mat', round(etas(1)*100));
	end
	save(output_name, 'description', 'results', ...
		'mode', 'etas', 'numberflow', 'weight', 'options', 'node_opt', ...
		'link_opt', 'VNF_opt', 'slice_opt', 'type', 'NUM_EVENT');
end
%% Test
