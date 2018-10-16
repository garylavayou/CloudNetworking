%% Hybrid Slice Reconfiguration
% Task 1: evaluating the influcence of serveral parameters to the
%		reconfiguration ratio of HSR and HSR-RSV, compared with FSR and a
%		baseline scheme that does not consider reconfiguration cost.
%	Task 2: evaluting the performance of HSR and HSR-RSV.
%
% Task 3: evaluating the efficacy of the parallel computing algorithm
%		(dual-ADMM) for FSR.

global DEBUG;
DEBUG = true;
if nargin == 0
	task = 'PerfEval';
end
if nargin >= 3
	b_save = bitand(operations, hex2dec('0001'));
	b_plot = bitand(operations, hex2dec('0002'));
else
	b_save = 0;
	b_plot = 0;
end

pre_slice_reconfigure;
type.Index = [144; 154; 164; 174; 184];
type.Permanent = 4;
type.Static = [1; 2; 3];
type.StaticCount = [1; 2; 2];
type.StaticClass = {'SimpleSlice'};
if strcmpi(task, 'PerfEval')
	mode = 'var-eta'; etas = 1; numberflow = 100; weight = 10;
elseif strcmpi(task, 'Vars')
	switch mode
		case 'var-eta'
			etas = [1/32 1/16 1/8 1/4 1/2 1 2 4 8]; numberflow = 100; weight = 10; %#ok<*NASGU>
		case 'var-weight'
			weight = 10:10:80; etas = 1;  numberflow = 100;
		case 'var-number'
			numberflow = 30:30:240; etas = 1; weight = 10;
		otherwise
			error('error: invalid mode [%s].', task);
	end
elseif strcmpi(task, 'PerfEvalFSR')
	
end
if exist('test_mode', 'var') && ~isempty(task)
	switch task
		case 'var-eta'
			etas = [1/32 1/16 1/8 1/4 1/2 1 2 4 8]; numberflow = 100; weight = 10; %#ok<*NASGU>
		case 'var-weight'
			weight = 10:10:80; etas = 1;  numberflow = 100;
		case 'var-number'
			numberflow = 30:30:240; etas = 1; weight = 10;
		otherwise
			error('error: invalid mode [%s].', task);
	end
	mode = task;
else
	mode = 'var-eta'; etas = 1; numberflow = 100; weight = 10;
end
b_dimbaseline = true;       % Baseline
b_dimconfig0 = true;        % HSR
b_dimconfig = true;         % HSR-RSV
b_fastconfig0 = true;				% FSR
NUM_EVENT = 401;            % {200,400} the trigger-interval is set to 50. 
EXPNAME = sprintf('EXP04');
runexp04xxx;

%% Output
% # plot figure
if nargin >= 3 && ~isempty(b_plot) && b_plot
	data_plot3;
	data_plot31;
	dataplot3s;
	dataplot3sa;
end
% # Save Results
if nargin>= 2 && ~isempty(b_save) && b_save
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
