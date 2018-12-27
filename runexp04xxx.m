%% Run 
% used by <runexp04> and <runtest04>.

%% Preparing Simluation
% (1) Regularize the method names;
% (2) count the number of runs.
num_methods = length(test_methods);
title_names = cell(num_methods, 1);
TOTAL_NUM = 0;
for i = num_methods:-1:1
	invoke_methods(i) = ReconfigMethod.Enumerate(test_methods{i});
	title_names{i} = invoke_methods(i).FullName;
	if ismember(invoke_methods(i),[ReconfigMethod.Baseline, ReconfigMethod.DimBaseline])...
			&& strcmpi(mode, 'vareta')
		TOTAL_NUM = TOTAL_NUM + NUM_EVENT;
	else
		TOTAL_NUM = TOTAL_NUM + NUM_EVENT*num_vars;
	end
end
setwaitbar;

%% Run Simulation
clear user_opt;
if ~strcmpi(mode, 'vareta') && numel(etas) ~= 1
	warning('%s: ''eta'' is not a scalar.', calledby);
end
if ~strcmpi(mode, 'varweight') && numel(weight) ~= 1
	warning('%s: ''weight'' is not a scalar.', calledby);
end
if ~strcmpi(mode, 'varnumber') && numel(numberflow) ~= 1
	warning('%s: ''numberflow'' is not a scalar.', calledby);
end
for j = 1:length(invoke_methods)
	options.ReconfigMethod = invoke_methods(j);
	progress_bar.Name = horzcat(EXPNAME, ' - ', title_names{j});
	pause(0.01);
	for i = 1:num_vars
		GlobalState.Initialize();
		seed_dynamic = SEED; 
		switch mode
			case 'varweight'
				user_opt.Weight = weight(i);
				if exist('numberflow', 'var')
					user_opt.NumberFlows = numberflow(1);
				end
				if exist('etas', 'var')
					options.UnitReconfigureCost = etas(1);
				end
			case 'varnumber'
				if exist('weight', 'var')
					user_opt.Weight = weight(1);
				end
				user_opt.NumberFlows = numberflow(i);
				if exist('etas', 'var')
					options.UnitReconfigureCost = etas(1);
				end
			case 'vareta'
				if exist('weight', 'var')
					user_opt.Weight = weight(1);
				end
				if exist('numberflow', 'var')
					user_opt.NumberFlows = numberflow(1);
				end
				options.UnitReconfigureCost = etas(i);
		end
		SingleSliceReconfiguration;
		if i == 1
			results.(invoke_methods(j).char) = {g_results};
		else
			results.(invoke_methods(j).char){i,1} = g_results;
		end
		if ismember(invoke_methods(j), [ReconfigMethod.Baseline,ReconfigMethod.DimBaseline])...
				&& strcmpi(mode, 'vareta')
			break;
		end
	end
end
close(progress_bar);