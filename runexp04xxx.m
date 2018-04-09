%% Run script
% single slice reconfiguration
%%
invoke_methods = [];
title_names = cell(0);
if ~exist('mode', 'var')
    mode = 'var-eta';
end
switch mode
    case 'var-weight'
        num_vars = length(weight);
    case 'var-eta'
        num_vars = length(etas);
    case 'var-number'
        num_vars = length(numberflow);
    case 'var-penalty'
        if isempty(penalty)
            num_vars = 1;
        else
            num_vars = length(penalty);
        end
end
TOTAL_NUM = 0;
%% Fast Reconfiguration with Resource Reservation
if b_fastconfig 
    title_names = [title_names, 'Fast Reconfiguration'];
    invoke_methods = [invoke_methods; ReconfigMethod.FastconfigReserve];
    TOTAL_NUM = TOTAL_NUM + NUM_EVENT*num_vars;
end
if b_fastconfig0
    title_names = [title_names, 'Fast Reconfiguration 0'];
    invoke_methods = [invoke_methods; ReconfigMethod.Fastconfig];
    TOTAL_NUM = TOTAL_NUM + NUM_EVENT*num_vars;
end
if b_dimconfig
    title_names = [title_names, 'Hybrid Slicing Scheme'];
    invoke_methods = [invoke_methods; ReconfigMethod.DimconfigReserve];
    TOTAL_NUM = TOTAL_NUM + NUM_EVENT*num_vars;
end
if b_dimconfig0
    title_names = [title_names, 'Hybrid Slicing Scheme 0'];
    invoke_methods = [invoke_methods; ReconfigMethod.Dimconfig];
    TOTAL_NUM = TOTAL_NUM + NUM_EVENT*num_vars;
end
if b_baseline
    % Only need to compute once for different reconfiguration cost efficient, since the
    % optimization procedure is independent on the cosefficient, and the reconfiguration
    % cost with other coefficient can be derieved from the one results by using the
    % coefficient.
    title_names = [title_names, 'Baseline Reconfiguration'];
    invoke_methods = [invoke_methods; ReconfigMethod.Baseline];
    if strcmpi(mode, 'var-eta')
        TOTAL_NUM = TOTAL_NUM + NUM_EVENT;
    else
        TOTAL_NUM = TOTAL_NUM + NUM_EVENT*num_vars;
    end
end
if b_dimbaseline
    title_names = [title_names, 'Baseline Reconfiguration With Dimensioning'];
    invoke_methods = [invoke_methods; ReconfigMethod.DimBaseline];
    if strcmpi(mode, 'var-eta')
        TOTAL_NUM = TOTAL_NUM + NUM_EVENT;
    else
        TOTAL_NUM = TOTAL_NUM + NUM_EVENT*num_vars;
    end
end
%-------------------------------------------------------------------------
if b_fastconfig2
    title_names = [title_names, 'Fast Reconfiguration 2'];
    invoke_methods = [invoke_methods; ReconfigMethod.Fastconfig2];
    TOTAL_NUM = TOTAL_NUM + NUM_EVENT*num_vars;
end
if b_dimconfig2
    title_names = [title_names, 'Hybrid Slicing Scheme 2'];
    invoke_methods = [invoke_methods; ReconfigMethod.Dimconfig2];
    TOTAL_NUM = TOTAL_NUM + NUM_EVENT*num_vars;
end
%%
global total_iter_num;
total_iter_num = 0;
if exist('progress_bar', 'var') && isvalid(progress_bar)
    close(progress_bar);
end
progress_bar = waitbar(total_iter_num/TOTAL_NUM, ...
    sprintf('Simulation Progress: %d/%d', total_iter_num, TOTAL_NUM));
jframe=getJFrame(progress_bar);
jframe.setAlwaysOnTop(1);

%%
clear user_opt;
if ~strcmpi(mode, 'var-eta') && exist('etas', 'var') && ~isscalar(etas)
    warning('%s: ''eta'' is not a scalar.', calledby);
end
if ~strcmpi(mode, 'var-weight') && exist('weight', 'var') && ~isscalar(weight)
    warning('%s: ''weight'' is not a scalar.', calledby);
end
if ~strcmpi(mode, 'var-number') && exist('numberflow', 'var') && ~isscalar(numberflow)
    warning('%s: ''numberflow'' is not a scalar.', calledby);
end
if ~strcmpi(mode, 'var-penalty') && exist('penalty', 'var') && ~isscalar(numberflow)
    warning('%s: ''penalty'' is not a scalar.', calledby);
end
for j = 1:length(invoke_methods)
    options.ReconfigMethod = invoke_methods(j);
    progress_bar.Name = horzcat(EXPNAME, ' - ', title_names{j});
    pause(0.01);
    for i = 1:num_vars
        GlobalState.Initialize();
        seed_dynamic = SEED;
        switch mode
            case 'var-weight'
                user_opt.Weight = weight(i);
                if exist('numberflow', 'var')
                    user_opt.NumberFlows = numberflow(1);
                end
                if exist('etas', 'var')
                    options.UnitReconfigureCost = etas(1);
                end
            case 'var-number'
                if exist('weight', 'var')
                    user_opt.Weight = weight(1);
                end
                user_opt.NumberFlows = numberflow(i);
                if exist('etas', 'var')
                    options.UnitReconfigureCost = etas(1);
                end
            case 'var-eta'
                if exist('weight', 'var')
                    user_opt.Weight = weight(1);
                end
                if exist('numberflow', 'var')
                    user_opt.NumberFlows = numberflow(1);
                end
                options.UnitReconfigureCost = etas(i);
            case 'var-penalty'
                if ~isempty(penalty)
                    user_opt.penalty = penalty(i);
                end
                if exist('weight', 'var')
                    user_opt.Weight = weight(1);
                end
                if exist('numberflow', 'var')
                    user_opt.NumberFlows = numberflow(1);
                end
                if exist('etas', 'var')
                    options.UnitReconfigureCost = etas(1);
                end
        end
        SingleSliceReconfiguration;
        if i == 1
            results.(invoke_methods(j).char) = {g_results};
        else
            results.(invoke_methods(j).char){i,1} = g_results;
        end
        if strcmpi(mode, 'var-eta') && ...
                (invoke_methods(j)==ReconfigMethod.Baseline || ...
                invoke_methods(j)==ReconfigMethod.DimBaseline)
            break;
        end
    end
end

close(progress_bar);
