%% Run script
% single slice reconfiguration
NUM_TEST = length(etas);
TOTAL_NUM= NUM_EVENT*(NUM_TEST*(...
    b_fastconfig+b_fastconfig2+b_dimconfig+b_dimconfig1+b_dimconfig2+b_dimconfig0)+b_reconfig);
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
if b_reconfig
    % Only need to compute once for different reconfiguration cost efficient, since the
    % optimization procedure is independent on the cosefficient, and the reconfiguration
    % cost with other coefficient can be derieved from the one results by using the
    % coefficient.
    options.ReconfigMethod = 'reconfig';
    progress_bar.Name = horzcat(EXPNAME, ' - ', 'Reconfiguration');
    pause(0.01);
    GlobalState.Initialize();
    seed_dynamic = SEED;
    options.UnitReconfigureCost = 1;
    SingleSliceReconfiguration;
    results.Reconfig = g_results;
end

%%
if b_fastconfig
    options.ReconfigMethod = 'fastconfig';    % {'reconfig', 'fastconfig', 'dimension', 'fastconfig2'}
    progress_bar.Name = horzcat(EXPNAME, ' - ', 'Fast Reconfiguration');
    pause(0.01);
    for i = 1:length(etas)
        GlobalState.Initialize();
        seed_dynamic = SEED;
        options.UnitReconfigureCost = etas(i);
        SingleSliceReconfiguration;
        if i == 1
            results.Fastconfig = {g_results};
        else
            results.Fastconfig{i,1} = g_results;
        end
    end
end

%%
if b_fastconfig2
    options.ReconfigMethod = 'fastconfig2';    % {'reconfig', 'fastconfig', 'dimension', 'fastconfig2'}
    progress_bar.Name = horzcat(EXPNAME, ' - ', 'Fast Reconfiguration 2');
    pause(0.01);
    for i = 1:length(etas)
        GlobalState.Initialize();
        seed_dynamic = SEED;
        options.UnitReconfigureCost = etas(i);
        SingleSliceReconfiguration;
        if i == 1
            results.Fastconfig2 = {g_results};
        else
            results.Fastconfig2{i,1} = g_results;
        end
    end
end

%%
if b_dimconfig
    options.ReconfigMethod = 'dimconfig';
    progress_bar.Name = horzcat(EXPNAME, ' - ', 'Hybrid Slicing Scheme');
    pause(0.01);
    for i = 1:length(etas)
        GlobalState.Initialize();
        seed_dynamic = SEED;
        options.UnitReconfigureCost = etas(i);
        SingleSliceReconfiguration;
        if i == 1
            results.Dimconfig = {g_results};
        else
            results.Dimconfig{i,1} = g_results;
        end
    end
end

%%
if b_dimconfig1
    options.ReconfigMethod = 'dimconfig1';
    progress_bar.Name = horzcat(EXPNAME, ' - ', 'Hybrid Slicing Scheme 1');
    pause(0.01);
    for i = 1:length(etas)
        GlobalState.Initialize();
        seed_dynamic = SEED;
        options.UnitReconfigureCost = etas(i);
        SingleSliceReconfiguration;
        if i == 1
            results.Dimconfig1 = {g_results};
        else
            results.Dimconfig1{i,1} = g_results;
        end
    end
end

%%
if b_dimconfig1
    options.ReconfigMethod = 'dimconfig2';
    progress_bar.Name = horzcat(EXPNAME, ' - ', 'Hybrid Slicing Scheme 2');
    pause(0.01);
    for i = 1:length(etas)
        GlobalState.Initialize();
        seed_dynamic = SEED;
        options.UnitReconfigureCost = etas(i);
        SingleSliceReconfiguration;
        if i == 1
            results.Dimconfig2 = {g_results};
        else
            results.Dimconfig2{i,1} = g_results;
        end
    end
end

%%
if b_dimconfig0
    options.ReconfigMethod = 'dimconfig0';
    progress_bar.Name = horzcat(EXPNAME, ' - ', 'Hybrid Slicing Scheme 0');
    pause(0.01);
    for i = 1:length(etas)
        GlobalState.Initialize();
        seed_dynamic = SEED;
        options.UnitReconfigureCost = etas(i);
        SingleSliceReconfiguration;
        if i == 1
            results.Dimconfig0 = {g_results};
        else
            results.Dimconfig0{i,1} = g_results;
        end
    end
end

close(progress_bar);
