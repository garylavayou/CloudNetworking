%% Run script
% single slice reconfiguration: a warm-up phase is reused.
NUM_TEST = length(thetas);
TOTAL_NUM= options.NumberEventWarmUp*NUM_TEST + ...
    NUM_EVENT*(NUM_TEST*(b_fastconfig+b_fastconfig2+b_dimconfig+b_dimconfig2)+b_reconfig);
global total_iter_num;
total_iter_num = 0;
if exist('progress_bar', 'var') && isvalid(progress_bar)
    close(progress_bar);
end
progress_bar = waitbar(total_iter_num/TOTAL_NUM, ...
    sprintf('Simulation Progress: %d/%d', total_iter_num, TOTAL_NUM));
jframe=getJFrame(progress_bar);
jframe.setAlwaysOnTop(1);
% WindowAPI(progress_bar, 'topmost');

%% Start phase
progress_bar.Name = horzcat(EXPNAME, ' - ', 'Warm-up Phase');
pause(0.01);
PNs = creatempty(options.NetworkType, length(thetas), 0);
Dispachers = creatempty('SliceFlowEventDispatcher', length(thetas), 0);
global_state(length(thetas),1) = GlobalState;
old_num_event = NUM_EVENT;
NUM_EVENT = options.NumberEventWarmUp; %#ok<NASGU>
options.Method = 'dimconfig';
for i = 1:length(thetas)
    GlobalState.Initialize();
    seed_dynamic = SEED;
    DynamicSlice.THETA(thetas(i));
    SingleSliceReconfiguration;
    %% Create Network
    PNs(i) = PN;
    Dispachers(i) = SFED;
    global_state(i).Save;
    if i == 1
        results.Warmup = {g_results};
    else
        results.Warmup{i} = g_results;
    end
end
NUM_EVENT = old_num_event;

PNsnew = creatempty(options.NetworkType, length(thetas), 0);
Dispachersnew = creatempty('SliceFlowEventDispatcher', length(thetas), 0);
for i = 1:length(thetas)
    PNsnew(i) = PNs(i).copy;
    Dispachersnew(i) = Dispachers(i).copy;
end
delete(PNs);
delete(Dispachers);
PNs = PNsnew;
Dispachers = Dispachersnew;
%%
if b_fastconfig
    options.Method = 'fastconfig';    % {'reconfig', 'fastconfig', 'dimension', 'fastconfig2'}
    progress_bar.Name = horzcat(EXPNAME, ' - ', 'Fast Reconfiguration');
    pause(0.01);
    for i = 1:length(thetas)
        global_state(i).Restore();
        PN = PNs(i).copy;
        SFED = Dispachers(i).copy;
        DynamicSlice.THETA(thetas(i));  % should be reset each iteration.
        RepeatSliceReconfiguration;
        if i == 1
            results.Fastconfig = {g_results};
        else
            results.Fastconfig{i,1} = g_results;
        end
    end
end

%%
if b_fastconfig2
    options.Method = 'fastconfig2';    % {'reconfig', 'fastconfig', 'dimension', 'fastconfig2'}
    progress_bar.Name = horzcat(EXPNAME, ' - ', 'Fast Reconfiguration 2');
    pause(0.01);
    for i = 1:length(thetas)
        global_state(i).Restore();
        PN = PNs(i).copy;
        SFED = Dispachers(i).copy;
        DynamicSlice.THETA(thetas(i));
        RepeatSliceReconfiguration;
        if i == 1
            results.Fastconfig2 = {g_results};
        else
            results.Fastconfig2{i,1} = g_results;
        end
    end
end

%%
if b_reconfig
    % Only need to compute once for different reconfiguration cost efficient, since the
    % optimization procedure is independent on the cosefficient, and the reconfiguration
    % cost with other coefficient can be derieved from the one results by using the
    % coefficient.
    options.Method = 'reconfig';
    progress_bar.Name = horzcat(EXPNAME, ' - ', 'Reconfiguration');
    pause(0.01);
    global_state(1).Restore();
    PN = PNs(1).copy;
    SFED = Dispachers(1).copy;
    DynamicSlice.THETA(1);
    RepeatSliceReconfiguration;
    results.Reconfig = g_results;
end

%%
if b_dimconfig
    options.Method = 'dimconfig';
    progress_bar.Name = horzcat(EXPNAME, ' - ', 'Hybrid Slicing Scheme');
    pause(0.01);
    for i = 1:length(thetas)
        global_state(i).Restore();
        PN = PNs(i).copy;
        SFED = Dispachers(i).copy;
        DynamicSlice.THETA(thetas(i));
        RepeatSliceReconfiguration;
        if i == 1
            results.Dimconfig = {g_results};
        else
            results.Dimconfig{i,1} = g_results;
        end
    end
end

%%
if b_dimconfig2
    options.Method = 'dimconfig2';
    progress_bar.Name = horzcat(EXPNAME, ' - ', 'Hybrid Slicing Scheme 2');
    pause(0.01);
    for i = 1:length(thetas)
        global_state(i).Restore();
        PN = PNs(i).copy;
        SFED = Dispachers(i).copy;
        DynamicSlice.THETA(thetas(i));
        RepeatSliceReconfiguration;
        if i == 1
            results.Dimconfig2 = {g_results};
        else
            results.Dimconfig2{i,1} = g_results;
        end
    end
end

close(progress_bar);
