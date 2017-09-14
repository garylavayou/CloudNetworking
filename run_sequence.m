%% Process a sequence of network slice requests

if (~exist('output_results', 'var') || ~isstruct(output_results))
    if exist('output_file_name', 'file')
        load('output_file_name', 'output_results');
    else
        output_results = struct;
    end
end

%% Set progress bar
current_length = 0;
if exist('progress_bar', 'var') && isvalid(progress_bar)
    close(progress_bar);
end
progress_bar = waitbar(current_length/experiment_length, ...
    sprintf('%d/%d',current_length, experiment_length));
progress_bar.Name = 'Computing Progress';
pause(0.01);


%% Construct Network
% Initialize substrate network
% if b_single_optimal || b_price_adjust1 || b_price_adjust2 || b_dual_decomp || ...
%         b_resource_part || b_part_price
if isfield(net_opt, 'NetworkClass') && strcmp(net_opt.NetworkClass, 'CloudNetwork')
    PN = CloudNetwork(node_opt, link_opt, VNF_opt, net_opt);
else
    PN = PhysicalNetwork(node_opt, link_opt, VNF_opt, net_opt);
end
PN.slice_template = Slice.loadSliceTemplate(type.Index);

if b_static_slice
    if isfield(net_opt, 'NetworkClass') && strcmp(net_opt.NetworkClass, 'CloudNetwork')
        PN_static = CloudNetwork(node_opt, link_opt, VNF_opt, net_opt);
    else
        PN_static = PhysicalNetwork(node_opt, link_opt, VNF_opt, net_opt);
    end
    PN_static.slice_template = Slice.loadSliceTemplate(type.Index);
    link_price = PN_static.getLinkCost * (1 + options.PricingFactor);
    node_price = PN_static.getNodeCost * (1 + options.PricingFactor);
    PN_static.setLinkField('Price', link_price);
    PN_static.setDataCenterField('Price', node_price);
    static_opts = options;
    static_opts.Method = 'slice-price';
    clear link_uc node_uc link_price node_price;
end

%% Network Slice as a Poisson Arrival Process
% Two types of network slice as a compound Poisson process.
if isfield(arrival, 'StopTime')
    num_events = arrival.StopTime;
else
    num_events = arrival.Number*2;
end
num_type = length(type.Index);
num_fix_type = length(type.Fixed);
event_set = struct;
e = 1;
for i = (1+num_fix_type):num_type
    event_set(e).ServiceInterval = PN.slice_template(i).ServiceInterval;
    event_set(e).Probability = PN.slice_template(i).Probability;
    e = e + 1;
end
RE = RequestEvent(event_set, arrival, seed_dynamic);
clear e i event_set;

seed_static = seed_dynamic;
for t = 1:num_fix_type
    for s = 1:type.FixedCount(t)
        slice_opt = PN.slice_template(type.Fixed(t));
        slice_opt.RandomSeed = seed_dynamic;
        seed_dynamic = seed_dynamic + 1;
        PN.AddSlice(slice_opt);
    end
end
clear s slice_opt;

%% Statistics Initialization
if b_single_optimal
    [stat_optimal, slice_stat_optimal] = ...
        createStatTable(num_events, num_type, 'optimal-spp');
    if ~isfield(net_opt, 'AdmitPolicy')
        net_opt.AdmitPolicy = 'reject-flow';
    end
end
if b_dual_decomp
    [stat_dual, slice_stat_dual] = ...
        createStatTable(num_events, num_type, 'dynamic-price');
end
if b_price_adjust1
    [stat_price1, slice_stat_price1] = ...
        createStatTable(num_events, num_type, 'dynamic-price');
end
if b_price_adjust2
    [stat_price2, slice_stat_price2] = ...
        createStatTable(num_events, num_type, 'dynamic-price');
end
if b_resource_part
    [stat_part, slice_stat_part] = ...
        createStatTable(num_events, num_type, 'dynamic-price');
end
if b_part_price
    [stat_partprice, slice_stat_partprice] = ...
        createStatTable(num_events, num_type, 'dynamic-price');
end
if b_static_slice
    [stat_static, slice_stat_static] = ...
        createStatTable(num_events, num_type, 'static');
    % TODO: the persistant slices can be optimized.
    for t = 1:num_fix_type
        for s = 1:type.FixedCount(t)
            slice_opt = PN_static.slice_template(type.Fixed(t));
            slice_opt.RandomSeed = seed_static;
            seed_static = seed_static + 1;
            slice_opt.method = 'static-slicing';
            slice_opt.admit_ploicy = net_opt.AdmitPolicy;
            sl = PN_static.AddSlice(slice_opt);
            if ~isempty(sl)
                PN_static.staticSlicing(sl, static_opts);
            end
        end
    end
end
clear t s slice_opt;
number_slices = zeros(num_events,num_type);
stat.times = zeros(num_events, 1);
%%
RE.reset;
for i = 1:num_events
    e = RE.nextEvent;
    %%% Count of slice types
    % NOTE: countCurrentType counts the type of events. If no slice is rejected, this
    % information equals to the number of slices. Otherwise, it cannot be used as the
    % count of arrival slices.
    T = RE.countCurrentType;
    disp(T);
    stat.times(i) = e.Time;
    s_type = e.Type+num_fix_type;
    if strcmp(e.Description, 'arrival')
        %%%
        % |e.Type+num_fix_type | is the index of slice type in the array.
        % |e.Type| is the index of event type.
        fprintf('\tSlice type %d.\n', s_type);
        %         if b_single_optimal || b_price_adjust1 || b_price_adjust2 || b_dual_decomp || ...
        %                 b_resource_part || b_part_price
        slice_opt = PN.slice_template(s_type);
        slice_opt.RandomSeed = seed_dynamic;
        seed_dynamic = seed_dynamic + 1;
        slice_opt.Identifier = e.Identifier;
        PN.AddSlice(slice_opt);
        %         end
        number_slices(i,:) = PN.CountSlices;
        
        if b_static_slice
            slice_opt = PN_static.slice_template(s_type);
            slice_opt.RandomSeed = seed_static;
            seed_static = seed_static + 1;
            slice_opt.Identifier = e.Identifier;
            slice_opt.method = 'static-slicing';
            slice_opt.admit_ploicy = net_opt.AdmitPolicy;
            sl = PN_static.AddSlice(slice_opt);
        end
    else
        fprintf('\tDeparting slice type %d.\n', e.Type+num_fix_type);
        %         if b_single_optimal || b_price_adjust1 || b_price_adjust2 || b_dual_decomp || ...
        %                 b_resource_part || b_part_price
        PN.RemoveSlice(e.Id);
        %         end
        number_slices(i,:) = PN.CountSlices;
        if b_static_slice
            sl = PN_static.RemoveSlice(e.Id);
        end
    end
    if PN.NumberSlices == 0         %% might not invalid, since PN may not exist
        break;
    end
    if b_static_slice
        if strcmp(e.Description, 'arrival')
            if i > 1
                stat_static{i, 'NumberReject'} = stat_static{i-1, 'NumberReject'};
                stat_static{i, 'NumberPartialReject'}...
                    = stat_static{i-1, 'NumberPartialReject'};
            end
            if isempty(sl)
                stat_static{i, 'NumberReject'}(s_type)...
                    = stat_static{i, 'NumberReject'}(s_type) + 1;
                stat_static{i, 'NumberPartialReject'}(s_type)...
                    = stat_static{i, 'NumberPartialReject'}(s_type) + 1;
            else
                if sl.NumberFlows < PN_static.slice_template(s_type).NumberFlows
                    stat_static{i, 'NumberPartialReject'}(s_type)...
                        = stat_static{i, 'NumberPartialReject'}(s_type) + 1;
                end
            end
            
            tic;
            [output_static] = PN_static.staticSlicing(sl, static_opts);
            rt = toc;
            [tb, stbs] = saveStatTable(PN_static, output_static, rt, type.Index, ...
                'static');
            stat_static(i, tb.Properties.VariableNames) = tb;
            for j = 1:num_type
                slice_stat_static{j}(i, stbs.Properties.VariableNames) = stbs(j,:);
            end
        else
            if isempty(sl)
                % Nothing should be changed.
                cprintf('comments', 'Info: slice %d have been rejected earlier.\n', e.Id);
                rt = 0;
                if i>1
                    stat_static(i,:) = stat_static(i-1,:);
                    for j = 1:num_type
                        slice_stat_static{j}(i, :) = slice_stat_static{j}(i-1, :);
                    end
                end
            else
                % When a slice is removed, staticSlicing update the network state.
                tic;
                [output_static] = PN_static.staticSlicing([], static_opts);
                rt = toc;
                %%%
                % It is determined that i > 1.
                [tb, stbs] = saveStatTable(PN_static, output_static, rt, ...
                    type.Index, 'static');
                stat_static(i, tb.Properties.VariableNames) = tb;
                for j = 1:num_type
                    slice_stat_static{j}(i, stbs.Properties.VariableNames) = stbs(j,:);
                end
            end
        end
        current_length = current_length + 1;
        waitbar(current_length/experiment_length, progress_bar, ...
            sprintf('%d/%d',current_length, experiment_length));
    end
    
    if b_single_optimal
        options.Method = 'normal';
        [output_optimal, rt] = PN.singleSliceOptimization(options);
        [tb, stbs] = saveStatTable(PN, output_optimal, rt, type.Index, ...
            'optimal-spp');
        stat_optimal(i, tb.Properties.VariableNames) = tb;
        for j = 1:num_type
            slice_stat_optimal{j}(i, stbs.Properties.VariableNames) = stbs(j,:);
        end
        current_length = current_length + 1;
        waitbar(current_length/experiment_length, progress_bar, ...
            sprintf('%d/%d',current_length, experiment_length));
    end
    
    if b_dual_decomp
        tic;
        output_dual = PN.optimizeNetSocialWelfare1(options);
        rt = toc;
        [tb, stbs] = saveStatTable(PN, output_dual, rt, type.Index, ...
            'dynamic-price');
        stat_price1(i, tb.Properties.VariableNames) = tb;
        for j = 1:num_type
            slice_stat_price1{j}(i, stbs.Properties.VariableNames) = stbs(j,:);
        end
        current_length = current_length + 1;
        waitbar(current_length/experiment_length, progress_bar, ...
            sprintf('%d/%d',current_length, experiment_length));
    end
    
    if b_price_adjust1
        [output_price, rt] = PN.optimizeResourcePrice([], options);
        [tb, stbs] = saveStatTable(PN, output_price, rt, type.Index, ...
            'dynamic-price');
        stat_price1(i, tb.Properties.VariableNames) = tb;
        for j = 1:num_type
            slice_stat_price1{j}(i, stbs.Properties.VariableNames) = stbs(j,:);
        end
        current_length = current_length + 1;
        waitbar(current_length/experiment_length, progress_bar, ...
            sprintf('%d/%d',current_length, experiment_length));
    end
    
    if b_price_adjust2
        [output_price, rt] = PN.optimizeResourcePriceNew([], options);
        [tb, stbs] = saveStatTable(PN, output_price, rt, type.Index, ...
            'dynamic-price');
        stat_price2(i, tb.Properties.VariableNames) = tb;
        for j = 1:num_type
            slice_stat_price2{j}(i, stbs.Properties.VariableNames) = stbs(j,:);
        end
        current_length = current_length + 1;
        waitbar(current_length/experiment_length, progress_bar, ...
            sprintf('%d/%d',current_length, experiment_length));
    end
    
    if b_resource_part
        tic;
        [output_part] = PN.resourcePartitionOptimization([], options);
        rt = toc;
        [tb, stbs] = saveStatTable(PN, output_price, rt, type.Index, ...
            'dynamic-price');
        stat_part(i, tb.Properties.VariableNames) = tb;
        for j = 1:num_type
            slice_stat_part{j}(i, stbs.Properties.VariableNames) = stbs(j,:);
        end
        current_length = current_length + 1;
        waitbar(current_length/experiment_length, progress_bar, ...
            sprintf('%d/%d',current_length, experiment_length));
    end
    
    if b_part_price
        tic;
        [output_partprice] = PN.partitionResourcePricing([], options);
        rt = toc;
        [tb, stbs] = saveStatTable(PN, output_price, rt, type.Index, ...
            'dynamic-price');
        stat_partprice(i, tb.Properties.VariableNames) = tb;
        for j = 1:num_type
            slice_stat_partprice{j}(i, stbs.Properties.VariableNames) = stbs(j,:);
        end
        current_length = current_length + 1;
        waitbar(current_length/experiment_length, progress_bar, ...
            sprintf('%d/%d',current_length, experiment_length));
    end
end
disp('Total arrival statistics:')
disp(RE.countArriveType);
clear e i T r p j;
close(progress_bar);