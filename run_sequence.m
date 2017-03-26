%% Process a sequence of network slice requests

%% Construct Network
% Initialize substrate network
PN1 = PhysicalNetwork(node_opt, link_opt, VNF_opt, net_opt);
PN1.slice_template = Slice.loadSliceTemplate(type.Index);
if b_static_slice
    PN2 = PN1.copy;
end
%% Network Slice as a Poisson Arrival Process
% Two types of network slice as a compound Poisson process.
num_events = arrival.Number*2;
num_type = length(type.Index);
num_fix_type = length(type.Fixed);
event_set = struct;
e = 1;
for i = (1+num_fix_type):num_type
    event_set(e).ServiceInterval = PN1.slice_template(i).ServiceInterval;
    event_set(e).Probability = PN1.slice_template(i).Probability;
    e = e + 1;
end
RE = RequestEvent(event_set, arrival, seed);
clear e i event_set;

seed2 = seed;
for t = 1:num_fix_type
    for s = 1:type.FixedCount(t)
        slice_opt = PN1.slice_template(type.Fixed(t));
        slice_opt.RandomSeed = seed;
        seed = seed + 1;
        PN1.AddSlice(slice_opt);
    end
end
clear s slice_opt;

%% Statistics Initialization
if b_single_optimal
    profit_approx.optimal = zeros(num_events,1);
    profit_accurate.optimal = zeros(num_events,1);
    utilization.optimal = zeros(num_events,1);
    runtime.optimal = zeros(num_events,2);
    rate_stat.optimal = cell(num_type,1);
    profit_stat.optimal = cell(num_type,1);
    for t = 1:num_type
        profit_stat.optimal{t} = table([],[],[],'VariableNames',{'Average','Max','Min'});
        rate_stat.optimal{t} = table([],[],[],'VariableNames',{'Average','Max','Min'});
    end
    profit_stat.optimal{num_type+1} = zeros(num_events,1);
end
if b_dual_decomp
    profit_approx.dual = zeros(num_events,1); %#ok<UNRCH>
    profit_accurate.dual = zeros(num_events,1);
    utilization.dual = zeros(num_events,1);
    runtime.dual = zeros(num_events,1);
    avg_rate.dual = zeros(num_events,1);
    rate_stat.dual = cell(num_type,1);
    profit_stat.dual = cell(num_type,1);
    for t = 1:num_type
        profit_stat.dual{t} = table([],[],[],'VariableNames',{'Average','Max','Min'});
        rate_stat.dual{t} = table([],[],[],'VariableNames',{'Average','Max','Min'});
    end
    profit_stat.optimal{num_type+1} = zeros(num_events,1);
end
if b_price_adjust
    profit_approx.price = zeros(num_events,1);
    profit_accurate.price = zeros(num_events,1);
    utilization.price = zeros(num_events,1);
    runtime.price = zeros(num_events,2);
    avg_rate.price = zeros(num_events,1);
    rate_stat.price = cell(num_type,1);
    profit_stat.price = cell(num_type,1);
    for t = 1:num_type
        profit_stat.price{t} = table([],[],[],'VariableNames',{'Average','Max','Min'});
        rate_stat.price{t} = table([],[],[],'VariableNames',{'Average','Max','Min'});
    end
    profit_stat.price{num_type+1} = zeros(num_events,1);
end
if b_price_adjust2
    profit_approx.price2 = zeros(num_events,1);
    profit_accurate.price2 = zeros(num_events,1);
    utilization.price2 = zeros(num_events,1);
    runtime.price2 = zeros(num_events,1);
    avg_rate.price2 = zeros(num_events,1);
    rate_stat.price2 = cell(num_type,1);
    profit_stat.price2 = cell(num_type,1);
    for t = 1:num_type
        profit_stat.price2{t} = table([],[],[],'VariableNames',{'Average','Max','Min'});
        rate_stat.price2{t} = table([],[],[],'VariableNames',{'Average','Max','Min'});
    end
    profit_stat.price2{num_type+1} = zeros(num_events,1);
end
if b_resource_part
    profit_approx.part = zeros(num_events,1);
    profit_accurate.part = zeros(num_events,1);
    utilization.part = zeros(num_events,1);
    runtime.part = zeros(num_events,1);
    rate_stat.part = cell(num_type,1);
    profit_stat.part = cell(num_type,1);
    for t = 1:num_type
        profit_stat.part{t} = table([],[],[],'VariableNames',{'Average','Max','Min'});
        rate_stat.part{t} = table([],[],[],'VariableNames',{'Average','Max','Min'});
    end
    profit_stat.part{num_type+1} = zeros(num_events,1);
end
if b_part_price
    profit_approx.partprice = zeros(num_events,1);
    profit_accurate.partprice = zeros(num_events,1);
    utilization.partprice = zeros(num_events,1);
    runtime.partprice = zeros(num_events,1);
    rate_stat.partprice = cell(num_type,1);
    profit_stat.partprice = cell(num_type,1);
    for t = 1:num_type
        profit_stat.partprice{t} = table([],[],[],'VariableNames',{'Average','Max','Min'});
        rate_stat.partprice{t} = table([],[],[],'VariableNames',{'Average','Max','Min'});
    end
    profit_stat.partprice{num_type+1} = zeros(num_events,1);
end
if b_static_slice
    profit_approx.static = zeros(num_events,1);
    profit_accurate.static = zeros(num_events,1);
    utilization.static = zeros(num_events,1);
    runtime.static = zeros(num_events,1);
    rate_stat.static = cell(num_type,1);
    profit_stat.static = cell(num_type+1,1);
    for t = 1:num_type
        profit_stat.static{t} = table([],[],[],'VariableNames',{'Average','Max','Min'});
        rate_stat.static{t} = table([],[],[],'VariableNames',{'Average','Max','Min'});
    end
    profit_stat.static{num_type+1} = zeros(num_events,1);
    stat.number_flows_static = zeros(num_events,1);
    stat.number_reject = zeros(num_events,3);
    stat.number_part_reject = zeros(num_events,3);
    stat.number_slices_static = zeros(num_events,3);

    link_uc = PN2.getLinkField('UnitCost');
    node_uc = PN2.getNodeField('UnitCost');
    link_price = (link_uc + PN2.phis_l) * (1 + options.PricingFactor);
    node_price = (node_uc + PN2.phis_n) * (1 + options.PricingFactor);
    PN2.setLinkField('Price', link_price);
    PN2.setNodeField('Price', node_price);
    % TODO: the persistant slices can be optimized.
    static_opts = options;
    static_opts.Method = 'slice-price';
    for t = 1:num_fix_type
        for s = 1:type.FixedCount(t)
            slice_opt = PN2.slice_template(type.Fixed(t));
            slice_opt.RandomSeed = seed2 + s;
            seed2 = seed2 + 1;
            slice_opt.method = 'static-slicing';
            slice_opt.admit_ploicy = 'reject-flow';
            sl = PN2.AddSlice(slice_opt);
            if ~isempty(sl)
                PN2.staticSlicing(sl, static_opts);
                stat.number_slices_static(1,1) = stat.number_slices_static(1,1) + 1;
            end
        end
    end
    clear link_uc node_uc link_price node_price slice_opt s seed2;
end
clear t s;
stat.number_slices = zeros(num_events,3);
for t = 1:num_fix_type
    stat.number_slices(:,type.Fixed(t)) = type.FixedCount(t);
end
stat.number_flows = zeros(num_events,1);
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
    stat.number_slices(i,T.Value+num_fix_type) = ...
        stat.number_slices(i,T.Value+num_fix_type) + T.Count';
    stat.times(i) = e.Time;
    if strcmp(e.Description, 'arrival')
        %%%
        % |e.Type+num_fix_type | is the index of slice type in the array.
        % |e.Type| is the index of event type.
        fprintf('\tSlice type %d.\n', e.Type+num_fix_type);
        slice_opt = PN1.slice_template(e.Type+num_fix_type);
        slice_opt.RandomSeed = seed;
        seed = seed + 1;
        slice_opt.Identifier = e.Identifier;
        PN1.AddSlice(slice_opt);
        stat.number_flows(i) = PN1.NumberFlows;
        if b_static_slice
            slice_opt = PN2.slice_template(e.Type+num_fix_type);
            slice_opt.RandomSeed = seed-1;
            slice_opt.Identifier = e.Identifier;
            slice_opt.method = 'static-slicing';
            slice_opt.admit_ploicy = 'reject-flow';
            sl = PN2.AddSlice(slice_opt);
            stat.number_flows_static(i) = PN2.NumberFlows;
            if i > 1
                stat.number_slices_static(i,:) = stat.number_slices_static(i-1,:);
            end
            if isempty(sl)
                stat.number_reject(i, e.Type+num_fix_type) = 1;
                stat.number_part_reject(i, e.Type+num_fix_type) = 1;
            else
                if sl.NumberFlows < PN2.slice_template(e.Type+num_fix_type).NumberFlows
                    stat.number_part_reject(i, e.Type+num_fix_type) = 1;
                end
                stat.number_slices_static(i,e.Type+num_fix_type) = ...
                    stat.number_slices_static(i,e.Type+num_fix_type) + 1;
            end
            
            tic;
            [output_static] = PN2.staticSlicing(sl, static_opts);
            runtime.static(i) = toc;
            profit_approx.static(i) = output_static.welfare_approx;
            profit_accurate.static(i) = output_static.welfare_accurate;
            utilization.static(i) = PN2.utilizationRatio;
            for j = 1:num_type
                [p,r] = PN2.statSlice(type.Index(j), output_static.profit.ApproximatePrice);
                profit_stat.static{j}(end+1,:) = p;
                rate_stat.static{j}(end+1,:) = r;
            end
            profit_stat.static{end}(i,:) = output_static.profit.ApproximatePrice(end);
        end
        if b_single_optimal
            [output_optimal, rt] = PN1.singleSliceOptimization(options);
            runtime.optimal(i,2) = rt.Serial;
            runtime.optimal(i,1) = rt.Parallel;
            profit_approx.optimal(i) = output_optimal.welfare_approx;
            profit_accurate.optimal(i) = output_optimal.welfare_accurate;
            utilization.optimal(i) = PN1.utilizationRatio;
            for j = 1:num_type
                [p,r] = PN1.statSlice(type.Index(j), output_optimal.profit.ApproximatePrice);
                profit_stat.optimal{j}(end+1,:) = p;
                rate_stat.optimal{j}(end+1,:) = r;
            end
            profit_stat.optimal{end}(i,:) = output_optimal.profit.ApproximatePrice(end);
        end
        if b_dual_decomp
            tic;
            output_dual = PN1.optimizeNetSocialWelfare1(options);
            runtime.dual(i) = toc;
            profit_approx.dual(i) = output_dual.welfare_approx;
            profit_accurate.dual(i) = output_dual.welfare_accurate;
            utilization.dual(i) = PN1.utilizationRatio;
            for j = 1:num_type
                [p,r] = PN1.statSlice(type.Index(j), output_dual.profit.ApproximatePrice);
                profit_stat.dual{j}(end+1,:) = p;
                rate_stat.dual{j}(end+1,:) = r;
            end
            profit_stat.dual{end}(i,:) = output_dual.profit.ApproximatePrice(end);
        end
        if b_price_adjust
            [output_price, rt] = PN1.optimizeResourcePrice([], options);
            runtime.price(i,2) = rt.Serial;
            runtime.price(i,1) = rt.Parallel;
            profit_approx.price(i) = output_price.welfare_approx;
            profit_accurate.price(i) = output_price.welfare_accurate;
            utilization.price(i) = PN1.utilizationRatio;
            for j = 1:num_type
                [p,r] = PN1.statSlice(type.Index(j), output_price.profit.ApproximatePrice);
                profit_stat.price{j}(end+1,:) = p;
                rate_stat.price{j}(end+1,:) = r;
            end
            profit_stat.price{end}(i,:) = output_price.profit.ApproximatePrice(end);
        end
        if b_price_adjust2
            tic;
            [output_price] = PN1.optimizeResourcePrice2([], options);
            runtime.price2(i) = toc;
            profit_approx.price2(i) = output_price.welfare_approx;
            profit_accurate.price2(i) = output_price.welfare_accurate;
            utilization.price2(i) = PN1.utilizationRatio;
            for j = 1:num_type
                [p,r] = PN1.statSlice(type.Index(j), output_price.profit.ApproximatePrice);
                profit_stat.price2{j}(end+1,:) = p;
                rate_stat.price2{j}(end+1,:) = r;
            end
            profit_stat.price2{end}(i,:) = output_price.profit.ApproximatePrice(end);
        end
        if b_resource_part
            tic;
            [output_part] = PN1.resourcePartitionOptimization([], options);
            runtime.part(i) = toc;
            profit_approx.part(i) = output_part.welfare_approx;
            profit_accurate.part(i) = output_part.welfare_accurate;
            utilization.part(i) = PN1.utilizationRatio;
            for j = 1:num_type
                [p,r] = PN1.statSlice(type.Index(j), output_part.profit.ApproximatePrice);
                profit_stat.part{j}(end+1,:) = p;
                rate_stat.part{j}(end+1,:) = r;
            end
            profit_stat.part{end}(i,:) = output_part.profit.ApproximatePrice(end);
        end
        if b_part_price
            tic;
            [output_partprice] = PN1.partitionResourcePricing([], options);
            runtime.partprice(i) = toc;
            profit_approx.partprice(i) = output_partprice.welfare_approx;
            profit_accurate.partprice(i) = output_partprice.welfare_accurate;
            utilization.partprice(i) = PN1.utilizationRatio;
            for j = 1:num_type
                [p,r] = PN1.statSlice(type.Index(j), output_partprice.profit.ApproximatePrice);
                profit_stat.partprice{j}(end+1,:) = p;
                rate_stat.partprice{j}(end+1,:) = r;
            end
            profit_stat.partprice{end}(i,:) = output_partprice.profit.ApproximatePrice(end);
        end
    else
        fprintf('\tDeparting slice type %d.\n', e.Type+num_fix_type);
        PN1.RemoveSlice(e.Id);
        stat.number_flows(i) = PN1.NumberFlows;
        if PN1.NumberSlices > 0
            if b_single_optimal
                [output_optimal, rt] = PN1.singleSliceOptimization(options);
                runtime.optimal(i,2) = rt.Serial;
                runtime.optimal(i,1) = rt.Parallel;
                profit_approx.optimal(i) = output_optimal.welfare_approx;
                profit_accurate.optimal(i) = output_optimal.welfare_accurate;
                utilization.optimal(i) = PN1.utilizationRatio;
                for j = 1:num_type
                    [p,r] = PN1.statSlice(type.Index(j), output_optimal.profit.ApproximatePrice);
                    profit_stat.optimal{j}(end+1,:) = p;
                    rate_stat.optimal{j}(end+1,:) = r;
                end
                profit_stat.optimal{end}(i,:) = output_optimal.profit.ApproximatePrice(end);
            end
            if b_dual_decomp
                tic;
                [output_dual] = PN1.optimizeNetSocialWelfare1(options);
                runtime.dual(i) = toc;
                profit_approx.dual(i) = output_dual.welfare_approx;
                profit_accurate.dual(i) = output_dual.welfare_accurate;
                utilization.dual(i) = PN1.utilizationRatio;
                for j = 1:num_type
                    [p,r] = PN1.statSlice(type.Index(j), output_dual.profit.ApproximatePrice);
                    profit_stat.dual{j}(end+1,:) = p;
                    rate_stat.dual{j}(end+1,:) = r;
                end
                profit_stat.dual{end}(i,:) = output_dual.profit.ApproximatePrice(end);
            end
            if b_price_adjust
                [output_price, rt] = PN1.optimizeResourcePrice([], options);
                runtime.price(i,2) = rt.Serial;
                runtime.price(i,1) = rt.Parallel;
                profit_approx.price(i) = output_price.welfare_approx;
                profit_accurate.price(i) = output_price.welfare_accurate;
                utilization.price(i) = PN1.utilizationRatio;
                for j = 1:num_type
                    [p,r] = PN1.statSlice(type.Index(j), output_price.profit.ApproximatePrice);
                    profit_stat.price{j}(end+1,:) = p;
                    rate_stat.price{j}(end+1,:) = r;
                end
                profit_stat.price{end}(i,:) = output_price.profit.ApproximatePrice(end);
            end
            if b_price_adjust2
                tic;
                [output_price] = PN1.optimizeResourcePrice3([], options);
                runtime.price2(i) = toc;
                profit_approx.price2(i) = output_price.welfare_approx;
                profit_accurate.price2(i) = output_price.welfare_accurate;
                utilization.price2(i) = PN1.utilizationRatio;
                for j = 1:num_type
                    [p,r] = PN1.statSlice(type.Index(j), output_price.profit.ApproximatePrice);
                    profit_stat.price2{j}(end+1,:) = p;
                    rate_stat.price2{j}(end+1,:) = r;
                end
                profit_stat.price2{end}(i,:) = output_price.profit.ApproximatePrice(end);
            end
            if b_resource_part
                tic;
                [output_part] = PN1.resourcePartitionOptimization([], options);
                runtime.part(i) = toc;
                profit_approx.part(i) = output_part.welfare_approx;
                profit_accurate.part(i) = output_part.welfare_accurate;
                utilization.part(i) = PN1.utilizationRatio;
                for j = 1:num_type
                    [p,r] = PN1.statSlice(type.Index(j), output_part.profit.ApproximatePrice);
                    profit_stat.part{j}(end+1,:) = p;
                    rate_stat.part{j}(end+1,:) = r;
                end
                profit_stat.part{end}(i,:) = output_part.profit.ApproximatePrice(end);
            end
            if b_part_price
                tic;
                [output_partprice] = PN1.partitionResourcePricing([], options);
                runtime.partprice(i) = toc;
                profit_approx.partprice(i) = output_partprice.welfare_approx;
                profit_accurate.partprice(i) = output_partprice.welfare_accurate;
                utilization.partprice(i) = PN1.utilizationRatio;
                for j = 1:num_type
                    [p,r] = PN1.statSlice(type.Index(j), output_partprice.profit.ApproximatePrice);
                    profit_stat.partprice{j}(end+1,:) = p;
                    rate_stat.partprice{j}(end+1,:) = r;
                end
                profit_stat.partprice{end}(i,:) = output_partprice.profit.ApproximatePrice(end);
            end

        end
        sl = PN2.RemoveSlice(e.Id);       
        if PN2.NumberSlices > 0 && b_static_slice
            stat.number_slices_static(i,:) = stat.number_slices_static(i-1,:);
            if isempty(sl)
                % Nothing should be changed.
                fprintf('Info: slice %d have been rejected earlier.\n', e.Id);
                runtime.static(i) = 0;
            else
                % When a slice is removed, staticSlicing update the network state.
                tic;
                [output_static] = PN2.staticSlicing([], static_opts);
                runtime.static(i) = toc;
                %%%
                % It is determined that i > 1.
                stat.number_slices_static(i,e.Type+num_fix_type) = ...
                    stat.number_slices_static(i-1,e.Type+num_fix_type) - 1;
            end
            stat.number_flows_static(i) = PN2.NumberFlows;
            profit_approx.static(i) = output_static.welfare_approx;
            profit_accurate.static(i) = output_static.welfare_accurate;
            utilization.static(i) = PN2.utilizationRatio;
            for j = 1:num_type
                [p,r] = PN2.statSlice(type.Index(j), output_static.profit.ApproximatePrice);
                profit_stat.static{j}(end+1,:) = p;
                rate_stat.static{j}(end+1,:) = r;
            end
            profit_stat.static{end}(i,:) = output_static.profit.ApproximatePrice(end);
        end
    end
end
disp('Total arrival statistics:')
disp(RE.countArriveType);
clear e i T r p j;