%%
% Giving a combination of different types of slices, calculate the performace metrics.
% See also <run_test1.m>

%% NOTE
% The data structres that save the output results are deprecated. Now we use table to save
% the output results (see <run_exp41.m> and <run_exp42.m>)

%% Configurable parameters in this experiment
link_opt.CostUnit = 100;        % 75 | 100 | 150
node_opt.CostUnit = 300;        % 250 | 300 | 500
options.PricingFactor = 2;      % 3 | 2 | 1 ,

%% Specification of VNFs and Network Slices
link_opt.delay = LinkDelayOption.Random;
link_opt.cost = LinkCostOption.CapacityInverse;
link_opt.CapacityFactor = 30;
link_opt.RandomSeed = 20170423;
link_opt.delay = LinkDelayOption.Random;
net_opt.delta = 0.7;

node_opt.model = NetworkModel.Sample1;
node_opt.capacity = NodeCapacityOption.BandwidthProportion;
node_opt.cost = NodeCostOption.CapacityInverse;
node_opt.capacity_factor = 1.5;     % [0.3; 0.5; 0.8; 1]
net_opt.AdmitPolicy = 'reject-flow';
net_opt.NetworkClass = 'PhysicalNetwork';

VNF_opt.Number = 4;            % number of VNF type
VNF_opt.Model = VNFIntegrateModel.AllInOne;
VNF_opt.RandomSeed = [20161101 0];

%% Slice Combination
slice_config.Type = [1; 2; 3];
% slice_config.RatioSet = {[1 1 1],[4 1 1],[1 4 1], [1 1 4]};
slice_config.RatioSet = {[1 1 1]};
slice_config.Count = 4:4:36;

%% Algorithm options
options.ProfitType = {'ApproximatePrice','AccuratePrice'};
options.WelfareType = {'Accurate', 'Approximate'};

%% Experiment Control
options.Display = 'off';
b_single_optimal = true;
b_price_adjust1 = true;
b_price_adjust2 = true;
b_static_slice = true;
experiment_length = length(slice_config.RatioSet) * length(slice_config.Count)*...
    (b_single_optimal+b_price_adjust1+b_price_adjust2+b_static_slice);
current_length = 0;
progress_bar = waitbar(current_length/experiment_length, ...
    sprintf('%d/%d',current_length, experiment_length));
progress_bar.Name = 'Computing Progress';
pause(0.01);

num_point = length(slice_config.Count);
num_type = length(slice_config.Type);
num_config = length(slice_config.RatioSet);
output_results = cell(num_config,1);
for exp_id = 1:num_config
    slice_config.Ratio = slice_config.RatioSet{exp_id};
    seed = 20170430;
    
    %% Statistics Initialization
    if ~isfield(net_opt, 'AdmitPolicy')
        net_opt.AdmitPolicy = 'reject-flow';
    end
    Average = zeros(num_point, 1);
    Max = zeros(num_point, 1);
    Min = zeros(num_point, 1);
    StdDev = zeros(num_point, 1);
    if b_single_optimal
        profit_approx.optimal = zeros(num_point,2); % upper bound and with fixed price
        profit_accurate.optimal = zeros(num_point,2);
        cost.optimal = zeros(num_point,1);
        utilization.optimal = zeros(num_point,1);
        node_utilization.optimal = table(Average, Max, Min, StdDev);
        link_utilization.optimal = table(Average, Max, Min, StdDev);
        runtime.optimal = zeros(num_point,2);
        rate_stat.optimal = cell(num_type,1);
        profit_stat.optimal = cell(num_type+1,1);
        for t = 1:num_type
            profit_stat.optimal{t} = table(Average, Max, Min, StdDev);
            rate_stat.optimal{t} = table(Average, Max, Min, StdDev);
        end
        profit_stat.optimal{num_type+1} = zeros(num_point,1);
    end
    if b_price_adjust1
        profit_approx.price1 = zeros(num_point,1);
        profit_accurate.price1 = zeros(num_point,1);
        cost.price1 = zeros(num_point,1);
        utilization.price1 = zeros(num_point,1);
        node_utilization.price1 = table(Average, Max, Min, StdDev);
        link_utilization.price1 = table(Average, Max, Min, StdDev);
        runtime.price1 = zeros(num_point,2);
        avg_rate.price1 = zeros(num_point,1);
        rate_stat.price1 = cell(num_type,1);
        profit_stat.price1 = cell(num_type+1,1);
        for t = 1:num_type
            profit_stat.price1{t} = table(Average, Max, Min, StdDev);
            rate_stat.price1{t} = table(Average, Max, Min, StdDev);
        end
        profit_stat.price1{num_type+1} = zeros(num_point,1);
    end
    if b_price_adjust2
        profit_approx.price2 = zeros(num_point,1);
        profit_accurate.price2 = zeros(num_point,1);
        cost.price2 = zeros(num_point,1);
        utilization.price2 = zeros(num_point,1);
        node_utilization.price2 = table(Average, Max, Min, StdDev);
        link_utilization.price2 = table(Average, Max, Min, StdDev);
        runtime.price2 = zeros(num_point,2);
        avg_rate.price2 = zeros(num_point,1);
        rate_stat.price2 = cell(num_type,1);
        profit_stat.price2 = cell(num_type+1,1);
        for t = 1:num_type
            profit_stat.price2{t} = table(Average, Max, Min, StdDev);
            rate_stat.price2{t} = table(Average, Max, Min, StdDev);
        end
        profit_stat.price2{num_type+1} = zeros(num_point,1);
    end
    if b_static_slice
        profit_approx.static = zeros(num_point,1);
        profit_accurate.static = zeros(num_point,1);
        cost.static = zeros(num_point,1);
        utilization.static = zeros(num_point,1);
        node_utilization.static = table(Average, Max, Min, StdDev);
        link_utilization.static = table(Average, Max, Min, StdDev);
        runtime.static = zeros(num_point,1);
        rate_stat.static = cell(num_type,1);
        profit_stat.static = cell(num_type+1,1);
        for t = 1:num_type
            profit_stat.static{t} = table(Average, Max, Min, StdDev);
            rate_stat.static{t} = table(Average, Max, Min, StdDev);
        end
        profit_stat.static{num_type+1} = zeros(num_point,1);
        stat.number_flows_static = zeros(num_point,1);
        stat.number_reject = zeros(num_point,3);
        stat.number_part_reject = zeros(num_point,3);
        stat.number_slices_static = zeros(num_point,3);
        static_opts = options;
        static_opts.Method = 'slice-price';
    end
    %%
    slice_sequence = zeros(slice_config.Count(end),1);
    slice_type = zeros(slice_config.Count(end),1);
    total_num_slices_0 = 0;
    stat.number_slices = zeros(num_point,3);
    for point_id = 1:num_point
        seed_static = seed;
        seed_dynamic = seed;
        fprintf('(%s)Iteration %d, experiment number %d.\n', datestr(now), exp_id, point_id);
        %% Specify the combination of slices
        % The combination of slices of the current point is based on that of the previous
        % point. So we should guarantee that the common part of the two point should be
        % the same, i.e. the common slices should have the same sequence and flow demand.
        total_num_slices = slice_config.Count(point_id);
        for j = 1:(num_type-1)
            ct = round(total_num_slices * slice_config.Ratio(j)/sum(slice_config.Ratio));
            if ct > total_num_slices - sum(stat.number_slices(point_id,:))
                ct = ct - 1;
            end
            stat.number_slices(point_id,j) = ct;
        end
        stat.number_slices(point_id,num_type) = total_num_slices - sum(stat.number_slices(point_id,:));
        if point_id == 1
            delta_number_slices = stat.number_slices(point_id,:);
        else
            delta_number_slices = stat.number_slices(point_id,:) - stat.number_slices(point_id-1,:);
            delta_number_slices(delta_number_slices<0) = 0;
        end
        rng(seed + point_id);
        delta_total = sum(delta_number_slices);
        delta_slice_sequence = unique_randi(delta_total, delta_total, 'stable');
        slice_sequence((total_num_slices_0+1):total_num_slices) = ...
            delta_slice_sequence + total_num_slices_0;
        delta_slice_type = zeros(delta_total,1);
        t1 = 0;
        for j = 1:num_type
            t2 = t1 + delta_number_slices(j);
            delta_slice_type((t1+1):t2) = j;
            t1 = t2;
        end
        slice_type((total_num_slices_0+1):total_num_slices) = delta_slice_type;
        total_num_slices_0 = total_num_slices;
        %     slice_type = zeros(total_num_slices,1);
        %     slice_count = stat.number_slices(i,:);
        %     j = 1;
        %     while nnz(slice_count) > 0
        %         for k = 1:num_type
        %             if slice_count(k)>0
        %                 slice_type(j) = k;
        %                 j = j + 1;
        %                 slice_count(k) = slice_count(k) - 1;
        %             end
        %         end
        %     end
        
        if b_static_slice
            if isfield(net_opt, 'NetworkClass') && ...
                    strcmp(net_opt.NetworkClass, 'CloudNetwork')
                PN = CloudNetwork(node_opt, link_opt, VNF_opt, net_opt);
            else
                PN = PhysicalNetwork(node_opt, link_opt, VNF_opt, net_opt);
            end
            PN.slice_template = Slice.loadSliceTemplate(slice_config.Type);
            link_price = PN.getLinkCost * (1 + options.PricingFactor);
            node_price = PN.getNodeCost * (1 + options.PricingFactor);
            PN.setLinkField('Price', link_price);
            PN.setDataCenterField('Price', node_price);
            rt = zeros(total_num_slices,1);
            for j = 1:total_num_slices
                s_type = slice_type(slice_sequence(j));
                slice_opt = PN.slice_template(s_type);
                slice_opt.RandomSeed = seed_static;
                seed_static = seed_static + 1;
                slice_opt.method = 'static-slicing';
                slice_opt.admit_ploicy = net_opt.AdmitPolicy;
                sl = PN.AddSlice(slice_opt);
                stat.number_slices_static(point_id,s_type) = ...
                    stat.number_slices_static(point_id,s_type) + 1;
                if isempty(sl)
                    stat.number_reject(point_id, s_type) = stat.number_reject(point_id, s_type) + 1;
                    stat.number_part_reject(point_id, s_type) = stat.number_part_reject(point_id, s_type) + 1;
                elseif sl.NumberFlows < PN.slice_template(s_type).NumberFlows
                    stat.number_part_reject(point_id, s_type) = stat.number_part_reject(point_id, s_type) + 1;
                end
                fprintf('(%s)Static slicing: adding a slice %d.\n', datestr(now), ...
                    slice_opt.Type);
                tic;
                [output_static] = PN.staticSlicing(sl, static_opts);
                rt(j) = toc;
            end
            runtime.static(point_id) = max(rt(rt~=0));
            stat.number_flows_static(point_id) = PN.NumberFlows;
            profit_approx.static(point_id) = output_static.welfare_approx;
            profit_accurate.static(point_id) = output_static.welfare_accurate;
            utilization.static(point_id) = PN.utilizationRatio;
            [r1, r2, r3, r4] = PN.nodeUtilization;
            node_utilization.static{point_id,:} = [r1, r2, r3, r4];
            [r1, r2, r3, r4] = PN.linkUtilization;
            link_utilization.static{point_id,:} = [r1, r2, r3, r4];
            cost.static(point_id) = PN.getNetworkCost([],[]);
            for j = 1:num_type
                [p,r] = PN.statSlice(slice_config.Type(j), output_static.profit.ApproximatePrice);
                profit_stat.static{j}{point_id,:} = p;
                rate_stat.static{j}{point_id,:} = r;
            end
            profit_stat.static{end}(point_id) = output_static.profit.ApproximatePrice(end);
            current_length = current_length + 1;
            waitbar(current_length/experiment_length, progress_bar, ....
                sprintf('%d/%d',current_length, experiment_length));
        end
        
        if b_single_optimal || b_price_adjust1 || b_price_adjust2
            if isfield(net_opt, 'NetworkClass') && strcmp(net_opt.NetworkClass, 'CloudNetwork')
                PN = CloudNetwork(node_opt, link_opt, VNF_opt, net_opt);
            else
                PN = PhysicalNetwork(node_opt, link_opt, VNF_opt, net_opt);
            end
            PN.slice_template = Slice.loadSliceTemplate(slice_config.Type);
            fprintf('(%s)Non-static slicing: adding all slices.\n', datestr(now));
            for j = 1:total_num_slices
                s_type = slice_type(slice_sequence(j));
                slice_opt = PN.slice_template(s_type);
                slice_opt.RandomSeed = seed_dynamic;
                seed_dynamic = seed_dynamic + 1;
                PN.AddSlice(slice_opt);
            end
            stat.number_flows(point_id) = PN.NumberFlows;
        end
        %%
        if b_single_optimal
            fprintf('(%s)Global SPP.\n', datestr(now));
            options.Method = 'normal';
            [output_optimal, rt] = PN.singleSliceOptimization(options);
            runtime.optimal(point_id,2) = rt.Serial;
            runtime.optimal(point_id,1) = rt.Parallel;
            profit_approx.optimal(point_id,:) = ...
                [output_optimal.welfare_approx_optimal, output_optimal.welfare_approx];
            profit_accurate.optimal(point_id,:) = ...
                [output_optimal.welfare_accurate_optimal, output_optimal.welfare_accurate];
            utilization.optimal(point_id) = PN.utilizationRatio;
            [r1, r2, r3, r4] = PN.nodeUtilization;
            node_utilization.optimal{point_id,:} = [r1, r2, r3, r4];
            [r1, r2, r3, r4] = PN.linkUtilization;
            link_utilization.optimal{point_id,:} = [r1, r2, r3, r4];
            cost.optimal(point_id) = PN.getNetworkCost([],[]);
            for j = 1:num_type
                [p,r] = PN.statSlice(slice_config.Type(j), output_optimal.profit.ApproximatePrice);
                profit_stat.optimal{j}{point_id,:} = p;
                rate_stat.optimal{j}{point_id,:} = r;
            end
            profit_stat.optimal{end}(point_id) = output_optimal.profit.ApproximatePrice(end);
            current_length = current_length + 1;
            waitbar(current_length/experiment_length, progress_bar, ...
                sprintf('%d/%d',current_length, experiment_length));
        end
        %%
        if b_price_adjust1
            fprintf('(%s)Pricing-1.\n', datestr(now));
            [output_price, rt] = PN.optimizeResourcePrice([], options);
            runtime.price1(point_id,2) = rt.Serial;
            runtime.price1(point_id,1) = rt.Parallel;
            profit_approx.price1(point_id) = output_price.welfare_approx;
            profit_accurate.price1(point_id) = output_price.welfare_accurate;
            utilization.price1(point_id) = PN.utilizationRatio;
            [r1, r2, r3, r4] = PN.nodeUtilization;
            node_utilization.price1{point_id,:} = [r1, r2, r3, r4];
            [r1, r2, r3, r4] = PN.linkUtilization;
            link_utilization.price1{point_id,:} = [r1, r2, r3, r4];
            cost.price(point_id) = PN.getNetworkCost([],[]);
            for j = 1:num_type
                [p,r] = PN.statSlice(slice_config.Type(j), output_price.profit.ApproximatePrice);
                profit_stat.price1{j}{point_id,:} = p;
                rate_stat.price1{j}{point_id,:} = r;
            end
            profit_stat.price1{end}(point_id) = output_price.profit.ApproximatePrice(end);
            current_length = current_length + 1;
            waitbar(current_length/experiment_length, progress_bar, ...
                sprintf('%d/%d',current_length, experiment_length));
        end
        %%
        if b_price_adjust2
            fprintf('(%s)Pricing-2.\n', datestr(now));
            [output_price, rt] = PN.optimizeResourcePriceNew([], options);
            runtime.price2(point_id,2) = rt.Serial;
            runtime.price2(point_id,1) = rt.Parallel;
            profit_approx.price2(point_id) = output_price.welfare_approx;
            profit_accurate.price2(point_id) = output_price.welfare_accurate;
            utilization.price2(point_id) = PN.utilizationRatio;
            [r1, r2, r3, r4] = PN.nodeUtilization;
            node_utilization.price2{point_id,:} = [r1, r2, r3, r4];
            [r1, r2, r3, r4] = PN.linkUtilization;
            link_utilization.price2{point_id,:} = [r1, r2, r3, r4];
            cost.price2(point_id) = PN.getNetworkCost([],[]);
            for j = 1:num_type
                [p,r] = PN.statSlice(slice_config.Type(j), output_price.profit.ApproximatePrice);
                profit_stat.price2{j}{point_id,:} = p;
                rate_stat.price2{j}{point_id,:} = r;
            end
            profit_stat.price2{end}(point_id) = output_price.profit.ApproximatePrice(end);
            current_length = current_length + 1;
            waitbar(current_length/experiment_length, progress_bar, ...
                sprintf('%d/%d',current_length, experiment_length));
        end
    end
    
    %%
    % if ~exist('output_results', 'var')
    %     output_results = cell;
    % end
    if isequal(slice_config.Ratio, [1 1 1])
        output_results{1} = {profit_accurate, profit_approx, profit_stat, rate_stat, ...
            runtime, stat, utilization, slice_config, ...
            node_utilization, link_utilization, cost};
    end
    if isequal(slice_config.Ratio, [4 1 1])
        output_results{2} = {profit_accurate, profit_approx, profit_stat, rate_stat, ...
            runtime, stat, utilization, slice_config, ...
            node_utilization, link_utilization, cost};
    end
    if isequal(slice_config.Ratio, [1 4 1])
        output_results{3} = {profit_accurate, profit_approx, profit_stat, rate_stat, ...
            runtime, stat, utilization, slice_config, ...
            node_utilization, link_utilization, cost};
    end
    if isequal(slice_config.Ratio, [1 1 4])
        output_results{4} = {profit_accurate, profit_approx, profit_stat, rate_stat, ...
            runtime, stat, utilization, slice_config, ...
            node_utilization, link_utilization, cost};
    end
end
close(progress_bar);
%% Figure
cost_limit = [1000 6000];
profit_limit = [0 9000];

%%
% description = 'Method2, max SP profit, low cost, low price';
% save('output_4301.mat', 'output_results', 'description');

% description = 'Method2, average profit ratio, low cost, low price';
% save('output_4302.mat', 'output_results', 'description');

% description = 'Method2, average profit ratio, middle cost, middle price';
% save('output_4303.mat', 'output_results', 'description');

% description = 'Method2, average profit ratio, high cost, high price';
% save('output_4304.mat', 'output_results', 'description');

% description = 'Method2, average profit ratio, high cost, middle price';
% save('output_4305.mat', 'output_results', 'description');

% description = 'Method2, average profit ratio, high cost, low price';
% save('output_4306.mat', 'output_results', 'description');

% description = 'Method2, average profit ratio, middle cost, high price';
% save('output_4307.mat', 'output_results', 'description');

% description = 'Method2, average profit ratio, middle cost, low price';
% save('output_4308.mat', 'output_results', 'description');
 
% description = 'Method2, average profit ratio, low cost, high price';
% save('output_4309.mat', 'output_results', 'description');

% description = 'Method2, average profit ratio, low cost, middle price';
% save('output_4310.mat', 'output_results', 'description');

% description = 'Method2, max SP profit, low cost, middle price';
% save('output_4311.mat', 'output_results', 'description');

% description = 'Method2, max SP profit, low cost, high price';
% save('Results\output_4312.mat', 'output_results', 'description', 'link_opt', 'node_opt',...
%     'net_opt', 'VNF_opt', 'slice_config', 'options');