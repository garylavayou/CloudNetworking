%%
% Giving a combination of different types of slices, calculate the
% performace metrics. 
%
% Network Model: Sample1/Sample2/SD_RAN

%% NOTE
% Compared with <run_exp301.m>, in this program, each iteration will not
% repeatly add slices to the network, only new slices will be added to the
% network based on the last one iteration. Thus, the computation time of
% static slicing and path calculation can be saved.

%% Specification of Substrate Network
clearvars -except progress_bar;
clear global;
global DEBUG;
DEBUG = false;
node_opt.Model = NetworkModel.Sample1; % NetworkModel.Sample2

%% Configurable parameters in this experiment
switch node_opt.Model
	case NetworkModel.Sample1
		link_opt.CostUnit = 100;        % 75 | 100 | 150
		node_opt.CostUnit = 300;        % 250 | 300 | 500
	case NetworkModel.SD_RAN
		link_opt.CostUnit = 0.2;
		node_opt.CostUnit = 0.1;
end
net_opt.PricingFactor = 2;      % 1 | 2 | 3 ,
net_opt.Threshold = 'average';
% Experiment Control
b_single_optimal = true;
b_price_adjust1 = true;
b_price_adjust2 = true;
b_static_slice = true;
b_save = false;
b_plot = false;

%% Specification of VNFs and Network Slices
link_opt.RandomSeed = 20170423;
link_opt.delay = LinkDelayOption.Random;
switch node_opt.Model
	case NetworkModel.Sample1
		link_opt.cost = LinkCostOption.CapacityInverse;
		link_opt.CapacityFactor = 30;
		node_opt.capacity = NodeCapacityOption.BandwidthProportion;
		node_opt.cost = NodeCostOption.CapacityInverse;
		node_opt.CapacityFactor = 1.5;     % [0.3; 0.5; 0.8; 1]
		net_opt.ClassName = 'CloudNetwork';
		VNF_opt.Number = 4;            % number of VNF type
	case NetworkModel.SD_RAN
		net_opt.ClassName = 'CloudAccessNetwork';
end
net_opt.AdmitPolicy = 'reject-flow';
VNF_opt.Model = VNFIntegrateModel.AllInOne;
VNF_opt.RandomSeed = [20161101 0];

%% Slice Combination
switch node_opt.Model
	case NetworkModel.Sample1
slice_config.Type = [11; 21; 31];
	case NetworkModel.SD_RAN
end
slice_config.RatioSet = {[1;1;1],[4; 1; 1],[1; 4; 1], [1; 1; 4]};
slice_config.Count = 4:4:36;

%% Initialization
experiment_length = length(slice_config.RatioSet) * length(slice_config.Count)*...
    (b_single_optimal+b_price_adjust1+b_price_adjust2+b_static_slice);
current_length = 0;
progress_bar = waitbar(current_length/experiment_length, ...
    sprintf('%d/%d',current_length, experiment_length));
progress_bar.Name = 'Computing Progress';
pause(0.01);

num_config = length(slice_config.RatioSet);
if (~exist('output_results', 'var') || ~isstruct(output_results))
    if exist('output_file_name', 'file')
        load('output_file_name', 'output_results');
    else
        output_results = struct;
    end
end
num_point = length(slice_config.Count);
num_type = length(slice_config.Type);
for exp_id = 1:num_config
    slice_config.Ratio = slice_config.RatioSet{exp_id};
    seed = 20170430;
    
    %% Statistics Initialization
    if b_single_optimal
        [stat_optimal, slice_stat_optimal] = ...
            CloudNetwork.createStatTable(num_point, num_type, 'optimal-spp'); 
    end
    if b_price_adjust1
        [stat_price1, slice_stat_price1] = ...
            CloudNetwork.createStatTable(num_point, num_type, 'dynamic-price'); 
    end
    if b_price_adjust2
        [stat_price2, slice_stat_price2] = ...
            CloudNetwork.createStatTable(num_point, num_type, 'dynamic-price'); 
    end
    if b_static_slice
        [stat_static, slice_stat_static] = ...
            CloudNetwork.createStatTable(num_point, num_type, 'static');  
        net_opt.SlicingMethod = SlicingMethod.StaticPricing;
        PN_static = instantiateclass(net_opt.ClassName, ...
            node_opt, link_opt, VNF_opt, net_opt);
        PN_static.slice_template = Slice.loadSliceTemplate(slice_config.Type);
        link_price = PN_static.getLinkCost * (1 + net_opt.PricingFactor);
        node_price = PN_static.getNodeCost * (1 + net_opt.PricingFactor);
        PN_static.writeLink('Price', link_price);
        PN_static.writeDataCenter('Price', node_price);
    end
    if b_single_optimal || b_price_adjust1 || b_price_adjust2
        PN = instantiateclass(net_opt.ClassName, ...
            node_opt, link_opt, VNF_opt, net_opt);
        PN.slice_template = Slice.loadSliceTemplate(slice_config.Type);
    end
    %%
    seed_static = seed;
    seed_dynamic = seed;
    number_slices = zeros(num_point,num_type);
    for point_id = 1:num_point
        fprintf('(%s)Iteration %d, experiment number %d.\n', datestr(now), exp_id, point_id);
        %% Specify the combination of slices
        % The combination of slices of the current point is based on that of the previous
        % point. So we should guarantee that the common part of the two point should be
        % the same, i.e. the common slices should have the same sequence and flow demand.
        total_num_slices = slice_config.Count(point_id);
        for j = 1:(num_type-1)
            ct = round(total_num_slices * slice_config.Ratio(j)/sum(slice_config.Ratio));
            if ct > total_num_slices - sum(number_slices(point_id,:))
                ct = ct - 1;
            end
            number_slices(point_id,j) = ct;
        end
        number_slices(point_id,num_type) = total_num_slices - sum(number_slices(point_id,:));
        if point_id == 1
            delta_number_slices = number_slices(point_id,:);
        else
            delta_number_slices = number_slices(point_id,:) - number_slices(point_id-1,:);
            delta_number_slices(delta_number_slices<0) = 0;
        end
        rng(seed + point_id);
        delta_total = sum(delta_number_slices);
        delta_slice_sequence = unique_randi(delta_total, delta_total, 'stable');
        delta_slice_type = zeros(delta_total,1);
        t1 = 0;
        for j = 1:num_type
            t2 = t1 + delta_number_slices(j);
            delta_slice_type((t1+1):t2) = j;
            t1 = t2;
        end
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
            rt = zeros(delta_total,1);
            for j = 1:delta_total
                s_type = delta_slice_type(delta_slice_sequence(j));
                slice_opt = PN_static.slice_template(s_type);
                slice_opt.RandomSeed = seed_static;
                seed_static = seed_static + 1;
                sl = PN_static.AddSlice(slice_opt);
                if point_id>1 && j == 1
                    stat_static{point_id, 'NumberSlicesStatic'}...
                        = stat_static{point_id-1, 'NumberSlicesStatic'};
                    stat_static{point_id, 'NumberReject'}...
                        = stat_static{point_id-1, 'NumberReject'};
                    stat_static{point_id, 'NumberPartialReject'}...
                        = stat_static{point_id-1, 'NumberPartialReject'};
                end
                if isempty(sl)
                    stat_static{point_id, 'NumberReject'}(s_type)...
                        = stat_static{point_id, 'NumberReject'}(s_type) + 1;
                    stat_static{point_id, 'NumberPartialReject'}(s_type)...
                        = stat_static{point_id, 'NumberPartialReject'}(s_type) + 1;
                else
                    if sl.NumberFlows < PN_static.slice_template(s_type).NumberFlows
                        stat_static{point_id, 'NumberPartialReject'}(s_type)...
                            = stat_static{point_id, 'NumberPartialReject'}(s_type) + 1;
                    end
                    stat_static{point_id, 'NumberSlicesStatic'}(s_type)...
                        = stat_static{point_id, 'NumberSlicesStatic'}(s_type) + 1;
                end
                fprintf('(%s)Static slicing: adding a slice %d.\n', datestr(now), ...
                    slice_opt.Type);
                tic;
                [output_static] = PN_static.staticSlicing(sl);
                rt(j) = toc;
            end
            [tb, stbs] = saveStatTable(PN_static, output_static, rt, ...
							slice_config.Type, 'static');
            stat_static(point_id, tb.Properties.VariableNames) = tb;
            for j = 1:num_type
                slice_stat_static{j}(point_id, stbs.Properties.VariableNames) = stbs(j,:);
            end
            current_length = current_length + 1;
            waitbar(current_length/experiment_length, progress_bar, ....
                sprintf('%d/%d',current_length, experiment_length));
        end
        
        if b_single_optimal || b_price_adjust1 || b_price_adjust2
            fprintf('(%s)Non-static slicing: adding all slices.\n', datestr(now));
            for j = 1:delta_total
                s_type = delta_slice_type(delta_slice_sequence(j));
                slice_opt = PN.slice_template(s_type);
                slice_opt.RandomSeed = seed_dynamic;
                seed_dynamic = seed_dynamic + 1;
                PN.AddSlice(slice_opt);
            end
        end
        %%
        if b_single_optimal
            fprintf('(%s)Global SPP.\n', datestr(now));
						PN.setOptions({'SlicingMethod', 'PricingFactor'}, ...
							{SlicingMethod.SingleNormal, net_opt.PricingFactor});
            [output_optimal, rt] = PN.singleSliceOptimization(net_opt);
            [tb, stbs] = saveStatTable(PN, output_optimal, rt, ...
							slice_config.Type, 'optimal-spp');
            stat_optimal(point_id, tb.Properties.VariableNames) = tb;
            for j = 1:num_type
                slice_stat_optimal{j}(point_id, stbs.Properties.VariableNames) = stbs(j,:);
            end
            current_length = current_length + 1;
            waitbar(current_length/experiment_length, progress_bar, ...
                sprintf('%d/%d',current_length, experiment_length));
        end
        %%
        if b_price_adjust1
            fprintf('(%s)Pricing-1.\n', datestr(now));
            PN.setOptions('SlicingMethod', SlicingMethod.AdjustPricing);
            [output_price, rt] = PN.optimizeResourcePrice();
            [tb, stbs] = saveStatTable(PN, output_price, rt, ...
							slice_config.Type, 'dynamic-price');
            stat_price1(point_id, tb.Properties.VariableNames) = tb;
            for j = 1:num_type
                slice_stat_price1{j}(point_id, stbs.Properties.VariableNames) = stbs(j,:);
            end
            current_length = current_length + 1;
            waitbar(current_length/experiment_length, progress_bar, ...
                sprintf('%d/%d',current_length, experiment_length));
        end
        %%
        if b_price_adjust2
            fprintf('(%s)Pricing-2.\n', datestr(now));
            PN.setOptions('SlicingMethod', SlicingMethod.AdjustPricing);
            [output_price, rt] = PN.optimizeResourcePriceNew();
            [tb, stbs] = saveStatTable(PN, output_price, rt, ...
							slice_config.Type, 'dynamic-price');
            stat_price2(point_id, tb.Properties.VariableNames) = tb;
            for j = 1:num_type
                slice_stat_price2{j}(point_id, stbs.Properties.VariableNames) = stbs(j,:);
            end
            current_length = current_length + 1;
            waitbar(current_length/experiment_length, progress_bar, ...
                sprintf('%d/%d',current_length, experiment_length));
        end
    end
    
    if b_single_optimal
        stat_optimal.NumberSlices = number_slices;
        output_results(exp_id).Optimal.Stat = stat_optimal; %#ok<SAGROW>
        output_results(exp_id).Optimal.SliceStat = slice_stat_optimal; %#ok<SAGROW>
    end
    if b_static_slice
        output_results(exp_id).Static.Stat = stat_static; %#ok<SAGROW>
        output_results(exp_id).Static.SliceStat = slice_stat_static; %#ok<SAGROW>
    end
    if b_price_adjust1
        stat_price1.NumberSlices = number_slices;
        output_results(exp_id).Price1.Stat = stat_price1; %#ok<SAGROW>
        output_results(exp_id).Price1.SliceStat = slice_stat_price1; %#ok<SAGROW>
    end
    if b_price_adjust2
        stat_price2.NumberSlices = number_slices;
        output_results(exp_id).Price2.Stat = stat_price2; %#ok<SAGROW>
        output_results(exp_id).Price2.SliceStat = slice_stat_price2; %#ok<SAGROW>
    end

end
close(progress_bar);
%% Figure

%% Output
% plot figures
if b_plot
	data_plot2;
end

if b_save
		switch net_opt.Model
		case NetworkModel.Sample1
			description = '2-1';
			expnum = '21';
		case NetworkModel.Sample2
			description = '2-2';
			expnum = '22';
		case NetworkModel.SD_RAN
			description = '2-3';
			expnum = '23';
	end
	description = sprintf('Experiment(No. %s): Different Consitition of Network Slices\n',...
		description);
	description = sprintf('%sNetwork Topology = %s\n', description, net_opt.Model.char);
	description = sprintf('%sThreshold = %s\n', desciption, net_opt.Threshold);
	description = sprintf('%sLink Cost Unit = %g\nNode Cost Unit%g\n', ...
		desciption, link_opt.CostUnit, node_opt.CostUnit);
	description = sprintf('%sPricingFactor = %g\n', net_opt.PricingFactor);
	if net_opt.Model == NetworkModel.SD_RAN
		cost_desc = '_LC--_NC--';
	else
		cost_desc = strcat('_LC', replace(num2str(link_opt.CostUnit), '.', ''),...
			'_NC', replace(num2str(node_opt.CostUnit), '.', ''))'
	end
	output_file_name = strcat('Results\EXP02', exp_num, cost_desc, ...
		'_PF', replace(num2str(net_opt.PricingFactor, '%02d'), '.', ''),...
		'_TH', net_opt.Threshold(1:3),...
		'.mat');
	save(output_filename, 'output_results', 'description', 'link_opt', 'node_opt', ...
		'net_opt', 'VNF_opt', 'slice_config');
end
