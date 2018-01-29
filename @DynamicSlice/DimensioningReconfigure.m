function [profit,cost, exitflag, fidx] = DimensioningReconfigure( this, action, new_opts )
global DEBUG INFO; 
if nargin <= 2
    new_opts = struct;
end
new_opts.PricingPolicy = 'quadratic-price';
new_opts = structmerge(new_opts, ...
    getstructfields(this.Parent.options, 'PricingPolicy', 'ignore'));

%% 
[omega, sigma_o] = this.utilizationRatio();
this.sh_data.index = this.sh_data.index + 1;
a = this.sh_options.alpha;
TL = this.sh_options.trend_length;
SL = this.sh_options.series_length;
omega_new = a*this.sh_data.omegas(end) + (1-a)*omega;
if length(this.sh_data.omegas) < SL
    this.sh_data.omegas(end+1) = omega_new;
else
    this.sh_data.omegas = [this.sh_data.omegas(2:end), omega_new];
end
trend = {'ascend','descend'};
for i = 1:2
    this.sh_data.omega_trend.(trend{i})= ...
        update_trend(this.sh_data.omega_trend.(trend{i}), ...
        struct('value', omega_new, 'index', this.sh_data.index), ...
        struct('trend', trend{i}, 'length', SL));
end
this.b_dim = false;
this.options.Reserve = 1;
if omega > this.sh_options.omega_upper && ...
        (length(this.sh_data.omega_trend.ascend)>TL ||...
        length(this.sh_data.profit_trend.ascend)>TL)
    this.b_dim = true;
    this.options.Reserve = this.sh_options.omega_lower+0.1;
elseif omega < this.sh_options.omega_lower && ...
        (length(this.sh_data.omega_trend.descend)>TL ||...
        length(this.sh_data.profit_trend.descend)>TL)
    this.b_dim = true;
    this.options.Reserve = this.sh_options.omega_upper-0.05;
else
    this.options.Reserve = omega;
    if sigma_o > 0.25 || ~isempty(fieldnames(this.net_changes))
        this.b_dim = true;
    else 
        switch this.getOption('Trigger')
            case 'TimeBased'
                if this.time.Current - this.time.LastDimensioning >= this.options.TimeInterval
                    this.b_dim = true;
                end
            case 'EventBased'
                if mod(this.event.RecentCount, this.options.EventInterval) == 0
                    this.b_dim = true;
                end
            case 'ProfitBased'
            otherwise
        end
        if this.NumberFlows > 0 && sum(this.VirtualLinks.Capacity) == 0
            % If the slice has no resource but new flows arrive, perform dimensioning at once.
            this.b_dim = true;
        end
        if this.NumberFlows == 0
            % the dimension method will handle the idl slice, see also
            % <DynamicCloudNetwork.optimizeResourcePriceNew>
            this.b_dim = true;
        end
        if isfield(new_opts, 'b_dim') && new_opts.b_dim == true
            this.b_dim = true;
        end
        if this.b_dim
            update_reconfig_cost;
        end
    end
end

if ~this.b_dim
    switch this.options.ReconfigMethod
        case {'dimconfig', 'dimconfig1'}
            [profit,cost] = this.fastReconfigure(action, new_opts);
        case {'dimconfig2'}
            [profit,cost] = this.fastReconfigure2(action, new_opts);
        case {'dimconfig0'}
            new_opts.CostModel = 'fixcost';
            new_opts.SlicingMethod = 'slice-price';
            this.prices.Link = this.VirtualLinks.Price;
            this.prices.Node = this.VirtualDataCenters.Price;
            [profit,cost] = this.optimalFlowRate(new_opts);
            return;
        otherwise
    end
    
    fidx = find(this.FlowTable.Rate<=0.1*median(this.FlowTable.Rate));
    if isempty(fidx)
        exitflag = 1;
        this.sh_data.reserve_dev = this.sh_data.reserve_dev/2;
    else
        exitflag = 0;
        fidx0 = this.FlowTable.Rate<=0;
        if ~isempty(setdiff(this.FlowTable.Identifier(fidx0), this.reject_index))
            message = sprintf('%s: [%s] reject flows due to zero rate.', ...
                calledby, this.options.ReconfigMethod);
            if ~isempty(DEBUG) && DEBUG
                warning(message); %#ok<SPWRN>
            else
                if ~isempty(INFO) && INFO
                    cprintf('SystemCommands', '%s\n', message);
                end
            end
        end
        % due to reconfiguration cost, some arriving flows might be directly rejected. In this
        % case we need to perform slice redimensioning.
        if ~isempty(setdiff(this.FlowTable.Identifier(fidx), this.reject_index))
            message = sprintf('%s: [%s] reconfigure flows due to low rate.\n', ...
                calledby, this.options.ReconfigMethod);
            if ~isempty(DEBUG) && DEBUG
                warning(message); %#ok<SPWRN>
            else
                if ~isempty(INFO) && INFO
                    cprintf('SystemCommands', '%s\n', message);
                end
            end
            this.b_dim = true;
            update_reconfig_cost;
            this.sh_data.reserve_dev = this.sh_data.reserve_dev + 1;
            this.options.Reserve = this.options.Reserve - this.sh_data.reserve_dev*0.02;
        end
    end
    this.reject_index = this.FlowTable.Identifier(fidx);
end
%% TODO
% update vnf_reconfig_cost, since the number of nodes has changed.
% To make sure that the 'Load' field of physical links and nodes is updated.
% ISSUE-1: If a resource has no load, then the reconfiguration cost is zero.
% ISSUE-2: The reconfiguration cost when dimensioning is based on total system load, while
%          intra-slice reconfiguration cost is based on slice load.
% After the dimensioning process, the reconfiguration cost will be updated in _finalize_.
%% 
% As the link and node reconfiguration cost is changed, we should update the coefficient.

%% Perform dimensioning
% resource pricing with consideration on reconfiguration cost.
% updated price only applied to those involved slices. (If applied to other slices,
% more reconfiguration might be trigger, since their profit changes due to price changes.) 

% Update old variables, see also <OnAddingFlow>.
% The corresponding reconfiguration cost has been updated when update the adhoc flows.see also
% <OnAddingFlow> and <OnRemovingFlow>. 
if this.b_dim
    if isfield(this.changed_index, 'v')
        % vnf instance will change only when adding new adhoc flows.
        % when removing flow, the vnf instance might be released after slice dimensioning
        num_vnfs = length(this.changed_index.v);
        this.old_variables.v = zeros(num_vnfs,1);
        this.old_variables.v(~this.changed_index.v) = this.Variables.v;
        % this.changed_index.v = [];  % should be reset, check if there are other uses.
    else
        this.old_variables.v = this.Variables.v;
    end
    %% Resource Reservation
    % [optional function]: when we perform slice dimensioning, the slice will try to release
    %   all idle resources to lower the operational cost. However, certain amount of redundant
    %   resources can accomodate the traffic dynamics and reduce the reconfiguration cost.
    %   Therefore, we should not release all idle resources at dimensioning, unless the
    %   utilization of a resource is lower than the specified threshold.
    %
    % We will calculate the lower bound of resource reservation, given the threshold 'th'.
    %   1. If a resource's utilization ratio 'u=0': the low-bound as c/2;
    %   2. Else if a resource's utilization ratio 'u' is low than 'th':
    %       we set the lower-bound as 'uc/[(1+th)/2]';
    %   3. Otherwise, the lower-bound is 'c'.
    %
    % * Computing Resource Reservation *
    % With 'dimconfig', we cannot reconfigure the VNF instance capacity after dimensioning.
    % Therefore, we need to reserve resource for VNF instances. On the other hand, with
    % 'dimconfig2' and 'dimconfig0', we only need reserved resources for the nodes.
    %% TODO differentiate dimconfig0 from other reconfiguration methods
    if isfield(this.options, 'bReserve') && this.options.bReserve
        this.lower_bounds = struct;
        this.lower_bounds.link = setlowerbounds(...
            this.old_state.link_load, this.old_state.link_capacity, 0.7);
        switch this.options.ReconfigMethod
            case {'dimconfig1'}
                this.lower_bounds.VNF = setlowerbounds(...
                    this.old_state.vnf_load, this.old_state.vnf_capacity, 0.7);       % getVNFCapacity should be renamed as this.getVNFLoad.
            case {'dimconfig', 'dimconfig0', 'dimconfig2'}
                this.lower_bounds.node = setlowerbounds(...
                    this.old_state.node_load, this.old_state.node_capacity, 0.7);
        end
    else
        this.lower_bounds = struct([]);
    end
    % dimensioning will be perform by <DynamicCloudNetwork>
    this.b_derive_vnf = false;      % VNF capapcity is derived from variables v, instead of z.
    notify(this, 'RequestDimensioning', EventData());
    
    if this.Results.Value == 0
        profit = this.Results.Profit(1);
    else
        profit = [];
    end
    if nargout>=2
        if this.Results.Value == 0
            cost = this.getSliceCost(new_opts.PricingPolicy);
        else
            cost =[];
        end
    end
    %% Reset statistics
    % update the period for performing dimensioning by exponential moving average
    this.time.DimensionInterval = this.time.DimensionInterval * this.a + (this.time.Current-this.time.LastDimensioning)*(1-this.a);
    this.time.LastDimensioning = this.time.Current;
    this.event.RecentCount = 0;
    this.lower_bounds = struct([]);
end

%% update trend
profit_new = a*this.sh_data.profits(end) + (1-a)*profit;
if length(this.sh_data.profits) < this.sh_options.series_length
    this.sh_data.profits(end+1) = profit_new;
else
    this.sh_data.profits = [this.sh_data.profits(2:end), profit_new];
end
type = {'omega_trend', 'profit_trend'};
for i = 1:2
    this.sh_data.profit_trend.(trend{i})= ...
        update_trend(this.sh_data.profit_trend.(trend{i}), ...
        struct('value', profit_new, 'index', this.sh_data.index), ...
        struct('trend', trend{i}, 'length', this.sh_options.series_length));
end
if this.sh_data.index == intmax
    idx_offset = min([this.sh_data.omega_trend.ascend(1).index,...
        this.sh_data.omega_trend.descend(1).index,...
        this.sh_data.profit_trend.ascend(1).index,...
        this.sh_data.profit_trend.desend(1).index]);
    this.sh_data.index = 0;
    for i = 1:2
        for j = 1:2
            [this.sh_data.(type{j}).(trend{i})] = ...
                update_trend(this.sh_data.(type{j}).(trend{i}), ...
                idx_offset, struct('bResetIndex', true));
            this.sh_data.index = max(this.sh_data.index,...
                this.sh_data.(type{j}).(trend{i})(end).index);
        end
    end
end
if this.b_dim
    for i = 1:2
        for j = 1:2
            if this.options.Reserve ~= omega
                [this.sh_data.(type{j}).(trend{i})] = this.sh_data.(type{j}).(trend{i})(end);
            else
                [this.sh_data.(type{j}).(trend{i})] = ...
                    update_trend(this.sh_data.(type{j}).(trend{i}), ...
                    this.sh_options.trend_length, struct('bPurge', true));
            end
        end
    end
end

%% inline
    function update_reconfig_cost()
        %% Period re-dimensioing
        % The number of virtual nodes/links is not change in the slice, as well as the
        % number of VNF instances.
        % When performing re-dimensioning, the reconfiguration cost is larger than that of
        % fast reconfiguration, so we need to update the reconfiguration cost.
        this.Parent.updateRedimensionCost(this);
        if strcmpi(action, 'add')
            this.x_reconfig_cost = (this.I_edge_path)' * this.VirtualLinks.ReconfigCost;
            this.z_reconfig_cost = repmat(this.VirtualDataCenters.ReconfigCost, ...
                this.NumberPaths*this.NumberVNFs, 1);
            if this.ENABLE_DYNAMIC_NORMALIZER
                beta =this.getbeta();
                this.topts.x_reconfig_cost = beta.x*this.x_reconfig_cost;
                this.topts.z_reconfig_cost = beta.z*this.z_reconfig_cost;
            else
                this.topts.x_reconfig_cost = this.options.ReconfigScaler*this.x_reconfig_cost;
                this.topts.z_reconfig_cost = this.options.ReconfigScaler*this.z_reconfig_cost;
            end
        else
            %%
            % Since we need maintain information of the deleted variables, we should use
            % the old state to calculate the information.
            this.x_reconfig_cost = ...
                (this.old_state.I_edge_path)' * this.VirtualLinks.ReconfigCost;
            old_num_paths = length(this.old_state.path_owner);
            this.z_reconfig_cost = repmat(this.VirtualDataCenters.ReconfigCost, ...
                old_num_paths*this.NumberVNFs, 1);
            if this.ENABLE_DYNAMIC_NORMALIZER
                beta =this.getbeta();
                this.topts.x_reconfig_cost = ...
                    beta.x*this.x_reconfig_cost(~this.changed_index.x);
                this.topts.z_reconfig_cost = ...
                    beta.z*this.z_reconfig_cost(~this.changed_index.z);
            else
                this.topts.x_reconfig_cost = ...
                    this.options.ReconfigScaler*this.x_reconfig_cost(~this.changed_index.x);
                this.topts.z_reconfig_cost = ...
                    this.options.ReconfigScaler*this.z_reconfig_cost(~this.changed_index.z);
            end
        end
        this.vnf_reconfig_cost = this.Parent.options.VNFReconfigCoefficient * ...
            repmat(this.VirtualDataCenters.ReconfigCost, this.NumberVNFs, 1);
        if this.ENABLE_DYNAMIC_NORMALIZER
            beta = this.getbeta();
            this.topts.vnf_reconfig_cost = beta.v*this.vnf_reconfig_cost;
        else
            this.topts.vnf_reconfig_cost = ...
                this.options.ReconfigScaler*this.vnf_reconfig_cost;
        end
    end

end

function lbs = setlowerbounds(load, capacity, threshold)
    lbs = zeros(size(load));
    idx = find(capacity~=0);
    util = load(idx)./capacity(idx);
    idx1 = idx(util<=threshold/2);
    idx2 = idx((util>threshold/2) & (util<=threshold));
    idx3 = idx(util>threshold);
    lbs(idx1) = capacity(idx1)/2;
    lbs(idx2) = load(idx2)/(threshold);
    lbs(idx3) = capacity(idx3);
end
