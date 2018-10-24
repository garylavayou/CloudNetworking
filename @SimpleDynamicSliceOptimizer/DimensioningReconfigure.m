function [profit,cost, exitflag, fidx] = DimensioningReconfigure( this, action, new_opts )
global DEBUG INFO; 
if nargin <= 2
    new_opts = struct;
end
new_opts.PricingPolicy = 'quadratic-price';
new_opts = structmerge(new_opts, ...
    getstructfields(this.hs.Parent.options, 'PricingPolicy', 'ignore'));

%% Slice Dimensioning Scheduling
[omega, sigma_o] = this.hs.utilizationRatio();
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
if this.options.bReserve
    %% How to set the current utilization constraint
    % General principles applied to 'remove' flow: the target utilization ratio shoull be
    % greater than the current value. The greater value of 'omega' could help utilize
    % redundant resources.
    % (a) conservative: omega;
    % (b) bound-trend: max(omega, omega_new);
    % (c) bound: a*omega+(1-a)*omega_upper;
    % For 'add' flow: it is dependent on the fast reconfiguration method how to reserve idle
    % resources.
    omega_m = mean([this.sh_options.omega_upper,this.sh_options.omega_lower]);
    if omega < omega_m
        % when resource utilization ratio is low, set the target ratio higher to use idle
        % resources.
        this.options.Reserve = omega_m + (omega_m - omega);
    else
        this.options.Reserve = omega;
    end
end
if isempty(fieldnames(this.hs.net_changes))
    this.b_dim = false;
    if omega > this.sh_options.omega_upper && ...
            (length(this.sh_data.omega_trend.ascend)>TL ||...
            length(this.sh_data.profit_trend.ascend)>TL)
        this.b_dim = true;
        if this.options.bReserve
            this.options.Reserve = 0.5*this.sh_options.omega_lower+...
                0.5*this.sh_options.omega_upper;     % reserve resources for future traffic rising
        end
    elseif omega < this.sh_options.omega_lower && ...
            (length(this.sh_data.omega_trend.descend)>TL ||...
            length(this.sh_data.profit_trend.descend)>TL)
        this.b_dim = true;
        if this.options.bReserve
            this.options.Reserve = 0.75*this.sh_options.omega_upper+...
                0.25*this.sh_options.omega_lower;     % release redundant resources due to traffic decline
        end
    else
        if~isempty(DEBUG) && DEBUG
            disp([sigma_o this.sh_options.sigma]);
        end
        if sigma_o > this.sh_options.sigma
            this.b_dim = true;
        end
    end
    if ~this.b_dim
        switch this.options.Trigger
            case 'TimeBased'
                if this.time.Current - this.hs.time.LastDimensioning >= this.options.TimeInterval
                    this.b_dim = true;
                end
            case 'EventBased'
                if mod(this.hs.event.RecentCount, this.options.EventInterval) == 0
                    this.b_dim = true;
                end
            case 'ProfitBased'
            otherwise
        end
        if this.hs.NumberFlows > 0 && sum(this.hs.Links.Capacity) == 0
            % If the slice has no resource but new flows arrive, perform dimensioning at once.
            this.b_dim = true;
        end
        if this.hs.NumberFlows == 0
            % the dimension method will handle the idl slice, see also
            % <DynamicCloudNetwork.optimizeResourcePriceNew>
            this.b_dim = true;
        end
    end
%     if this.b_dim
%         this.update_reconfig_costinfo(action, true);
%     end
else
    this.b_dim = true;
end
if ~this.b_dim 
    switch this.options.ReconfigMethod
        case {ReconfigMethod.Dimconfig, ReconfigMethod.DimconfigReserve, ReconfigMethod.Dimconfig1}
            [profit,cost] = this.fastReconfigure(action, new_opts);
        case ReconfigMethod.Dimconfig2
            [profit,cost] = this.fastReconfigure2(action, new_opts);
        case ReconfigMethod.DimBaseline
            new_opts.CostModel = 'fixcost';
            new_opts.SlicingMethod = SlicingMethod.AdjustPricing;
            this.prices.Link = this.hs.Links.Price;
            this.prices.Node = this.hs.ServiceNodes.Price;
            [profit,cost] = this.optimalFlowRate(new_opts);
        otherwise
            error('%s: invalid reconfiguration method.', calledby);
    end
    
    fidx = find(this.hs.FlowTable.Rate<=0.1*median(this.hs.FlowTable.Rate));
    if isempty(fidx)
        exitflag = 1;
        this.sh_data.reserve_dev = this.sh_data.reserve_dev/2;
    else
        exitflag = 0;
        fidx0 = this.hs.FlowTable.Rate<=0;
        if ~isempty(setdiff(this.hs.FlowTable.Identifier(fidx0), this.hs.reject_index))
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
        if ~isempty(setdiff(this.hs.FlowTable.Identifier(fidx), this.hs.reject_index))
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
%             this.update_reconfig_costinfo(action, true);
            this.temp_vars = struct();
            this.diff_state = struct([]);
            this.prices.Link = [];
            this.prices.Node = [];      
            if this.options.bReserve
                this.sh_data.reserve_dev = this.sh_data.reserve_dev + 1;
                this.options.Reserve = this.options.Reserve - this.sh_data.reserve_dev*0.02;
            end
        end
    end
    this.hs.reject_index = this.hs.FlowTable.Identifier(fidx);
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
        % If the topology is augmented with new data centers, |this.old_variables.v| will be
        % modified to match the size of the new vector.
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
    % [optional function]: when we perform slice dimensioning, the slice will try to
    % release all idle resources to lower the operational cost. However, certain amount of
    % redundant resources can accomodate the traffic dynamics and reduce the
    % reconfiguration cost. 
    %   Therefore, we should not release all idle resources at dimensioning, unless the
    %   utilization of a resource is lower than the specified threshold.
    %
    % We will calculate the lower bound of resource reservation, given the threshold 'th'.
    %   1. If a resource's utilization ratio 'u<=th/2': the low-bound as c/2;
    %   2. Else if a resource's utilization ratio 'u' is low than 'th':
    %       we set the lower-bound as 'uc/th';
    %   3. Otherwise, the lower-bound is 'c'.
    % We set 'th=(omega_upper+omega_lower)/2', that is, 
    %   If a single resource utilization is higher than the overall utilization, we
    %   reserve all its current resources and more capacity might be added to this
    %   resource;   
    %   If the utilization is lower than 'omega' but higher than 'omega/2', then partial
    %   capacity could be released.
    %   Otherwise, if the utilization is lower than 'omega/2', we half its capacity (not
    %   so rapidly release as (2)).
    % * Computing Resource Reservation *
    % With 'dimconfig1', we cannot reconfigure the VNF instance capacity after dimensioning.
    % Therefore, we need to reserve resource for VNF instances. On the other hand, with
    % 'dimconfig2' and 'dimconfig0', we only need reserved resources for the nodes.
    % With 'DimconfigReserve', can we both reserve for individual VNF instances or nodes?
    %% TODO differentiate dimconfig0 from other reconfiguration methods
    if this.options.bReserve
        this.lower_bounds = struct;
        this.upper_bounds = struct;
        [this.lower_bounds.link,  this.upper_bounds.link] = setbounds(...
            this.hs.old_state.link_load, this.hs.old_state.link_capacity, this.options.Reserve);
        switch this.options.ReconfigMethod
            case ReconfigMethod.Dimconfig1
                [this.lower_bounds.VNF, this.upper_bounds.VNF] = setbounds(...
                    this.hs.old_state.vnf_load, this.old_state.vnf_capacity, this.options.Reserve);       % getVNFCapacity should be renamed as this.getVNFLoad.
            case {ReconfigMethod.DimconfigReserve, ReconfigMethod.Dimconfig, ReconfigMethod.Dimconfig2}
                if this.options.bReserve == 3
                    [this.lower_bounds.VNF, this.upper_bounds.VNF] = setbounds(...
                        this.hs.old_state.vnf_load, this.old_state.vnf_capacity, this.options.Reserve);       % getVNFCapacity should be renamed as this.getVNFLoad.
                else
                    [this.lower_bounds.node, this.upper_bounds.node] = setbounds(...
                        this.hs.old_state.node_load, this.hs.old_state.node_capacity, this.options.Reserve);
                end
        end
    else
        this.lower_bounds = struct([]);
        this.upper_bounds = struct([]);
    end
    % dimensioning will be perform by <DynamicCloudNetwork>
    %     if isfield(this.options, 'Reserve') && this.options.Reserve>0
    %         this.options.Reserve = -this.options.Reserve;
    %     end
    notify(this.hs, 'RequestDimensioning', EventData());
    
    if this.hs.Results.Value == 0
        profit = this.hs.Results.Profit(1);
    else
        profit = [];
    end
    if nargout>=2
        if this.hs.Results.Value == 0
            cost = this.hs.getCost(new_opts.PricingPolicy);
        else
            cost =[];
        end
    end
    %% Reset statistics
    % update the period for performing dimensioning by exponential moving average
    new_dim_interval = this.hs.time.DimensionInterval * this.a + ...
        (this.hs.time.Current-this.hs.time.LastDimensioning)*(1-this.a);
    this.hs.time.DimensionInterval = max(this.hs.time.MinDimensionInterval, new_dim_interval);
    this.hs.time.LastDimensioning = this.hs.time.Current;
    this.hs.event.RecentCount = 0;
    this.lower_bounds = struct([]);
end

%% update trend
profit_new = a*this.sh_data.profits(end) + (1-a)*profit;
if length(this.sh_data.profits) < SL
    this.sh_data.profits(end+1) = profit_new;
else
    this.sh_data.profits = [this.sh_data.profits(2:end), profit_new];
end
type = {'omega_trend', 'profit_trend'};
for i = 1:2
    this.sh_data.profit_trend.(trend{i})= ...
        update_trend(this.sh_data.profit_trend.(trend{i}), ...
        struct('value', profit_new, 'index', this.sh_data.index), ...
        struct('trend', trend{i}, 'length', SL));
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
omega_new = this.hs.utilizationRatio();
if this.b_dim
    for i = 1:2
        for j = 1:2
            if abs(omega_new-omega) < 0.05
                [this.sh_data.(type{j}).(trend{i})] = ...
                    update_trend(this.sh_data.(type{j}).(trend{i}), ...
                    TL, struct('bPurge', true));
            else
                [this.sh_data.(type{j}).(trend{i})] = this.sh_data.(type{j}).(trend{i})(end);
            end
        end
    end
end

end


function [lbs, ubs] = setbounds(load, capacity, threshold)
    lbs = zeros(size(load));
    ubs = zeros(size(load));
    idx = find(capacity~=0);
    util = load(idx)./capacity(idx);
    idx1 = idx(util<=threshold/2);
    idx2 = idx((util>threshold/2) & (util<=threshold));
    idx3 = idx(util>threshold);
    lbs(idx1) = capacity(idx1)/2;
    lbs(idx2) = load(idx2)/(threshold);
    lbs(idx3) = capacity(idx3);
    ubs(idx1) = (capacity(idx1)+lbs(idx1))/2;
end
