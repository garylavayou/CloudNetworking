function [profit,cost] = DimensioningReconfigure( this, action, new_opts )
global DEBUG; %#ok<NUSED>
if nargin <= 2
    new_opts = struct;
end
new_opts.PricingPolicy = 'quadratic-price';

this.b_dim = false;
if isempty(fieldnames(this.net_changes))
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
        this.b_dim = true;
    end
    if isfield(new_opts, 'b_dim') && new_opts.b_dim == true
        this.b_dim = true;
    end
    if this.b_dim
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
            this.topts.x_reconfig_cost = this.options.ReconfigScaler*this.x_reconfig_cost;
            this.topts.z_reconfig_cost = this.options.ReconfigScaler*this.z_reconfig_cost;
        else
            %%
            % Since we need maintain information of the deleted variables, we should use
            % the old state to calculate the information.
            this.x_reconfig_cost = ...
                (this.old_state.I_edge_path)' * this.VirtualLinks.ReconfigCost;
            old_num_paths = length(this.old_state.path_owner);
            this.z_reconfig_cost = repmat(this.VirtualDataCenters.ReconfigCost, ...
                old_num_paths*this.NumberVNFs, 1);
            this.topts.x_reconfig_cost = ...
                this.options.ReconfigScaler*this.x_reconfig_cost(~this.changed_index.x);
            this.topts.z_reconfig_cost = ...
                this.options.ReconfigScaler*this.z_reconfig_cost(~this.changed_index.z);
        end
        this.vnf_reconfig_cost = this.Parent.options.VNFReconfigCoefficient * ...
            repmat(this.VirtualDataCenters.ReconfigCost, this.NumberVNFs, 1);
        this.topts.vnf_reconfig_cost = ...
            this.options.ReconfigScaler*this.vnf_reconfig_cost;
    end
else
    this.b_dim = true;
end

if ~this.b_dim
    switch this.options.ReconfigMethod
        case {'dimconfig', 'dimconfig1'}
            [profit,cost] = this.fastReconfigure(action, new_opts);
        case {'dimconfig2'}
            [profit,cost] = this.fastReconfigure2(action, new_opts);
        case {'dimconfig0'}
            new_opts.CostModel = 'fixcost';
            new_opts.Method = 'slice-price';
            this.prices.Link = this.VirtualLinks.Price;
            this.prices.Node = this.VirtualDataCenters.Price;
            [profit,cost] = this.optimalFlowRate(new_opts);
        otherwise
    end
    return;
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
if isfield(new_opts, 'bReserve') && new_opts.bReserve
    this.lower_bounds = struct;
    this.lower_bounds.link = setlowerbounds(...
        this.old_state.link_load, this.old_state.link_capacity, 0.7);
    switch this.options.Method
        case {'dimconfig', 'dimconfig1'}
            this.lower_bounds.VNF = setlowerbounds(...
                this.old_state.vnf_load, this.old_state.vnf_capacity, 0.7);       % getVNFCapacity should be renamed as this.getVNFLoad.
        case {'dimconfig0', 'dimconfig2'}
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
