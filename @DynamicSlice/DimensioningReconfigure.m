function profit = DimensioningReconfigure( this, action, new_opts )
global DEBUG; %#ok<NUSED>
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
    if this.b_dim
        %% Period re-dimensioing
        % The number of virtual nodes/links is not change in the slice, as well as the number of VNF
        % instances.
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
    if strcmpi('dimconfig', this.options.Method)
        profit = this.fastReconfigure(action, new_opts);
    elseif strcmpi('dimconfig2', this.options.Method)
        profit = this.fastReconfigure2(action, new_opts);
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

% dimensioning will be perform by <DynamicCloudNetwork>
this.b_derive_vnf = false;      
notify(this, 'RequestDimensioning', EventData());

if this.Results.Value == 0
    profit = this.Results.Profit(1);
else
    profit = [];
end

%% Reset statistics
this.time.LastDimensioning = this.time.Current;
this.event.RecentCount = 0;
end
