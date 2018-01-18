function [exitflag, fidx] = executeMethod(this, action)
global g_results event_num DEBUG INFO; 
this.temp_vars = struct;

%%
this.getAs_res;         % TODO, incremental upate of As_res, Hrep, and Hdiag.
this.getHrep;
this.getHdiag;
if ~strcmpi(this.options.ReconfigMethod, 'fastconfig')
    this.vnf_reconfig_cost = this.Parent.options.VNFReconfigCoefficient * ...
        repmat(this.VirtualDataCenters.ReconfigCost, this.NumberVNFs, 1);
    if this.ENABLE_DYNAMIC_NORMALIZER
        beta = this.getbeta();
        this.topts.vnf_reconfig_cost = beta.v*this.vnf_reconfig_cost;
    else
        this.topts.vnf_reconfig_cost = this.options.ReconfigScaler*this.vnf_reconfig_cost;
    end
    %% Save the VNF capacity to the previous state, berfore optimization
    % After the slice is created, |VNFCapacity| is recorded. After each optimization, the
    % |VNFCapacity| is updated.
    % When perform 'dimconfig' or 'dimconfig2', if the topology is augmented with new data
    % centers, |this.old_variables.v| will be modified to match the size of the new
    % vector (See <DimensioningReconfigure>).
    this.old_variables.v = this.VNFCapacity(:);
end
switch this.options.ReconfigMethod
    case 'reconfig'
        % provide 'method' and 'model' to customize the <optimalFlowRate>
        % Since we adopt |FixedCost| model, the resource cost that is a
        % constant, is not include in the objective value |profit|. So, we
        % should exclude it from the |profit|, as well as the reconfiguration
        % cost.
        options.CostModel = 'fixcost';
        options.SlicingMethod = 'slice-price';
        this.prices.Link = this.VirtualLinks.Price;
        this.prices.Node = this.VirtualDataCenters.Price;
        [profit,cost] = this.optimalFlowRate(options);
    case 'fastconfig'
        [profit,cost] = this.fastReconfigure(action);
    case 'fastconfig2'
        [profit,cost] = this.fastReconfigure2(action);
    case {'dimconfig', 'dimconfig1', 'dimconfig2', 'dimconfig0'}
        [profit,cost] = this.DimensioningReconfigure(action);
    otherwise
        error('NetworkSlicing:UnsupportedMethod', ...
            'error: unsupported method (%s) for network slicing.', ...
            this.options.ReconfigMethod) ;
end
if contains(this.options.ReconfigMethod, {'dimconfig', 'dimconfig2'}) && isempty(profit)
    exitflag = -1;
    return;
end
% Check resource utilization
%{
idx = this.VirtualLinks.Capacity~=0;
link_util = this.VirtualLinks.Load(idx)./this.VirtualLinks.Capacity(idx);
sum_link_util = sum(this.VirtualLinks.Load(idx))./sum(this.VirtualLinks.Capacity(idx));
mean_link_util = mean(link_util);
%}
old_info = INFO;
INFO = true;
fidx = find(this.FlowTable.Rate<=0.1*median(this.FlowTable.Rate));
if isempty(fidx)
    exitflag = 1;
    if nargout >= 2
        fidx = [];
    end
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
    if contains(this.options.ReconfigMethod, {'dimconfig'}) && ...
            ~isempty(setdiff(this.FlowTable.Identifier(fidx), this.reject_index))
        message = sprintf('%s: [%s] reconfigure flows due to low rate.\n', ...
            calledby, this.options.ReconfigMethod);
        if ~isempty(DEBUG) && DEBUG
            warning(message); %#ok<SPWRN>
        else
            if ~isempty(INFO) && INFO
                cprintf('SystemCommands', '%s\n', message);
            end
        end
        options.b_dim = true;
        [profit,cost] = this.DimensioningReconfigure(action, options);
    end
end
this.reject_index = this.FlowTable.Identifier(fidx);
INFO = old_info;

% stat.Solution = this.Variables;
if this.ENABLE_DYNAMIC_NORMALIZER
    [stat, reconfig_cost] = this.get_reconfig_stat();
    if this.b_dim
        this.postl1normalizer(reconfig_cost, reconfig_cost.linear);
    else
        this.postl1normalizer(getstructfields(reconfig_cost, {'x','z'}),...
            getstructfields(reconfig_cost.linear, {'x','z'}));
    end
else
    stat = this.get_reconfig_stat();
end
switch this.options.ReconfigMethod
    case 'reconfig'
        stat.Profit = profit - stat.Cost;
        stat.ReconfigType = ReconfigType.Reconfigure;
    case {'fastconfig','fastconfig2'}
        stat.Profit = profit + stat.LinearCost - stat.Cost;
        stat.ReconfigType = ReconfigType.FastReconfigure;
    case{'dimconfig', 'dimconfig1', 'dimconfig2'}
        %% ISSUE: LinearCost does not include in the profit
        stat.Profit = profit - stat.Cost;      %  + stat.LinearCost 
        stat.ReconfigType = ReconfigType.Dimensioning;
        if ~this.b_dim
            stat.ReconfigType = ReconfigType.FastReconfigure;
            stat.Profit = stat.Profit + stat.LinearCost;
        elseif this.getOption('Adhoc')
            % If the slice do not support Adhoc flows, then we do not release resource
            % decriptors for the slice.
            this.release_resource_description();
        end
    case {'dimconfig0'}
        if this.b_dim
            stat.Profit = profit - stat.Cost;   % + stat.LinearCost 
            if this.getOption('Adhoc')
                this.release_resource_description();
            end
            stat.ReconfigType = ReconfigType.Dimensioning;
        else
            stat.Profit = profit - stat.Cost;
            stat.ReconfigType = ReconfigType.Reconfigure;
        end
end
stat.ResourceCost = cost;
stat.FairIndex = (sum(this.FlowTable.Rate))^2/(this.NumberFlows*sum(this.FlowTable.Rate.^2));
g_results(event_num,:) = stat;

%% Reset temporary variables
% These variables should be cleared. If it is used next time, it will be assigned with new
% values. Some methods will execute according to the state (if it is initialized) of the
% variables.
this.prices.Link = [];
this.prices.Node = [];
this.topts = struct();
this.changed_index = struct();
this.net_changes = struct();

end
