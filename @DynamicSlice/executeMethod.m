function [exitflag, fidx] = executeMethod(this, action)
global g_results event_num DEBUG InfoLevel; %#ok<NUSED>
this.temp_vars = struct;

%%
options = getstructfields(this.Parent.options, {'Form'});  % Form = {'compact'|'normal'}
options.PricingPolicy = 'quadratic-price';
this.getAs_res;         % TODO, incremental upate of As_res, Hrep, and Hdiag.
this.getHrep;
this.getHdiag;
if ~strcmpi(this.options.Method, 'fastconfig')
    this.vnf_reconfig_cost = this.Parent.options.VNFReconfigCoefficient * ...
        repmat(this.VirtualDataCenters.ReconfigCost, this.NumberVNFs, 1);
    this.topts.vnf_reconfig_cost = this.options.ReconfigScaler*this.vnf_reconfig_cost;
    %% Save the VNF capacity to the previous state, berfore optimization
    % After the slice is created, |VNFCapacity| is recorded. After each optimization, the
    % |VNFCapacity| is updated.
    % When perform 'dimconfig' or 'dimconfig2', if the topology is augmented with new data
    % centers, |this.old_variables.v| will be modified to match the size of the new
    % vector (See <DimensioningReconfigure>).
    this.old_variables.v = this.VNFCapacity(:);
end
switch this.options.Method
    case 'reconfig'
        % provide 'method' and 'model' to customize the <optimalFlowRate>
        % Since we adopt |FixedCost| model, the resource cost that is a
        % constant, is not include in the objective value |profit|. So, we
        % should exclude it from the |profit|, as well as the reconfiguration
        % cost.
        options.CostModel = 'fixcost';
        options.Method = 'slice-price';
        this.prices.Link = this.VirtualLinks.Price;
        this.prices.Node = this.VirtualDataCenters.Price;
        [profit,cost] = this.optimalFlowRate(options);
    case 'fastconfig'
        [profit,cost] = this.fastReconfigure(action, options);
    case 'fastconfig2'
        [profit,cost] = this.fastReconfigure2(action, options);
    case {'dimconfig', 'dimconfig1', 'dimconfig2', 'dimconfig0'}
        options.bReserve = true;  % resource reservation
        [profit,cost] = this.DimensioningReconfigure(action, options);
    otherwise
        error('NetworkSlicing:UnsupportedMethod', ...
            'error: unsupported method (%s) for network slicing.', ...
            this.options.Method) ;
end
if contains(this.options.Method, {'dimconfig', 'dimconfig2'}) && isempty(profit)
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
fidx = find(this.FlowTable.Rate<=0.1*median(this.FlowTable.Rate));
if isempty(fidx)
    exitflag = 1;
    if nargout >= 2
        fidx = [];
    end
else
    exitflag = 0;
    fidx0 = this.FlowTable.Rate<=0;
    if InfoLevel.UserModel>=DisplayLevel.Notify ||...
            InfoLevel.UserModelDebug>=DisplayLevel.Final
        if ~isempty(setdiff(this.FlowTable.Identifier(fidx0), this.reject_index))
            cprintf('SystemCommands', '[%s]: reject flows due to zero rate.\n', this.options.Method);
        end
    end
    % due to reconfiguration cost, some arriving flows might be directly rejected. In this
    % case we need to perform slice redimensioning.
    if contains(this.options.Method, {'dimconfig'}) && ...
            ~isempty(setdiff(this.FlowTable.Identifier(fidx), this.reject_index))
        cprintf('SystemCommands', '[%s]: reconfigure flows due to low rate.\n', this.options.Method);
        options.b_dim = true;
        [profit,cost] = this.DimensioningReconfigure(action, options);
    end
end
this.reject_index = this.FlowTable.Identifier(fidx);

% stat.Solution = this.Variables;
stat = this.get_reconfig_stat();
switch this.options.Method
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
