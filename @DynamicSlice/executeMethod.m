function [exitflag, fidx] = executeMethod(this, action)
global g_results event_num DEBUG INFO;  %#ok<NUSED>
this.temp_vars = struct();
this.diff_state = struct([]);
fidx = [];

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
        profit = this.optimalFlowRate(options);
    case 'fastconfig'
        profit = this.fastReconfigure(action);
    case 'fastconfig2'
        profit = this.fastReconfigure2(action);
    case {'dimconfig', 'dimconfig1', 'dimconfig2', 'dimconfig0'}
        profit = this.DimensioningReconfigure(action);
    otherwise
        error('NetworkSlicing:UnsupportedMethod', ...
            'error: unsupported method (%s) for network slicing.', ...
            this.options.ReconfigMethod) ;
end
if contains(this.options.ReconfigMethod, {'dimconfig', 'dimconfig1', 'dimconfig2'}) && isempty(profit)
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

% stat.Solution = this.Variables;
if this.ENABLE_DYNAMIC_NORMALIZER
    this.postl1normalizer();
end
stat = this.get_reconfig_stat();
stat.Profit = profit;
switch this.options.ReconfigMethod
    case 'reconfig'
        stat.ReconfigType = ReconfigType.Reconfigure;
    case {'fastconfig','fastconfig2'}
        stat.ReconfigType = ReconfigType.FastReconfigure;
    case{'dimconfig', 'dimconfig1', 'dimconfig2'}
        stat.ReconfigType = ReconfigType.Dimensioning;
        if ~this.b_dim
            stat.ReconfigType = ReconfigType.FastReconfigure;
        elseif this.getOption('Adhoc')
            % If the slice do not support Adhoc flows, then we do not release resource
            % decriptors for the slice.
            this.release_resource_description();
        end
    case {'dimconfig0'}
        if this.b_dim
            if this.getOption('Adhoc')
                this.release_resource_description();
            end
            stat.ReconfigType = ReconfigType.Dimensioning;
        else
            stat.ReconfigType = ReconfigType.Reconfigure;
        end
end
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
this.b_dim = false;
exitflag = 0;
end
