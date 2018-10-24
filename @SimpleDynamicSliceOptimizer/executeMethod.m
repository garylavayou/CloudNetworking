function [exitflag, fidx] = executeMethod(this, action)
global g_results event_num DEBUG INFO;  %#ok<NUSED>
this.temp_vars = struct();
this.diff_state = struct([]);
fidx = [];

%%
this.getAs_res;         % TODO, incremental upate of As_res, Hrep, and Hdiag.
this.getHrep;
this.getHdiag;
% if this.options.ReconfigMethod >= ReconfigMethod.Dimconfig
%     this.vnf_reconfig_cost = this.Parent.options.VNFReconfigCoefficient * ...
%         repmat(this.VirtualDataCenters.ReconfigCost, this.NumberVNFs, 1);
%     if this.ENABLE_DYNAMIC_NORMALIZER
%         beta = this.getbeta();
%         this.topts.vnf_reconfig_cost = beta.v*this.vnf_reconfig_cost;
%     else
%         this.topts.vnf_reconfig_cost = this.options.ReconfigScaler*this.vnf_reconfig_cost;
%     end
%     this.old_variables.v = this.VNFCapacity(:);
% end
tic
t1 = tic;
switch this.options.ReconfigMethod
    case ReconfigMethod.Baseline
        % provide 'method' and 'model' to customize the <optimalFlowRate>
        % Since we adopt |FixedCost| model, the resource cost that is a
        % constant, is not include in the objective value |profit|. So, we
        % should exclude it from the |profit|, as well as the reconfiguration
        % cost.
        options.CostModel = 'fixcost';
        options.SlicingMethod = SlicingMethod.AdjustPricing;
        this.prices.Link = this.hs.Links.Price;
        this.prices.Node = this.hs.ServiceNodes.Price;
        profit = this.optimalFlowRate(options);
    case {ReconfigMethod.Fastconfig,ReconfigMethod.FastconfigReserve}
        profit = this.fastReconfigure(action);
    case {ReconfigMethod.Fastconfig2, ReconfigMethod.Fastconfig2Reserve}
        profit = this.fastReconfigure2(action);
    case {ReconfigMethod.Dimconfig, ReconfigMethod.DimconfigReserve, ...
            ReconfigMethod.DimBaseline}
        profit = this.DimensioningReconfigure(action);
    otherwise
        error('NetworkSlicing:UnsupportedMethod', ...
            'error: unsupported method (%s) for network slicing.', ...
            this.options.ReconfigMethod.char) ;
end
t2 = toc(t1);
if this.options.ReconfigMethod>=ReconfigMethod.Dimconfig && isempty(profit)
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
stat = this.hs.get_reconfig_stat();
stat.Profit = profit;
if ~isempty(DEBUG) && DEBUG 
    disp(stat);
    fprintf('%s elapsed time: %ds\n',this.options.ReconfigMethod.char, t2);
end
switch this.options.ReconfigMethod
    case ReconfigMethod.Baseline
        stat.ReconfigType = ReconfigType.Reconfigure;
    case {ReconfigMethod.Fastconfig,ReconfigMethod.FastconfigReserve,...
            ReconfigMethod.Fastconfig2, ReconfigMethod.Fastconfig2Reserve}
        stat.ReconfigType = ReconfigType.FastReconfigure;
    case {ReconfigMethod.Dimconfig, ReconfigMethod.DimconfigReserve}
        stat.ReconfigType = ReconfigType.Dimensioning;
        if ~this.b_dim
            stat.ReconfigType = ReconfigType.FastReconfigure;
        elseif this.hs.options.Adhoc
            % If the slice do not support Adhoc flows, then we do not release resource
            % decriptors for the slice.
            this.hs.release_resource_description();
        end
    case ReconfigMethod.DimBaseline
        if this.b_dim
            if this.hs.options.Adhoc
                this.hs.release_resource_description();
            end
            stat.ReconfigType = ReconfigType.Dimensioning;
        else
            stat.ReconfigType = ReconfigType.Reconfigure;
        end
end
if this.ENABLE_DYNAMIC_NORMALIZER
    this.postl1normalizer();
    if ~isempty(DEBUG) && DEBUG
        disp(this.getbeta);
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
this.hs.net_changes = struct();
this.b_dim = false;
exitflag = 0;
end
