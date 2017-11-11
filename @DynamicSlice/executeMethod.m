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
% stat.Solution = this.Variables;
stat = this.get_reconfig_stat();
switch this.options.Method
    case 'reconfig'
        stat.Profit = profit - stat.Cost;
    case {'fastconfig','fastconfig2'}
        stat.Profit = profit + stat.LinearCost - stat.Cost;
    case{'dimconfig', 'dimconfig1', 'dimconfig2'}
        stat.Profit = profit + stat.LinearCost - stat.Cost;
        if this.b_dim && this.getOption('Adhoc')
            % If the slice do not support Adhoc flows, then we do not release resource
            % decriptors for the slice.
            this.release_resource_description();
        end
    case {'dimconfig0'}
        if this.b_dim
            stat.Profit = profit + stat.LinearCost - stat.Cost;
            if this.getOption('Adhoc')
                this.release_resource_description();
            end
        else
            stat.Profit = profit - stat.Cost;
        end
end
stat.ResourceCost = cost;
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

fidx = find(this.FlowTable.Rate<=0);
if isempty(fidx)
    exitflag = 1;
    if nargout >= 2
        fidx = [];
    end
else
    exitflag = 0;
    if InfoLevel.UserModel>=DisplayLevel.Notify ||...
            InfoLevel.UserModelDebug>=DisplayLevel.Final
        if ~isempty(setdiff(this.FlowTable.Identifier(fidx), this.reject_index))
            cprintf('SystemCommands', '[%s]: reject flows due to zero rate.\n', this.options.Method);
        end
    end
end
this.reject_index = this.FlowTable.Identifier(fidx);
end
