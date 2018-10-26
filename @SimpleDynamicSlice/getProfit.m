%% get profit
% This function is called during slice dimensioning.
% Profit of Slice Customer: utility - resource consumption payment - reconfiguration cost;
%
% The calculation methods are different during/after optimization:
%   (1) In the optimization, we use |temp_vars| to calculate resource consumption and
%       reconfiguration cost. During optimization, we use L1-norm to approximate the
%       reconfiguration of slices. Therefore, when we calculate the profit, we still use
%       the approximated reconfiguration cost. Although there is error from the true
%       reconfiguration cost, it is relatively small compared with the resource
%       consumption cost.
%   (2) After the optimization is completed, we use the final results (|Variables|) to
%       calculate the resource consumption and reconfiguration cost. Correspondingly, the
%       reconfiguration cost is calculated with the original formulation (according to the
%       difference of states).
%
% NOTE: If there is no resource reservation or reconfiguration cost constraint, we call
% the superclass method <getProfit> to calculate the profit. Due to constraint of
% reconfiguration cost or resource reservation, the link/VNF capacity might be larger than
% the demand.
%
% See also <Slice.getProfit>.
function profit = getProfit(this, options)
if nargin == 1
    options = struct;
end
if this.invoke_method == 0
    profit = getProfit@SimpleSlice(this, options);
    return;
end

options = structmerge(rmstructfields(options, 'bCompact'),...
    getstructfields(options, 'PricingPolicy', 'default', {'quadratic'}));
bFinal = this.isFinal();

options.bFinal = true;      % 'bFinal' is set to output the real profit (min -f => max f).
options.num_varx = length(this.temp_vars.x);
options.num_varz = length(this.temp_vars.z);
options.num_varv = length(this.temp_vars.v);
if bFinal
    this.prices.Link = this.Links.Price;
    this.prices.Node = this.ServiceNodes.Price;
    vars = [this.Variables.x; this.Variables.z; this.Variables.v;...
        this.Links.Capacity];
    profit = DynamicSlice.fcnProfitReserveSlicing(vars, this, options);
    if this.invoke_method >= 2 && this.b_dim   
        % only slices going through redimensioning have reconfiguration cost.
        profit = profit - this.get_reconfig_cost('const');
    end
else
    if nargin >= 2
        if isfield(options, 'LinkPrice')
            this.prices.Link = options.LinkPrice;
        end
        if isfield(options, 'NodePrice')
            this.prices.Node = options.NodePrice;
        end
    end
    if this.invoke_method == 1
        vars = [this.temp_vars.x; this.temp_vars.z; this.temp_vars.v;...
            this.temp_vars.c];
        profit = DynamicSlice.fcnProfitReserveSlicing(vars, this, options);
    else
        vars = [this.temp_vars.x; this.temp_vars.z; this.temp_vars.v;...
            this.temp_vars.tx; this.temp_vars.tz; this.temp_vars.tv;...
            this.temp_vars.c];
        profit = DynamicSlice.fcnProfitReconfigureSlicing(vars, this, options);
    end
end

this.prices.Node = [];
this.prices.Link = [];
end
