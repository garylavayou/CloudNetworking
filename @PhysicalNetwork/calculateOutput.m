function argout = calculateOutput(this, argin, options)
argout = argin;

argout.welfare_approx = 0;
argout.welfare_accurate = 0;
t = zeros(this.NumberSlices+1, 1);
argout.profit = table(t, t, t, t, 'VariableNames',...
    {'ApproximatePercent', 'ApproximatePrice', 'AccuratePercent', 'AccuratePrice'});
clear t;
argout.flow_rate = [];
options.Model = 'Accurate';
%% Calculate the net social welfare
% The resource price is also keep fixed since it related to the linear resouce dynamic
% cost, static cost and the prcing factor. Here, we use two model to assess the net
% social welfare,
%
% * *Approximate Model*: the net social welfare is the sum of raw profit of all slices.
% * *Accurate Model*: the net social welfare is the total utility less the total network
% cost.
for s = 1:this.NumberSlices
    sl = this.slices{s};
    argout.flow_rate = [argout.flow_rate; sl.FlowTable.Rate];
    
    var_x = [sl.Variables.x; sl.Variables.z];
    %% Calculate the profit
    % * *Proportional Net profit with approximate cost*
    %
    % For slices,
    %
    %      slice_profit = (utility - approximate_cost) * pf,
    %
    % where the approximate cost is computed by our approximation formula. The
    % proportional net profit means only one part of the raw net profit is attributed to
    % the slice, another part should be deliver to the substrate network. As a result,
    % for substrate network,
    %
    %      net_profit = sum(slice_profit) * (1-pf).
    %
    % The net profit of substrate network is the sum of proportion profit from each slice
    % computed by the _Approximate Model_. It is the same computing method, when using
    % _Accurate Model_ with proportion profit.
    %
    % *NOTE*: fcnNetProfit here evaluate the raw profit with resource cost.
    p = Slice.fcnNetProfit(var_x, sl);
    argout.welfare_approx = argout.welfare_approx + p;
    argout.welfare_accurate = argout.welfare_accurate...
        + sl.weight*sum(log(sl.FlowTable.Rate));
    argout.profit.ApproximatePercent(s) = options.PercentFactor * p;
    %%%
    % * *Proportional Net profit with accurate cost*
    %
    % For slices,
    %
    %      net_profit = (utility - accurate_cost) * pf.
    %
    % The accurate cost is based on the actual resource allocation after the optimization.
    % Slices which share a physical resource will also proportionally share the static
    % cost of the resource.
    %
    % *NOTE*: it is not accurate to calculate the static cost of each slice with this
    % method.
    p = Slice.fcnNetProfit(var_x, sl, options);
    if this.static_factor ~= 0
        nid = sl.VirtualNodes.PhysicalNode;
        idx = sl.VirtualNodes.Load>0;
        p = p - dot(sl.VirtualNodes.Load(idx)./this.getNodeField('Load', nid(idx)),...
            this.getNodeField('StaticCost', nid(idx)));
    end
    argout.profit.AccuratePercent(s) = options.PercentFactor * p;
    %%%
    % * *Net profit with offered price*
    %
    % For slices,
    %
    %      net_profit = utility - payment(price).
    %
    % The price is calculated accoding to optimization procedure.
    %
    % When the resource price is offered, the slice's profit is not related to how the
    % cost is computed. Therefore, |output.profit.ApproximatePrice| and
    % |output.profit.AccuratePercent| of slices is the same. The cost only affects the
    % substrate network's profit. For substrate network with _Approximate Model_ and
    % _Accurate Model_, the net profit is given by accordingly
    %
    %      net_profit = net_social_welfare - sum(slice_profit)
    %                 = sum(slice_payment(price))) - cost(approximate/accurate)
    %
    % *NOTE*: fcnProfit evalute the profit using offered price.
    argout.profit.ApproximatePrice(s) = Slice.fcnProfit(var_x, sl);
end
if ~isempty(this.eta)
    embed_profit_approx = this.eta*this.getNetworkCost;
    embed_profit_accurate = this.getNetworkCost([], [], options.Model);
else
    embed_profit_approx = 0;
    embed_profit_accurate = 0;
end
argout.welfare_approx = argout.welfare_approx + embed_profit_approx;
argout.welfare_accurate = argout.welfare_accurate - ...
    this.getNetworkCost([],[], options.Model) + embed_profit_accurate;

f = (1-options.PercentFactor)/options.PercentFactor;
argout.profit.ApproximatePercent(end) = ...
    sum(argout.profit.ApproximatePercent(1:(end-1))*f) + embed_profit_approx;
argout.profit.AccuratePercent(end) = ...
    sum(argout.profit.AccuratePercent(1:(end-1))*f) + embed_profit_accurate;
argout.profit.ApproximatePrice(end) = argout.welfare_approx ...
    - sum(argout.profit.ApproximatePrice(1:(end-1))) + embed_profit_approx;
argout.profit.AccuratePrice = argout.profit.ApproximatePrice;
argout.profit.AccuratePrice(end) = argout.welfare_accurate ...
    - sum(argout.profit.AccuratePrice(1:(end-1))) + embed_profit_accurate;
end

