%% Calculate the output
% The output variables includes,
%
% # Price and Load of physical nodes and links;
% # The net social welfare;
% # The profit of each slice and the substrate network .
% # Flow rate of all flows in the network.
%
% To calculate net social welfare, two methods might be used, i.e. |Approximate|,
% |Accurate|.
%
% the profit type with |Percent| has been deprecated.
%%
% |argin|: existing field for output.
function argout = calculateOutput(this, argin, options)
if ~isempty(argin)
    argout = argin;
end
if isfield(options, 'ProfitType')
    profit_type = options.ProfitType;
else
    warning('profit type is set as AccuratePrice.');
    profit_type = {'AccuratePrice'};
end
if isfield(options, 'WelfareType')
    welfare_type = options.WelfareType;
else
    warning('welfare type is set as Accurate.');
    welfare_type = {'Accurate'};
end
% Remove the unnecessary field in |options|.
if isfield(options, 'NodePrice')
    options = rmfield(options, {'NodePrice', 'LinkPrice'});
end
    
argout.link_price = this.getLinkField('Price');
% dc_node_index = find(this.Topology.Nodes.DataCenter);
argout.node_price = this.getNodeField('Price');            
argout.link_load = this.getLinkField('Load');
argout.node_load = this.getNodeField('Load');       
argout.flow_rate = [];

if cellstrfind(welfare_type, 'Accurate')
    argout.welfare_accurate = 0;
end
if cellstrfind(welfare_type, 'Approximate')
    argout.welfare_approx = 0;
end
n_profit_type = length(profit_type);
profit_table = zeros(this.NumberSlices+1, n_profit_type);

for s = 1:this.NumberSlices
    sl = this.slices{s};
    argout.flow_rate = [argout.flow_rate; sl.FlowTable.Rate];
    var_x = [sl.Variables.x; sl.Variables.z];
    
    %% Calculate the net social welfare
    % Here, we use two model to assess the net social welfare, 
    %
    % # *Approximate Model*: the net social welfare is the sum of net social welfare from
    %    all slices. For slices, the net social welfare is $slice\_welfare = (utility -
    %    approximate\_cost)$, where the approximate cost is computed by our approximation
    %    formula.
    % # *Accurate Model*: the net social welfare is the total utility less the total network
    %   cost.
    %
    % Without considering the static cost, the two methods are the same.
    if ~isempty(cellstrfind(welfare_type, 'Accurate'))
        argout.welfare_accurate = argout.welfare_accurate...
            + sl.weight*sum(fcnUtility(sl.FlowTable.Rate));
        %                    + sl.constant_profit;  => move to CloudNetworkEx
    end
    if ~isempty(cellstrfind(welfare_type, 'Approximate'))
        sw = Slice.fcnSocialWelfare(var_x, sl, 'Approximate');
        argout.welfare_approx = argout.welfare_approx + sw;
    end
    
    %% Calculate the profit of slices
    % To calculate profit, four profit methods might be used, i.e. |ApproximatePercent|,
    % |ApproximatePrice|, |AccuratePercent|, |AccuratePrice|.
    %
    % Note: using the |AccuratePercent| method is always not accurate, since it cannnot
    % precisely compute the static cost of a slice (the static cost should be shared
    % between slices). 
    %
    % Without considering the static cost, the |Approximate*| method and |Accurate*|
    % methods are the same. 
    for i = 1:length(profit_type)
        switch profit_type{i}
            case 'ApproximatePercent'   % <deprecated>
                profit_table(s,i) = options.PercentFactor * sw;
            case 'AccuratePercent'      % <deprecated>
                sw = Slice.fcnSocialWelfare(var_x, sl, 'Accurate');
%                 if this.static_factor ~= 0   => move to CloudNetworkEx
%                     nid = sl.DCNI;
%                     idx = sl.VirtualDataCenters.Load>0;
%                     sw = sw - ...
%                         dot(sl.VirtualDataCenters.Load(idx)./this.getNodeField('Load', nid(idx)),...
%                         this.getNodeField('StaticCost', nid(idx)));
%                 end
                profit_table(s,i) = options.PercentFactor * sw;
            case {'ApproximatePrice','AccuratePrice'}
                profit_table(s,i) = Slice.fcnProfit(var_x, sl, options);
            otherwise 
                error('error: invalid profit type.');
        end
    end
    %%%
    % * *Proportional profit with approximate cost*  <deprecated>
    %
    % For slices, the approximate profit is given by 
    %
    %      slice_profit = slice_welfare * pf,
    %
    %  The proportional profit means only part of the profit is assigned to the slice,
    %  another part should be deliver to the slice provider. As a result, for slice
    %  provider, 
    %
    %      net_profit = sum(slice_welfare) * (1-pf).
    %
    % The profit of slice provider is the sum of net social welfarefrom each slice
    % computed by the _Approximate Model_. It is the same computing method, when using
    % _Accurate Model_ with proportion profit.
    %
    % *NOTE*: fcnNetProfit here evaluate the raw profit with resource cost.
    %
    % * *Proportional profit with accurate cost* <deprecated>
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
    %
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
end
% if ~isempty(this.eta)
%     embed_profit_approx = this.eta*this.getNetworkCost;
%     embed_profit_accurate = this.getNetworkCost([], [], options.Model);
% else
%     embed_profit_approx = 0;
%     embed_profit_accurate = 0;
% end
% argout.welfare_approx = argout.welfare_approx + embed_profit_approx;
if cellstrfind(welfare_type, 'Accurate')
    argout.welfare_accurate = argout.welfare_accurate - ...
        this.getNetworkCost([],[], 'Accurate');
end

if cellstrfind(profit_type, 'Percent')
    f = (1-options.PercentFactor)/options.PercentFactor;
end
%%%
% calculate the profit of substrate network
for i = 1:length(profit_type)
    switch profit_type{i}
        case {'ApproximatePercent', 'AccuratePercent'}
            profit_table(end,i) = sum(profit_table(1:(end-1),i))*f;
        case 'ApproximatePrice'
            profit_table(end,i) = argout.welfare_approx - sum(profit_table(1:(end-1),i));
        case 'AccuratePrice'
            profit_table(end,i) = argout.welfare_accurate - sum(profit_table(1:(end-1),i));
        otherwise
            error('error: invalid profit type.');
    end
end

argout.profit = array2table(profit_table, 'VariableNames', profit_type);
end