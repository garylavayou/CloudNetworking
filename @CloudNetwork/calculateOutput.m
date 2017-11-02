%% Calculate the output
% The output variables includes,
%
% # Price and Load of physical nodes and links;
% # The net social welfare;
% # The profit of each slice and the substrate network .
% # Flow rate of all flows in the network.
%%
% |argin|: existing field for output.
function argout = calculateOutput(this, argin, options)
if nargin>= 2 && ~isempty(argin)
    argout = argin;
end
if nargin < 3
    options = struct;
end
options = getstructfields(options, {'Slices', 'PricingPolicy'}, 'default-ignore', ...
    {this.slices,[]});
if isempty(options.PricingPolicy)
    options.PricingPolicy = 'linear';
    warning('Slice.calculateOutput: <PricingPolicy> set to default [%s].', ...
        options.PricingPolicy);
end
options.bFinal = true;

argout.LinkPrice = this.getLinkField('Price');
argout.NodePrice = this.getNodeField('Price');            
argout.LinkLoad = this.getLinkField('Load');
argout.NodeLoad = this.getNodeField('Load');       
argout.FlowRate = [];
argout.Welfare = 0;

profit_table = zeros(this.NumberSlices+1, 1);

for s = 1:length(options.Slices)
    sl = options.Slices{s};
    argout.FlowRate = [argout.FlowRate; sl.FlowTable.Rate];
    %% Calculate the net social welfare
    %       The net social welfare is the total utility less the total network cost.
    argout.Welfare = argout.Welfare + sl.weight*sum(fcnUtility(sl.FlowTable.Rate));
    
    % Calculate the profit of slices
    profit_table(s) = sl.getProfit(options);
    %%%
    % * *Net profit with offered price*
    %
    % For slices,
    %
    %      net_profit = utility - payment(price).
    %
    % The price is calculated accoding to optimization procedure.
end
argout.Welfare = argout.Welfare - this.getNetworkCost;
% calculate the profit of substrate network
profit_table(end) = argout.Welfare - sum(profit_table(1:(end-1)));
argout.Profit = profit_table;
end