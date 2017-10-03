%% Calculate the output
% The output variables includes,
%
% # Price and Load of physical nodes and links;
% # The net social welfare;
% # The profit of each slice and the substrate network .
% # Flow rate of all flows in the network.
%%
% |argin|: existing field for output.
function argout = calculateOutput(this, argin, new_opts)
if ~isempty(argin)
    argout = argin;
end
    
argout.LinkPrice = this.getLinkField('Price');
argout.NodePrice = this.getNodeField('Price');            
argout.LinkLoad = this.getLinkField('Load');
argout.NodeLoad = this.getNodeField('Load');       
argout.FlowRate = [];
argout.Welfare = 0;

profit_table = zeros(this.NumberSlices+1, 1);

for s = 1:this.NumberSlices
    sl = this.slices{s};
    argout.FlowRate = [argout.FlowRate; sl.FlowTable.Rate];
    var_x = [sl.Variables.x; sl.Variables.z];
    
    %% Calculate the net social welfare
    %       The net social welfare is the total utility less the total network cost.
    argout.Welfare = argout.Welfare + sl.weight*sum(fcnUtility(sl.FlowTable.Rate));
    
    % Calculate the profit of slices
    profit_table(s) = Slice.fcnProfit(var_x, sl, ...
        getstructfields(new_opts, {'PricingPolicy'}));
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