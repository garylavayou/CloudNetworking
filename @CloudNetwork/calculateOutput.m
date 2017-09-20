%% Calculate the output
% The output variables includes,
%
% # Price and Load of physical nodes and links;
% # The net social welfare;
% # The profit of each slice and the substrate network .
% # Flow rate of all flows in the network.
%%
% |argin|: existing field for output.
%% TODO: move ProfitType and WelfareType to <CloudNetworkEx>
function argout = calculateOutput(this, argin, options)
if ~isempty(argin)
    argout = argin;
end
% Remove the unnecessary field in |options|.
if isfield(options, 'NodePrice')
    options = rmfield(options, {'NodePrice', 'LinkPrice'});
end
    
argout.link_price = this.getLinkField('Price');
argout.node_price = this.getNodeField('Price');            
argout.link_load = this.getLinkField('Load');
argout.node_load = this.getNodeField('Load');       
argout.flow_rate = [];
argout.welfare = 0;

profit_table = zeros(this.NumberSlices+1, 1);

for s = 1:this.NumberSlices
    sl = this.slices{s};
    argout.flow_rate = [argout.flow_rate; sl.FlowTable.Rate];
    var_x = [sl.Variables.x; sl.Variables.z];
    
    %% Calculate the net social welfare
    %       The net social welfare is the total utility less the total network cost.
    argout.welfare = argout.welfare + sl.weight*sum(fcnUtility(sl.FlowTable.Rate));
    
    % Calculate the profit of slices
    profit_table(s) = Slice.fcnProfit(var_x, sl, options);
    %%%
    % * *Net profit with offered price*
    %
    % For slices,
    %
    %      net_profit = utility - payment(price).
    %
    % The price is calculated accoding to optimization procedure.
end
argout.welfare = argout.welfare - this.getNetworkCost;
% calculate the profit of substrate network
profit_table(end) = argout.welfare - sum(profit_table(1:(end-1)));
argout.profit = profit_table;
end