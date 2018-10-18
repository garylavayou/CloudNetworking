%% Profit of the slice
% Evaluate the function value and gradient of the objective.
%
% When prices are specified (explicitly or implicitly), the profit is the total utility
% of the slice minus the total payment to the slice provider;

%% Input Arguments
% * |slice|: providing slice information.
%       *LinkPrice*:
%       *NodePrice*:
% * |options|:
%       *PricingPolicy*:
%       *bFinal*:
function profit = getProfit(slice, options)
% determine varriables.
if nargin >= 2 && isfield(options, 'bFinal') && options.bFinal
    vars = [slice.Variables.x; slice.Variables.z];
else
    vars = [slice.temp_vars.x; slice.temp_vars.z];
end
if nargin >= 2
    if isfield(options, 'LinkPrice')
        slice.prices.Link = options.LinkPrice;
    elseif isfield(options, 'bFinal') && options.bFinal 
        slice.prices.Link = slice.Links.Price;
    end
    if isfield(options, 'NodePrice')
        slice.prices.Node = options.NodePrice;
    elseif isfield(options, 'bFinal') && options.bFinal
        slice.prices.Node= slice.ServiceNodes.Price;
    end
end
if nargin < 2 || ~isfield(options, 'PricingPolicy')
    options.PricingPolicy = 'linear';
    warning('%s: <PricingPolicy> not specifed, set to [%s].', ...
        calledby, options.PricingPolicy);
end
%%
% Here, the invoked method must be <Slice.fcnProfit>.
% Use class name to avoid dynamic loading of class.
% Subclasses may override this method to define different ways to calculate profits.
options.bFinal = true;
profit = SimpleSlice.fcnProfit(vars, slice, options); 
slice.prices.Node = [];
slice.prices.Link = [];
end