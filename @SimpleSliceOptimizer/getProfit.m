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
function profit = getProfit(this, options)
defaultopts = getstructfields(this.hs.options, {'PricingPolicy'});
options = structmerge(defaultopts, options);

% determine varriables.
if nargin >= 2 && isfield(options, 'bFinal') && options.bFinal
    vars = [this.Variables.x; this.Variables.z];
else
    vars = [this.temp_vars.x; this.temp_vars.z];
end
if nargin >= 2
    if isfield(options, 'LinkPrice')
        this.prices.Link = options.LinkPrice;
    elseif isfield(options, 'bFinal') && options.bFinal 
        this.prices.Link = this.hs.Links.Price;
    end
    if isfield(options, 'NodePrice')
        this.prices.Node = options.NodePrice;
    elseif isfield(options, 'bFinal') && options.bFinal
        this.prices.Node= this.hs.ServiceNodes.Price;
    end
end
%%
% Here, the invoked method must be <Slice.fcnProfit>.
% Use class name to avoid dynamic loading of class.
% Subclasses may override this method to define different ways to calculate profits.
options.bFinal = true;
profit = SimpleSliceOptimizer.fcnProfit(vars, this, options); 
this.prices.Node = [];
this.prices.Link = [];
end