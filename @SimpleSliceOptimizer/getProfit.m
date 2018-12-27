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
if nargin <= 1
	options = Dictionary;
else
	options = Dictionary(options);
end
setdefault(options, this.hs.options, {'PricingPolicy'});

% determine varriables.
if this.hs.isFinal()
    vars = [this.Variables.x; this.Variables.z];
else
    vars = [this.temp_vars.x; this.temp_vars.z];
end
if isfield(options, 'LinkPrice')
	this.prices.Link = options.LinkPrice;
elseif this.hs.isFinal()
	this.prices.Link = this.hs.Links.Price;
end
if isfield(options, 'NodePrice')
	this.prices.Node = options.NodePrice;
elseif this.hs.isFinal()
	this.prices.Node= this.hs.ServiceNodes.Price;
end
%%
% Here, the invoked method must be <Slice.fcnProfit>.
% Use class name to avoid dynamic loading of class.
% Subclasses may override this method to define different ways to calculate profits.
options.bFinal = true;
options.bCompact = false;
profit = SimpleSliceOptimizer.fcnProfit(vars, this, options); 
this.setProblem('Price', []);
end