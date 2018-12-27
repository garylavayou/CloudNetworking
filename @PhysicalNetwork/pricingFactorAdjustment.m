%% Adjust the pricing factor
% If 'PricingFactor' is specified, use the given 'PricingFactor'. Otherwise
% (options.PricingFactor=0),  adjust the pricing factor by three-division method.
function [prices, runtime] = pricingFactorAdjustment(this, options)
defaultopts = getstructfields(this.options, {'PricingPolicy'}, 'error');
if nargin <= 1
	options = defaultopts;
else
	options = structmerge(defaultopts, options);
end
options.bInitialize = true;   % isFinalize = false;

if nargout == 3
    runtime.Serial = 0;
    runtime.Parallel = 0;
end
link_uc = this.getLinkCost;
node_uc = this.getNodeCost;
if this.op.options.PricingFactor ~= 0
    PricingFactor_h = options.PricingFactor;
else
    PricingFactor_h = 0.5;
    sp_profit = -inf;
end
PricingFactor_l = 0;
k = 0;
while true
    prices.Link = link_uc * (1 + PricingFactor_h);
    prices.Node = node_uc * (1 + PricingFactor_h);
    if nargout == 3
        t = priceIteration(this, prices, options);
        runtime.Serial = runtime.Serial + t.Serial;
        runtime.Parallel = runtime.Parallel + t.Parallel;
    else
        priceIteration(this, prices, options);
		end
		k = k + 1;
		if k == 1
			options.bInitialize = false;
		end
    if this.options.PricingFactor ~= 0
        break;
    end
    %%%
    % compute and compare the SP's profit
    sp_profit_new = this.getSliceProviderProfit([], prices, options);
    if sp_profit_new > sp_profit
        PricingFactor_l = PricingFactor_h;
        PricingFactor_h = PricingFactor_h * 2;
        sp_profit = sp_profit_new;
    else
        break;
    end
end
if this.options.PricingFactor == 0
    while PricingFactor_h-PricingFactor_l>=0.1
        sp_profit = zeros(2,1);
        m = zeros(2,1);
        for i = 1:2
            m(i) = (i*PricingFactor_h+(3-i)*PricingFactor_l)/3;
            prices.Link = link_uc * (1 + m(i));
            prices.Node = node_uc * (1 + m(i));
            if nargout == 3
                t = priceIteration(this, prices, options);
                runtime.Serial = runtime.Serial + t.Serial;
                runtime.Parallel = runtime.Parallel + t.Parallel;
            else
                priceIteration(this, prices, options);
            end
            sp_profit(i) = this.getSliceProviderProfit([], prices, options);
        end
        
        if sp_profit(1) > sp_profit(2)
            PricingFactor_h = m(2);
        else
            PricingFactor_l = m(1);
        end
    end
end
end
