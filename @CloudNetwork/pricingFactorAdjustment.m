%% Adjust the pricing factor
% three-division method
function [node_price, link_price, runtime] = pricingFactorAdjustment(this, options)
if nargout == 3
    runtime.Serial = 0;
    runtime.Parallel = 0;
end
link_uc = this.getLinkCost;
node_uc = this.getNodeCost;
if isfield(options, 'PricingFactor')
    PricingFactor.h = options.PricingFactor;
else
    PricingFactor.h = 0.5;
    sp_profit = -inf;
end
PricingFactor.l = 0;
while true
    link_price = link_uc * (1 + PricingFactor.h);
    node_price = node_uc * (1 + PricingFactor.h);
    if nargout == 3
        t = priceIteration(this, node_price, link_price, options);
        runtime.Serial = runtime.Serial + t.Serial;
        runtime.Parallel = runtime.Parallel + t.Parallel;
    else
        priceIteration(this, node_price, link_price, options);
    end
    if isfield(options, 'PricingFactor')
        break;
    end
    %%%
    % compute and compare the SP's profit
    sp_profit_new = this.getSliceProviderProfit(node_price, link_price);
    if sp_profit_new > sp_profit
        PricingFactor.l = PricingFactor.h;
        PricingFactor.h = PricingFactor.h * 2;
        sp_profit = sp_profit_new;
    else
        break;
    end
end
if ~isfield(options, 'PricingFactor')
    while PricingFactor.h-PricingFactor.l>=0.1
        sp_profit = zeros(2,1);
        m = zeros(2,1);
        for i = 1:2
            m(i) = (i*PricingFactor.h+(3-i)*PricingFactor.l)/3;
            link_price = link_uc * (1 + m(i));
            node_price = node_uc * (1 + m(i));
            if nargout == 3
                t = priceIteration(this, node_price, link_price, options);
                runtime.Serial = runtime.Serial + t.Serial;
                runtime.Parallel = runtime.Parallel + t.Parallel;
            else
                priceIteration(this, node_price, link_price, options);
            end
            sp_profit(i) = this.getSliceProviderProfit(node_price, link_price);
        end
        
        if sp_profit(1) > sp_profit(2)
            PricingFactor.h = m(2);
        else
            PricingFactor.l = m(1);
        end
    end
end
end