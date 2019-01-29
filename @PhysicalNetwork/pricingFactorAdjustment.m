%% Adjust the pricing factor
% Capacity is predetermined for each slice.
% If 'PricingFactor' is specified, use the given 'PricingFactor'. Otherwise
% (options.PricingFactor=0),  adjust the pricing factor by three-division method.
function [prices, results] = pricingFactorAdjustment(this, options)
defaultopts = getstructfields(this.options, {'PricingPolicy'}, 'error');
if nargin <= 1
	options = defaultopts;
else
	options = structmerge(defaultopts, options);
end
options.bInitialize = true;   % isFinalize = false;
if options.bCountTime
  prt = 0; srt = 0;
end

unitcost = options.UnitCost;
if this.options.PricingFactor ~= 0
    PricingFactor_h = this.options.PricingFactor;
else
    PricingFactor_h = 0.5;
    sp_profit = -inf;
end
PricingFactor_l = 0;
k = 0;
while true
    prices = unitcost * (1 + PricingFactor_h);
		[sp_profit_new, output] = priceIteration(this, prices, options); %#ok<ASGLU>
		RECORD_SUBROUTINE_TIME;
		k = k + 1;
		if k == 1
			options.bInitialize = false;
		end
    if this.options.PricingFactor ~= 0
        break;
    end
    %%%
    % compute and compare the SP's profit
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
					prices = unitcost * (1 + m(i));
					[sp_profit(i), output] = priceIteration(this, prices, options); %#ok<ASGLU>
					RECORD_SUBROUTINE_TIME;
				end
        if sp_profit(1) > sp_profit(2)
            PricingFactor_h = m(2);
        else
            PricingFactor_l = m(1);
        end
    end
end

if options.bCountTime
	results.runtime = struct('Parallel', prt, 'Serial', srt);
elseif narout>= 2
	results = struct;
end
end

%% ISSUE: VNFlist is not conmmonly shared.
function [sp_profit, output] = priceIteration(this, prices, options)
output = struct;
if options.bCountTime
	output.runtime = struct('Parallel', 0, 'Serial', 0);
end
prices = this.convertParameter(prices);

warning('TODO: parallel, see <SolveSCPCC>.');
for s = 1:this.NumberSlices
	sl = this.slices(s);
	sl.Optimizer.setProblem('Price', prices);
	%%%
	% optimize each slice with price and resource constraints.
	if options.bCountTime
		t1 = tic;
	end
	sl.Optimizer.optimalFlowRate(options);
	if options.bCountTime
		t2 = toc(t1);
		output.runtime.Serial = output.runtime.Serial + toc(t2);
		output.runtime.Parallel = max(output.runtime.Parallel, toc(t2));
	end
	sl.Optimizer.setProblem('Price', [])
end
sp_profit = this.getSliceProviderProfit([], prices, options);
end
