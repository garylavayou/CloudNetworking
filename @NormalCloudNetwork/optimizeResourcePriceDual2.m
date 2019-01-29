%% Optimize Resource Price
% Firstly, search a price 'αρ+λ' that close to the optimal price.
% "λ" is detemined by using "Dual-ADMM" or "Partial Inverse" method. The procedure is
% relatively slow.
% * Resource Cost Model: linear
% DATE: 2019-01-19
function varargout = optimizeResourcePriceDual2(this, varargin)
[slices, capacities, unitcosts, fields, options] = ...
	this.preOptimizeResourcePrice(varargin{:});
defaultopts = Dictionary('PricingPolicy', 'linear', 'OptimizeOrder', 0, 'unit', 1); 
options = structmerge(defaultopts, options, 'exclude');
options.Capacity = capacities;
if options.bCountTime
	t_start = tic; prt = 0; srt = 0;
end
if nargin <= 2 || ~isfield(options, 'TuningMethod')
	iter_method = 'DualADMM';		% {DualDecomposition, DualADMM, PartialInverse}
else
	iter_method = options.TuningMethod;
	options = rmfield(options, 'TuningMethod');
end
Ns = length(slices);
for i = 1:Ns
  slices(i).initialize();
	slices(i).Optimizer.initializeParallel('priceOptimalFlowRateDA', options);
end
if options.bCountTime
	t_stop = toc(t_start); prt = prt + t_stop; srt = srt + t_stop; t_start = tic;
end

%% Trail Prices
% Find the base price 'μ+αρ' that possibly maximize the profit of SP, while maximizing the
% social welfare.
% The slice provider continually increases the trail price based on the resource cost. The
% joint slice dimensioning problem then find out a social welfare price for slices.
sp_profit = -inf*ones(3,1); b_violate = true(3,1);
num_iters = 0;
delta_price = mean(unitcosts);
t = 0;
while true
	trial_price = unitcosts + delta_price * 2^t;
	b_violate(1:2) = b_violate(2:3);
	sp_profit(1:2) = sp_profit(2:3);
	switch iter_method
		case 'DualADMM'
			[sp_profit(3), b_violate(3), output] = ...
				this.SolveSCPDA(slices, trial_price, options);			% % return results
		case 'PartialInverse'
			[sp_profit(3), b_violate(3), output] = ...
				this.SolveSCPPP(slices, trial_price, options); %#ok<ASGLU>
		otherwise
			error('error: un-recognized method.');
	end
	RECORD_SUBROUTINE_TIME;
	num_iters = num_iters + 1;
	if num_iters == 1
		options.bInitialize = false;
	end
	if sp_profit(3) <= sp_profit(2)
		break;
	end
	t = t+1;
end
%% Joint Slice Dimensioning Problem: Second Stage
% Narrow down the price gap to find the price that could achieve
% higher profit for SP.
beta_L = 2^(t-2);
beta_M = 2^(t-1);
beta_R = 2^t;
epsilon = 10^-3;
beta_N = beta_L;
while abs(beta_N-beta_M)/beta_M > epsilon
	if sp_profit(1)<=sp_profit(3)
		beta_N = (beta_L+beta_M)/2;
	else
		beta_N = (beta_M+beta_R)/2;
	end
	trial_price = unitcosts + delta_price * beta_N;
	switch iter_method
		case 'DualADMM'
			[sp_profit_N, ~, output] = this.SolveSCPDA(slices, trial_price, options);			% return results
		case 'PartialInverse'
			[sp_profit_N, ~, output] = this.SolveSCPPP(slices, trial_price, options); 
		otherwise
			error('error: un-recognized method.');
	end
	RECORD_SUBROUTINE_TIME;
	if sp_profit(1)<=sp_profit(3) 
		if sp_profit_N <= sp_profit(2) 
			sp_profit(1) = sp_profit_N;
			beta_L = beta_N;
		else
			sp_profit(3) = sp_profit(2);
			sp_profit(2) = sp_profit_N;
			beta_R = beta_M;
			beta_M = beta_N;
			beta_N = beta_R;
		end
	else
		if sp_profit_N <= sp_profit(2) 
			sp_profit(3) = sp_profit_N;
			beta_R = beta_N;
		else
			sp_profit(1) = sp_profit(2);
			sp_profit(2) = sp_profit_N;
			beta_L = beta_M;
			beta_M = beta_N;
			beta_N = beta_L;
		end
	end
end
%% Finalize substrate network
if options.bCountTime
	t_start = tic;
end
this.finalize(this.convertParameter(output.Prices), slices);
if options.bCountTime
	t_stop = toc(t_start); prt = prt + t_stop/Ns; srt = srt + t_stop;
	this.op.runtime = struct('Parallel', prt, 'Serial', srt);
	this.op.iterations = num_iters;
end
[varargout{1:nargout}] = this.calculateOutput(slices, fields, ...
	struct('PricingPolicy', 'linear'));

end
