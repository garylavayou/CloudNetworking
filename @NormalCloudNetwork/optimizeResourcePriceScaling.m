function varargout = optimizeResourcePriceScaling(this, varargin)
global DEBUG; %#ok<NUSED>

%% Initialization
[slices, ~, unitcosts, fields, options] = this.preOptimizeResourcePrice(varargin{:});
defaultopts = Dictionary('PricingPolicy', 'linear', 'unit', 1);
options = structmerge(options, defaultopts, 'exclude');
if options.bCountTime
	t_start = tic; prt = 0; srt = 0;
end
Ns = length(slices);
for i = 1:Ns
	slices(i).initialize();
	slices(i).Optimizer.initializeParallel('priceOptimalFlowRate', options);
end
if options.bCountTime
	t_stop = toc(t_start); prt = prt + t_stop; srt = srt + t_stop; t_start = tic;
end

%% Trial Prices
% Find the price 'αp' that possibly maximize the profit of SP.
sp_profit = -inf*ones(3,1); b_violate = true(3,1);
num_iters = 0;
t = 0;
while true
	trial_price = unitcosts * 2^t;
	b_violate(1:2) = b_violate(2:3);
	sp_profit(1:2) = sp_profit(2:3);
	[sp_profit(3), b_violate(3), output] = ...
		this.SolveSCPCC(slices, trial_price, options); %#ok<ASGLU>
	RECORD_SUBROUTINE_TIME;
	num_iters = num_iters + 1;
	if num_iters == 1
		options.bInitialize = false;
	end
	if sp_profit(3) <= sp_profit(2) && ~b_violate(3)
		break;
	end
	t = t+1;
end
% statisfy condition: sp_profit(3) <= sp_profit(2) && ~b_violate(3)
beta_R = 2^t;
epsilon = 10^-3;
sp_profit_N = zeros(2,1);
if ~b_violate(2)
	% =>    sp_profit(2) >= sp_profit(1) && sp_profit(2) >= sp_profit(3).
	%
	% Both the last price and current price are feasible in terms of capacity.
	% serach with 'pure-tri-section' method.
	beta_L = 2^(t-2);
	while abs((sp_profit(2)-sp_profit(1))/sp_profit(1)) >= epsilon || ...
			abs((sp_profit(2)-sp_profit(3))/sp_profit(3)) >= epsilon
		% abs(beta_L-beta_R)/beta_L > epsilon
		beta_N = [beta_L*2+beta_R; beta_L+beta_R*2]/3;
		for i = 2:-1:1
			trial_price = unitcosts* beta_N(i);
			[sp_profit_N(i), ~, output] = ...
				this.SolveSCPCC(slices, trial_price, options); %#ok<ASGLU>
			RECORD_SUBROUTINE_TIME;
		end
		num_iters = num_iters + 2;
		if sp_profit_N(1) > sp_profit_N(2)
			beta_R = beta_N(2);	sp_profit(3) = sp_profit_N(2); sp_profit(2) = sp_profit_N(1);
		else
			beta_L = beta_N(1);	sp_profit(1) = sp_profit_N(1); sp_profit(2) = sp_profit_N(2);
		end
	end
	results.sp_profit = sp_profit_N(1);
	results.beta = beta_N(1);
	results.procedure = 'max';
else
	% =>	sp_profit(1) >= sp_profit(2) >= sp_profit(3)   or 
	% =>	sp_profit(2) >= sp_profit(1) && sp_profit(2) >= sp_profit(3).
	%
	% The last price is not feasible, while current price is feasible in terms of capacity.
	% serach with 'trim-tri-section' method.
	beta_L = 2^(t-1);
	while abs((sp_profit(2)-sp_profit(1))/sp_profit(1)) >= epsilon || ...
			abs((sp_profit(2)-sp_profit(3))/sp_profit(3)) >= epsilon
		% abs(beta_R-beta_L)/beta_L > epsilon
		beta_N = [beta_L*2+beta_R; beta_L+beta_R*2]/3;
		for i = 1:2
			trial_price = unitcosts * beta_N(i);
			[sp_profit_N(i), b_violate(i), output] = ...
				this.SolveSCPCC(slices, trial_price, options); %#ok<ASGLU>
			RECORD_SUBROUTINE_TIME;
		end
		num_iters = num_iters + 2;
		if b_violate(2)
			beta_L = beta_N(2); sp_profit(1) = sp_profit_N(2); sp_profit(2) = sp_profit_N(2);
		elseif b_violate(1)
			beta_L = beta_N(1); sp_profit(1) = sp_profit_N(1); sp_profit(2) = sp_profit_N(2);
			if sp_profit_N(1) > sp_profit_N(2)
				beta_R = beta_N(2); sp_profit(3) = sp_profit_N(2);
			end
		else
			if sp_profit_N(1) > sp_profit_N(2)
				beta_R = beta_N(2); sp_profit(3) = sp_profit_N(2); sp_profit(2) = sp_profit_N(1);
			else
				beta_L = beta_N(1); sp_profit(1) = sp_profit_N(1); sp_profit(2) = sp_profit_N(2);
			end
		end
	end
	results.sp_profit = sp_profit(3);
	results.beta = beta_R;
	results.procedure = 'feas-max';
end

%% Finalize substrate network
if options.bCountTime
	t_start = tic;
end
this.finalize(this.convertParameter(trial_price), slices);
if options.bCountTime
	t_stop = toc(t_start); prt = prt + t_stop/Ns; srt = srt + t_stop;
	this.op.runtime = struct('Parallel', prt, 'Serial', srt);
	this.op.iterations = num_iters;
end
[varargout{1:nargout}] = this.calculateOutput(slices, fields, ...
	struct('PricingPolicy', 'linear'));
if nargout >= 3
	varargout{3}.Beta = results.beta;
	varargout{3}.Procedure = results.procedure;
end

end

