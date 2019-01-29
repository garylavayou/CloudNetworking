%% Optimize Resource Price
% Firstly, search a price 'αρ' that close to the optimal price.
% Then use "Dual-ADMM" or "Partial Inverse" method to find the final solution.
% * Resource Cost Model: linear
% DATE: 2018-12-20
%		[stat, slice_stat, results] = optimizeResourcePriceDual(this, slices, options)
function varargout = optimizeResourcePriceDual(this, varargin)
[slices, capacities, unitcosts, fields, options] = ...
	this.preOptimizeResourcePrice(varargin{:});
defaultopts = Dictionary('PricingPolicy', 'linear', 'OptimizeOrder', 0, 'unit', 1); 
options = structmerge(defaultopts, options, 'exclude');
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
% Find the price 'αρ' that possibly maximize the profit of SP.
sp_profit = -inf*ones(3,1); b_violate = true(3,1);
num_iters = 0;
t = 0;
while true
	trial_price = unitcosts * 2^t;
	b_violate(1:2) = b_violate(2:3);
	sp_profit(1:2) = sp_profit(2:3);
	[sp_profit(3), b_violate(3), output] = ...
		this.SolveSCPCC(slices, trial_price, options);  %#ok<ASGLU> % with price model linear.
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

beta_L = 2^(t-2); beta_R = 2^t;
epsilon = 10^-3;
sp_profit_N = zeros(2,1); b_violate_N = zeros(2,1);
while abs((sp_profit(2)-sp_profit(1))/sp_profit(1)) >= epsilon || ...
		abs((sp_profit(2)-sp_profit(3))/sp_profit(3)) >= epsilon
	% abs(beta_L-beta_R)/beta_L > epsilon
	beta_N = [beta_L*2+beta_R; beta_L+beta_R*2]/3;
	for i = 1:2
		trial_price = unitcosts * beta_N(i);
		[sp_profit_N(i), b_violate_N(i), output]= this.SolveSCPCC(slices, trial_price, options); %#ok<ASGLU>
		RECORD_SUBROUTINE_TIME;
	end
	num_iters = num_iters + 2;
	if sp_profit_N(1) > sp_profit_N(2)
		beta_R = beta_N(2);	sp_profit(3) = sp_profit_N(2); 
		sp_profit(2) = max(sp_profit(2),sp_profit_N(1));
	else
		beta_L = beta_N(1);	sp_profit(1) = sp_profit_N(1); 
		sp_profit(2) = max(sp_profit(2),sp_profit_N(2));
	end
end
if ~b_violate_N(1)
	results.beta = beta_N(1);
	results.Prices = unitcosts * results.beta;
	results.Procedure = 'max';	
else
	%%	Serach the final price
	results.beta = (beta_L+beta_R)/2;
	trial_price = unitcosts * results.beta;
	options.Capacity = capacities;
	options.bInitialize = true;
	switch iter_method
		case 'DualADMM'
			options.unit = 1;
			[~, b_violate_N(1),	output] = this.SolveSCPDA(slices, trial_price, options);			% return results
		case 'PartialInverse'
			[~, b_violate_N(1),	output] = this.SolveSCPPP(slices, trial_price, options);
		case 'DualDecomposition'
			error('error: not supported!');
			output = this.SolveSCPDD(slices, trial_price, options); %#ok<UNRCH>
		otherwise
			error('error: un-recognized method.');
	end
	RECORD_SUBROUTINE_TIME;
	num_iters = num_iters + output.numiters;
	options.CapacityConstrained = true;
  options.Capacities = output.Capacities;
	options.bInitialize = true;
	results.Prices = output.Prices;
	results.Procedure = 'dual-iter';
	% [~, b, output] = this.SolveSCPCC(slices, results.Prices, options); %#ok<ASGLU>
	% RECORD_SUBROUTINE_TIME;
	% num_iters = num_iters + 1;
	% 	if b
	% 		warning('%s: Capacity violation.', calledby);
	% 	end
end

%% Finalize substrate network
if options.bCountTime
	t_start = tic;
end
this.finalize(this.convertParameter(results.Prices), slices);
if b_violate_N(1)
	loads = this.convertParameter(this.getNetworkLoad(slices));
	ratio = max(loads./capacities);
	if ratio > 1
		for s = 1:Ns
			sl = this.slices(s);
			op = sl.Optimizer;
			var_fields = fieldnames(op.Variables);
			for j = 1:length(var_fields)
				op.Variables.(var_fields{j}) = op.Variables.(var_fields{j})/ratio;
			end
			sl.FlowTable.Rate = sl.FlowTable.Rate/ratio;
			for fid=1:sl.NumberFlows
				path_list = sl.FlowTable{fid, 'Paths'};
				for j = 1:path_list.Width
					path_list{j}.MultiplyBandwidth(1/ratio);
				end
			end
			sl.ServiceNodes.Load = sl.ServiceNodes.Load/ratio;
			sl.Links.Load = sl.Links.Load/ratio;
			sl.ServiceNodes.Capacity = sl.ServiceNodes.Load;
			sl.Links.Capacity = sl.Links.Load;
		end
	end
end
if options.bCountTime
	t_stop = toc(t_start); prt = prt + t_stop/Ns; srt = srt + t_stop;
	this.op.runtime = struct('Parallel', prt, 'Serial', srt);
	this.op.iterations = num_iters;
end
[varargout{1:nargout}] = this.calculateOutput(slices, fields, ...
	struct('PricingPolicy', 'linear'));
if nargout >= 3
	varargout{3}.Beta = results.beta;
	varargout{3}.Procedure = results.Procedure;
end

end

