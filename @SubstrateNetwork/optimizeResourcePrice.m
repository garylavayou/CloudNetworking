function [output, runtime] = optimizeResourcePrice(this, slices)
global DEBUG; %#ok<NUSED>
%% Initialization
if nargin == 1
	slices = this.slices;       % all slices are involved in slice dimensioning
end
if nargout == 2 
	options.CountTime = true;
	t_start = tic;
else
	options.CountTime = false;
end
if strcmpi(this.options.Form, 'compact')
	options.bCompact = true;
end
options.OptimizeOrder = 0;
Ns = length(slices);
if nargin == 1
	node_capacity = this.readDataCenter('Capacity');
	link_capacity = this.readLink('Capacity');
else
	% residual capacity + reallocatable capacity.
	node_capacity = this.readDataCenter('ResidualCapacity');
	link_capacity = this.readLink('ResidualCapacity');
	for i = 1:NS
		sl = slices{i};
		node_capacity(sl.getDCPI) = node_capacity(sl.getDCPI) + ...
			sl.VirtualDataCenters.Capacity;
		link_capacity(sl.VirtualLinks.PhysicalLink) = ...
			link_capacity(sl.VirtualLinks.PhysicalLink) + sl.VirtualLinks.Capacity;
	end
end
st = rng;
rng(20180909);
link_uc = 0.5*this.getLinkCost.*((rand(this.NumberLinks,1)-0.5)/100+1);
node_uc = 0.5*this.getNodeCost.*((rand(this.NumberDataCenters,1)-0.5)/100+1);
rng(st);
t = 0;
sp_profit = -inf*ones(3,1);
% fmincon_opt = optimoptions(@fmincon);
% fmincon_opt.Algorithm = 'interior-point';
% fmincon_opt.SpecifyObjectiveGradient = true;
% fmincon_opt.Display = 'notify';     % iter
% fmincon_opt.OptimalityTolerance = 1e-5;
% fmincon_opt.CheckGradients = true;
% fmincon_opt.FiniteDifferenceType = 'central';
%% Parallel initialization
M = getParallelInfo();
array = cell(Ns, 1);
problem = cell(Ns, 1);
indices = cell(Ns, 1);
parfor (pj = 1:Ns,M)
	sl = slices{pj};
	[array{pj}, problem{pj}, indices{pj}] = sl.initializeProblem(options);
end
for j = 1:Ns
	slices{j}.setProblem(array{j}, problem{j}, indices{j}, link_capacity/Ns, node_capacity/Ns);
end
%% Joint Slice Dimensioning Problem: First Stage
% The slice provider continually increases the trail price based on
% the resource cost. The joint slice dimensioning problem then find
% out a social welfare price for slices.
iter_method = 'PartialInverse';		% {DualDecomposition, DualADMM}
while true
	trial_price_node = node_uc * 2^t;
	trial_price_link = link_uc * 2^t;
	switch iter_method
		case 'PartialInverse'
			results = this.SolveSCPPP(slices, trial_price_node, trial_price_link, options);
		case 'DualDecomposition'
			results = this.SolveSCPDD(slices, trial_price_node, trial_price_link, options);
		case 'DualADMM'
			results = this.SolveSCP(slices, trial_price_node, trial_price_link, options);			% return results
		otherwise
			error('error: un-recognized method.');
	end
	sp_profit(1:2) = sp_profit(2:3);
	sp_profit(3) = this.getSliceProviderProfit(results.Prices, struct('PricingPolicy', 'linear'));   % % CloudNetwork.getSliceProviderProfit
	
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
	trial_price_node = node_uc * beta_N;
	trial_price_link = link_uc * beta_N;
	switch iter_method
		case 'PartialInverse'
			results = this.SolveSCPPP(slices, trial_price_node, trial_price_link, options);
		case 'DualDecomposition'
			results = this.SolveSCPDD(slices, trial_price_node, trial_price_link, options);
		case 'DualADMM'
			results = this.SolveSCP(slices, trial_price_node, trial_price_link, options);			% return results
		otherwise
			error('error: un-recognized method.');
	end
	sp_profit_N = this.getSliceProviderProfit(results.Prices, struct('PricingPolicy', 'linear'));
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
output.NodePrice = results.NodePrice;
output.LinkPrice = results.LinkPrice;

%% Finalize substrate network
this.finalize(results.Prices, slices);
if nargout >= 1
	output = this.calculateOutput([], getstructfields(options, {'Slices','PricingPolicy'}));
end
if nargout == 2
	runtime = toc(t_start);
end

end
