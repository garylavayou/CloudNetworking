%% Net social welfare of a slice
% override <Slice.fcnSocialWelfare>.
% Slice resource cost is fixed, since the price is given and the resource amount is known.
% When using the |'fixcost'| option, the slice's resource cost is constant. Therefore, we
% do not include it in the objective as well as in the gradient, to reduce computation.
% However, when comparing the profit with other methods, the fixed resource cost should be
% considered.
%
% options:
%   # bCompact: 'true' for compact mode, which reduce the problem scale (not consider
%     those in-active variables corresponding to $b_np=0$).
%   # bFinal: set to 'true' return the real profit (max f => min -f).
function [profit, grad] = fcnSocialWelfare(vars, this, options)
if strcmp(options.CostModel, 'fixcost')
	if isfield(options, 'unit')
		unit = options.unit;
	else
		unit = 1;
	end
	if isempty(this.hs)
		pardata = this.pardata;
	else
		pardata = this.hs;
	end
	% 	if options.bCompact
	% 		full_vars = sparse(options.num_orig_vars,1);
	% 		full_vars(this.I_active_variables) = vars;
	% 		vars = full_vars;
	% 	end
	prbm = this.problem;
	var_r_index = sum(prbm.num_vars(1:2))+(1:sum(prbm.num_vars(3)));
	flow_rate = vars(var_r_index);
	profit = -pardata.Weight*sum(fcnUtility(unit*flow_rate));
	
	if options.bFinal
		profit = -profit;
	else
		grad = sparse(numel(vars),1);
		%% partial derviatives on r
		grad(var_r_index) = -pardata.Weight./(1/unit+flow_rate);
		% 		if options.bCompact
		% 			gd = gd(this.I_active_variables);
		% 		end
	end
else
	if nargout <= 1
		profit = fcnSocialWelfare@NormalSliceOptimizer(vars, this, options);
	else
		[profit, grad] = fcnSocialWelfare@NormalSliceOptimizer(vars, this, options);
	end
end
end