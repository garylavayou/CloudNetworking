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
function [profit, gd] = fcnSocialWelfare(vars, this, options)
if strcmp(options.CostModel, 'fixcost')
	if options.bCompact
		full_vars = sparse(options.num_orig_vars,1);
		full_vars(this.I_active_variables) = vars;
		vars = full_vars;
	end
	var_path = vars(1:this.pardata.NumberPaths);
	flow_rate = this.getFlowRate(var_path);
	profit = -sum(this.pardata.Weight*fcnUtility(flow_rate));
	
	if options.bFinal
		profit = -profit;
	else
		gd = spalloc(length(vars),1, this.pardata.NumberPaths);
		for p = 1:this.pardata.NumberPaths
			i = this.pardata.PathOwner(p);  % find the the owner (flow) of the path.
			% f = S.I_flow_path(:,p)~=0;
			gd(p) = -this.pardata.Weight/(1+this.I_flow_path(i,:)*var_path); %#ok<SPRIX>
		end
		if options.bCompact
			gd = gd(this.I_active_variables);
		end
	end
else
	if nargout <= 1
		profit = fcnSocialWelfare@SimpleSliceOptimizer(vars, this, options);
	else
		[profit, gd] = fcnSocialWelfare@SimpleSliceOptimizer(vars, this, options);
	end
end
end