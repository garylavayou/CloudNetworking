%% Net social welfare of a slice
% override <Slice.fcnSocialWelfare>.
% Slice resource cost is fixed, since the price is given and the resource amount is known.
% When using the |'fixcost'| option, the slice's resource cost is
% constant, and we do not consider it in the objective and thus in the gradient.
function [profit, grad] = fcnSocialWelfare(var_x, S, model)
if strcmp(model, 'fixcost')
    var_path = var_x(1:S.NumberPaths);
    flow_rate = S.getFlowRate(var_path);
    weight = S.weight*ones(S.NumberFlows, 1);
    profit = -sum(weight.*fcnUtility(flow_rate));

    if nargout <= 1
        profit = -profit;
    else
        grad = spalloc(length(var_x),1, S.NumberPaths);
        for p = 1:S.NumberPaths
            i = S.path_owner(p);
            f = S.I_flow_path(:,p)~=0;  % find the weight of the flow on the path.
            grad(p) = -weight(f)/(1+S.I_flow_path(i,:)*var_path); %#ok<SPRIX>
        end
    end
else
    if nargout <= 1
        profit = fcnSocialWelfare@Slice(var_x, S, model);
    else
        [profit, grad] = fcnSocialWelfare(var_x, S, model);
    end
end
end