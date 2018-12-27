function [exitflag, fidx] = executeMethod(this, action)
global g_results event_num DEBUG INFO;  %#ok<NUSED>
this.temp_vars = struct();
this.diff_state = struct([]);
fidx = [];
slice = this.hs;
tol_zero = this.options.NonzeroTolerance;
this.max_flow_rate = slice.FlowTable.Rate;
this.max_flow_rate(this.max_flow_rate<tol_zero*max(this.max_flow_rate)) = inf;

%% [TODO] reinitialize the state variables
t1 = tic;
switch this.options.ReconfigMethod
	case ReconfigMethod.Baseline
		% provide 'method' and 'model' to customize the <optimalFlowRate>
		% Since we adopt |FixedCost| model, the resource cost that is a
		% constant, is not include in the objective value |profit|. So, we
		% should exclude it from the |profit|, as well as the reconfiguration
		% cost.
		this.setProblem('LinkPrice', slice.Links.Price, ...
			'NodePrice', slice.ServiceNodes.Price);
		options = Dictionary(...
			'CostModel', 'fixcost', ...
			'isFinalize', true, ...
			'bInitialize', true);
		profit = this.optimalFlowRate(options);
	case {ReconfigMethod.Fastconfig,ReconfigMethod.FastconfigReserve}
		profit = this.fastReconfigure(action);
	case {ReconfigMethod.Fastconfig2, ReconfigMethod.Fastconfig2Reserve}
		error('error: not implemented!');
	case {ReconfigMethod.Dimconfig, ReconfigMethod.DimconfigReserve, ...
			ReconfigMethod.DimBaseline}
		profit = this.DimensioningReconfigure(action);
	otherwise
		error('NetworkSlicing:UnsupportedMethod', ...
			'error: unsupported method (%s) for network slicing.', ...
			this.options.ReconfigMethod.char) ;
end
t2 = toc(t1);
if this.options.ReconfigMethod>=ReconfigMethod.Dimconfig && isempty(profit)
	exitflag = -1;
	return;
end
this.max_flow_rate(this.max_flow_rate==inf) = 0;
this.max_flow_rate = max(slice.FlowTable.Rate, this.max_flow_rate);

% Check resource utilization
%{
idx = this.VirtualLinks.Capacity~=0;
link_util = this.VirtualLinks.Load(idx)./this.VirtualLinks.Capacity(idx);
sum_link_util = sum(this.VirtualLinks.Load(idx))./sum(this.VirtualLinks.Capacity(idx));
mean_link_util = mean(link_util);
%}

% stat.Solution = this.Variables;
stat = slice.get_reconfig_stat();   % see <DynamicSlice>.<get_reconfig_stat>
stat.Profit = profit;
if ~isempty(DEBUG) && DEBUG
	disp(stat);
	fprintf('%s elapsed time: %ds\n',this.options.ReconfigMethod.char, t2);
end
switch this.options.ReconfigMethod
	case ReconfigMethod.Baseline
		stat.ReconfigType = ReconfigType.Reconfigure;
	case {ReconfigMethod.Fastconfig,ReconfigMethod.FastconfigReserve,...
			ReconfigMethod.Fastconfig2, ReconfigMethod.Fastconfig2Reserve}
		stat.ReconfigType = ReconfigType.FastReconfigure;
	case {ReconfigMethod.Dimconfig, ReconfigMethod.DimconfigReserve}
		stat.ReconfigType = ReconfigType.Dimensioning;
		if ~this.b_dim
			stat.ReconfigType = ReconfigType.FastReconfigure;
		elseif slice.options.Adhoc
			% If the slice do not support Adhoc flows, then we do not release resource
			% descriptors for the slice.
			slice.release_resource_description();
		end
	case ReconfigMethod.DimBaseline
		if this.b_dim
			if slice.options.Adhoc
				slice.release_resource_description();
			end
			stat.ReconfigType = ReconfigType.Dimensioning;
		else
			stat.ReconfigType = ReconfigType.Reconfigure;
		end
end
if this.ENABLE_DYNAMIC_NORMALIZER
	this.postl1normalizer();
	if ~isempty(DEBUG) && DEBUG
		disp(this.getbeta);
	end
end
g_results(event_num,:) = stat;

%% Reset temporary variables
% These variables should be cleared. If it is used next time, it will be assigned with new
% values. Some methods will execute according to the state (if it is initialized) of the
% variables.
this.setProblem('Price', []);
this.topts = struct();
this.changed_index = struct();
slice.net_changes.erase();
this.b_dim = false;
exitflag = 0;
end
