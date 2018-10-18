classdef (Abstract, HandleCompatible) IDynamicReconfigure
	methods
		ds = diffstate(this, isfinal);
		[total_cost, reconfig_cost] = get_reconfig_cost(this, model, isfinal);
		stat = get_reconfig_stat(this, stat_names);
		s = get_state(this);
		set_state(this);
		[profit, cost] = handle_zero_flow(this, new_opts);
		[profit,cost] = fastReconfigure(this, action, options);
		identify_change(this, changed_path_index);
		update_reconfig_cost(this, action, bDim);
	end
	
end