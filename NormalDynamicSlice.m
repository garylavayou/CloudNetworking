classdef NormalDynamicSlice < NormalSlice & IDynamicSlice

	properties
		Property1
	end
	
	methods
		function this = NormalDynamicSlice(slice_data)
			this@NormalSlice(slice_data);
			this@IDynamicSlice(slice_data);
		end
		
		function finalize(this, prices)
			global DEBUG g_results;
			finalize@NormalSlice(this, prices);
			%% TODO
			
		end
		
		function c = getLinkCapacity(this, isfinal)
			if nargin == 1 || isfinal
				c = this.Links.Capacity;
			else
				if this.invoke_method == 0
					c = getLinkCapacity@NormalSlice(this, isfinal);
				else
					%% TODO c = this.temp_vars.c;
				end
			end
		end
		
		function c = getNodeCapacity(this, isfinal) % isfinal=true by default
			if nargin == 1 || isfinal
				c = this.ServiceNodes.Capacity;
			else
				if this.invoke_method == 0
					c = getNodeCapacity@NormalSlice(this, isfinal);
				else
					%% TODO
				end
			end
		end
		
		function ds = diffstate(this, isfinal)
			if isfinal && ~isempty(this.diff_state)
				if nargout == 1
					ds = this.diff_state;
				end
				return;
			end
			%% TODO
		end
		
		function [total_cost, reconfig_cost] = get_reconfig_cost(this, model, isfinal)
			if nargin <= 1
				model = 'const';
			end
			if strcmpi(model, 'const')
				isfinal = true;
			elseif strcmpi(model, 'linear') && nargin <= 2
				isfinal = false;
			end
			options = getstructfields(this.Parent.options, ...
				{'DiffNonzeroTolerance', 'NonzeroTolerance'});
			tol_vec = options.DiffNonzeroTolerance;
			s = this.diffstate(isfinal);
			
			%% TODO
		end
		
		function stat = get_reconfig_stat(this, stat_names)
			if nargin == 1
				stat_names = {'All'};
			elseif ischar(stat_names)
				stat_names = {stat_names};
			end
			options = getstructfields(this.Parent.options, ...
				{'DiffNonzeroTolerance', 'NonzeroTolerance'});
			s = this.diffstate(true);
			
			stat = table;
			tol_vec = options.DiffNonzeroTolerance;
			for i = 1:length(stat_names)
				%% TODO
			end
		end
		
		function s = save_state(this)
			%% TODO
		end
		
		function set_state(this)
			%% TODO
		end
		
		function [profit, cost] = handle_zero_flow(this, new_opts)
			%% TODO
		end
		
		[utility, node_load, link_load] = priceOptimalFlowRate(this, x0, new_opts);
		
		function [profit, cost] = optimalFlowRate( this, new_opts )
			%% TODO
		end
	
		[profit,cost] = fastReconfigure(this, action, options);
	end
end

