%% Dynamic Cloud Access Network 
%% TODO
% (1) simulate handover, which causes flow dynamics; Like the unexpected flows in
%     <DynamicNetwork>, the portion of handover flows should also be controlled with a
%     portion.
classdef SimpleDynamicAccessNetwork < IAccessNetwork & SimpleDynamicNetwork
    methods
        function this = SimpleDynamicAccessNetwork(node_opt, link_opt, VNF_opt, net_opt)
            this@SimpleDynamicNetwork(node_opt, link_opt, VNF_opt, net_opt);
            this@IAccessNetwork(node_opt, link_opt, VNF_opt, net_opt);
        end

        function sl = AddSlice(this, slice_opt, varargin)
            slice_opt = this.preAddingSlice(slice_opt);
            sl = AddSlice@SimpleDynamicNetwork(this, slice_opt, varargin{:});
            this.onAddingSlice(sl);
        end   
        
        function sl = RemoveSlice(this, arg1)
            sl = RemoveSlice@SimpleDynamicNetwork(this, arg1);
            this.optimizeResourcePrice();
        end
        
    end
    
    methods(Access=protected)
        function tf = onRemovingSlice(this) %#ok<MANU>
            tf = true;
        end
        function sl = createslice(this, slice_opt)
            this.slices(end+1) = DynamicAccessSlice(slice_opt);
            sl = this.slices(end);
						sl.getOptimizer(slice_opt);
        end
        function slice_opt = preAddingSlice(this, slice_opt)  
					slice_opt = preAddingSlice@SimpleDynamicNetwork(this, slice_opt);
					slice_opt = preAddingSlice@IAccessNetwork(this, slice_opt);
        end
        function sl = onAddingSlice(this, sl)          
            %             output2 = this.singleSliceOptimization(this.options);
            %             vnf2 = sl.VNFCapacity;
            output = this.optimizeResourcePrice();
            %             vnf1 = sl.VNFCapacity;
            if isempty(sl)
                this.RemoveSlice(sl);
                sl = sl.empty;
            else
                sl = this.slices(end);
                global g_results event_num; %#ok<TLEV>
                % At the beginning, the slice is added, without consideration of
                % reconfiguration cost.
                g_results.Profit(event_num,1) = output.Profit(1);
                g_results.Solution(event_num,1) = sl.Variables;
                [g_results.Cost(event_num,1),...
                    g_results.NumberReconfig(event_num,1),...
                    g_results.RatioReconfig(event_num,1)]...
                    = sl.Optimizer.get_reconfig_cost('const');
                g_results.NumberFlows(event_num,1) = sl.NumberFlows;
            end
        end
        %%
        % |finalize| should only be called when dimensiong network slices.
        function finalize(this, prices)
            finalize@SimpleDynamicNetwork(this, prices);
        end
    end
    
    methods(Access = protected)
        %%%
				% Override <DynamicNetwork>.<updateDynamicSliceOptions>, different from the
        % <NetworkOptimizer>.<updateSliceData> method.
        % Called by _createflow_ method of <DynamicNetwork>. 
        % |slice|: type of DynamicAccessSlice
        function slice_opt = updateDynamicSliceOptions(this, slice, slice_opt)
					switch slice.options.FlowPattern
						case {FlowPattern.RandomInterDataCenter, FlowPattern.RandomInterBaseStation}
							slice_opt.NodeSet = slice.options.NodeSet;
						case FlowPattern.RandomDataCenter2BaseStation
							slice_opt.BSNodeSet = slice.options.BSNodeSet;
							slice_opt.DCNodeSet = slice.options.DCNodeSet;
						otherwise
							error('error: cannot handle the flow pattern <%s>.', ...
								slice.options.FlowPattern.char);
					end
					structmerge(slice_opt, getstructfields(slice.options, 'MiddleNodes', ...
						'default-ignore', this.DataCenters.Node));
				end
    end
        
end

