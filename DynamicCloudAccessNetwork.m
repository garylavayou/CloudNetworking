%% Dynamic Cloud Access Network 
%% TODO
% (1) simulate handover, which causes flow dynamics; Like the unexpected flows in
%     <DynamicNetwork>, the portion of handover flows should also be controlled with a
%     portion.
classdef DynamicCloudAccessNetwork < CloudAccessNetwork & DynamicNetwork
    methods
        function this = DynamicCloudAccessNetwork(node_opt, link_opt, VNF_opt, net_opt)
            this@CloudAccessNetwork(node_opt, link_opt, VNF_opt, net_opt);
            this@DynamicNetwork([], [], [], net_opt);
            % <DynamicNetwork> has the same field with <PhysicalNetwork> is treated as an
            % interface. Thus, through <CloudAccessNetwork> the initilization has
            % finished.
        end

        function sl = AddSlice(this, slice_opt, varargin)
            slice_opt = this.preAddingSlice(slice_opt);
            sl = AddSlice@DynamicNetwork(this, slice_opt, varargin{:});
        end   
        
        function sl = RemoveSlice(this, arg1)
            sl = RemoveSlice@CloudNetwork(this, arg1);
            this.optimizeResourcePriceNew([], options);
        end
        
    end
    
    methods(Access=protected)
        function tf = onRemovingSlice(this) %#ok<MANU>
            tf = true;
        end
        function sl = createslice(this, slice_opt)
            this.slices{end+1} = DynamicAccessSlice(slice_opt);
            sl = this.slices{end};
        end
        
        function sl = onAddingSlice(this, sl)          
            %             output2 = this.singleSliceOptimization(this.options);
            %             vnf2 = sl.VNFCapacity;
            output = this.optimizeResourcePriceNew([]);
            %             vnf1 = sl.VNFCapacity;
            if isempty(sl)
                this.RemoveSlice(sl);
                sl = sl.empty;
            else
                sl = this.slices{end};
                global g_results event_num; %#ok<TLEV>
                % At the beginning, the slice is added, without consideration of
                % reconfiguration cost.
                g_results.Profit(event_num,1) = output.Profit(1);
                g_results.Solution(event_num,1) = sl.Variables;
                [g_results.Cost(event_num,1),...
                    g_results.NumberReconfig(event_num,1),...
                    g_results.RatioReconfig(event_num,1)]...
                    = sl.get_reconfig_cost;
                g_results.NumberFlows(event_num,1) = sl.NumberFlows;
            end
        end
        %%
        % |finalize| should only be called when dimensiong network slices.
        function finalize(this, node_price, link_price)
            finalize@CloudAccessNetwork(this, node_price, link_price);
            finalize@DynamicNetwork(this);
        end
    end
    
    methods(Static, Access = protected)
        %%%
        % subclass can override this method. Called by _createflow_ method of
        % <DynamicNetwork>. 
        % |slice|: type of DynamicAccessSlice
        % NOTE: this method is derived from <DynamicNetwork>, different from the
        % non-static _updateSliceData_ method, inheirited from <CloudNetwork>.
        function slice_opt = updateDynamicSliceOptions(slice, slice_opt)
            switch slice.Options.FlowPattern
                case FlowPattern.RandomInterDataCenter
                    slice_opt.NodeSet = slice.getDCNI;
                case FlowPattern.RandomInterBaseStation
                    slice_opt.NodeSet = slice.getBSNI;
                case FlowPattern.RandomDataCenter2BaseStation
                    slice_opt.BSNodeSet = slice.getBSNI;
                    slice_opt.DCNodeSet = slice.getDCNI;
                otherwise
                    error('error: cannot handle the flow pattern <%s>.', ...
                        slice.Options.FlowPattern.char);
            end
            slice_opt.MiddleNodes = slice.getDCNI;
        end
    end
    methods(Access=protected)          
       
    end
        
end

