classdef DynamicCloudNetwork < CloudNetwork & DynamicNetwork

    methods
        function this = DynamicCloudNetwork(node_opt, link_opt, VNF_opt, net_opt)
            this@CloudNetwork(node_opt, link_opt, VNF_opt, net_opt);
            this@DynamicNetwork([], [], [], net_opt);
            % <DynamicNetwork> has the same field with <PhysicalNetwork> is treated as an
            % interface (without data member). Thus, through <CloudAccessNetwork> the
            % initilization has finished.
        end

        function sl = AddSlice(this, slice_opt, varargin)
            slice_opt = this.preAddingSlice(slice_opt);
            % We select the method of <DynamicNetwork> to perform adding slice.
            %             AddSlice@CloudNetwork(this, slice_opt, varargin{:});
            sl = AddSlice@DynamicNetwork(this, slice_opt, varargin{:});
        end   
        
        function sl = RemoveSlice(this, arg1)
            sl = RemoveSlice@CloudNetwork(this, arg1);
            this.optimizeResourcePriceNew();
        end
        
    end
    
    methods(Access=protected)
        function tf = onRemovingSlice(this) %#ok<MANU>
            tf = true;
        end
        function sl = createslice(this, slice_opt)
            this.slices{end+1} = DynamicSlice(slice_opt);
            sl = this.slices{end};
        end
        
        % Override the inherited function.
        %         function slice_opt = preAddingSlice(this, slice_opt)
        %         end
        
        function sl = onAddingSlice(this, sl)          
            %             output2 = this.singleSliceOptimization(this.options);
            %             vnf2 = sl.VNFCapacity;
            output = this.optimizeResourcePriceNew;
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
            finalize@CloudNetwork(this, node_price, link_price);
            finalize@DynamicNetwork(this);
        end
    end
    
    methods(Static, Access = protected)
    end
    methods(Access=protected)    
        % overide the default action. 
        % (1) simulate flows originating from un-covered stations, which triggers resource
        %     reallocation of the slice. 
        %     (a) It should be controlled that only a portion of flows will originated
        %         from uncovered areas, other the reoource reallocation will be too
        %         frequent. The portion can be specified as a paramteter of the slice. We
        %         can monitor the actual portion of the unexpected flows, if the portion
        %         supercedes the specified the flow, and regenerate an expected flow.
        function ft = createflow(this, slice, numflow)
            
        end
    end
        
end

