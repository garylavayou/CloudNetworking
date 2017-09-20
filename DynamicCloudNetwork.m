classdef DynamicCloudNetwork < CloudNetwork & DynamicNetwork
    methods
        function this = DynamicCloudNetwork(node_opt, link_opt, VNF_opt, net_opt)
            this@CloudNetwork(node_opt, link_opt, VNF_opt, net_opt);
            %    this@DynamicNetwork(node_opt, link_opt, VNF_opt, net_opt);
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
            this.optimizeResourcePriceNew([], options);
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
            %             vnf2 = sl.Vnf;
            options = getstructfields(this.options, {'Display', 'Threshold'});
            output = this.optimizeResourcePriceNew([], options);
            %             vnf1 = sl.Vnf;
            if isempty(sl)
                this.RemoveSlice(sl);
                sl = sl.empty;
            else
                sl = this.slices{end};
                global g_results event_num; %#ok<TLEV>
                % At the beginning, the slice is added, without consideration of
                % reconfiguration cost.
                g_results.profit(event_num,1) = output.profit(1);
                g_results.solution(event_num,1) = sl.Variables;
                [g_results.cost(event_num,1),...
                    g_results.num_reconfig(event_num,1),...
                    g_results.rat_reconfig(event_num,1)]...
                    = sl.get_reconfig_cost;
                g_results.num_flows(event_num,1) = sl.NumberFlows;
            end
        end
        %%
        % |finalize| should only be called when dimensiong network slices.
        function finalize(this, node_price, link_price)
            finalize@CloudNetwork(this, node_price, link_price);
            for i = 1:this.NumberSlices
                sl = this.slices{i};
                sl.setVnfCapacity;
                %                 if strcmp(this.options.PricingPolicy, 'quadratic-price')
                [~, sl.edge_reconfig_cost] = sl.fcnLinkPricing(...
                    sl.VirtualLinks.Price, sl.VirtualLinks.Load);
                sl.edge_reconfig_cost = DynamicSlice.THETA*sl.edge_reconfig_cost;
                % here the |node_price| is the price of all data centers.
                [~, sl.node_reconfig_cost] = sl.fcnNodePricing(...
                    sl.VirtualDataCenters.Price, sl.VirtualDataCenters.Load);
                sl.node_reconfig_cost = DynamicSlice.THETA*sl.node_reconfig_cost;
                %                 else
                %                 end
            end
        end
    end
    
    methods(Static, Access = protected)
    end
    methods(Access=protected)                 
    end
        
end

