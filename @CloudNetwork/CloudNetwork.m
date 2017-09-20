%% Cloud Network
% Enable network resource allocation, by mechanisms such as pricing.
%
% In the network, forwarding nodes only take charge of packet processing. Forwarding nodes
% may connect to Data Center, and then VNF instances can be created in data center. The
% connection between forwarding nodes and data centers are assumed with infinity bandwidth
% and zero latency.
% 
% The topology of CloudNetwork is predefined.
%%
classdef CloudNetwork < PhysicalNetwork    
    methods
        %%
        % * *options*: 
        %       _PricingFactor_ is used for <singleSliceOptimization> and <staticSlicing>. 
        %       _Threshold_ is used for resource pricing. 
        %       _Method_ is used for selecting method to solve sub-problem.
        function this = CloudNetwork(varargin)
            this@PhysicalNetwork(varargin{:});
            if length(varargin)>=4
                this.options = structmerge(this.options, ...
                    getstructfields(varargin{4}, {'PricingFactor', 'Threshold', 'Method'}));
            end
        end
    end
    
    methods
        function sl = AddSlice(this, slice_opt, varargin)
            slice_opt = this.preAddingSlice(slice_opt);
            sl = AddSlice@PhysicalNetwork(this, slice_opt, varargin);
        end
        
        %%% compute the total link cost.
        function c = getTotalLinkCost(this, link_load)
            if nargin == 1
                c = dot(this.getLinkField('Load'), this.getLinkCost);
            else
                c = dot(link_load, this.getLinkCost);
            end
        end

        %%% compute the total node cost.
        function c = getTotalNodeCost(this, node_load)
            if nargin == 1
                c = dot(this.DataCenters.Load, this.getNodeCost);
            else
                c = dot(node_load, this.getNodeCost);
            end
        end
        
        %%%
        % * *Network Operation Cost*: 
        % There are two methods to calculate network cost,
        %
        % # Calculate with the approximate model, where the static node cost is computed
        % by the approximate formula.
        % # Calculate with the accurate model, where the static node cost is computed by
        % the solution of VNF deployment.
        %
        % When the network only include a single slice, this method equals to
        % _getSliceCost_ .
        %         function c = getNetworkCost(this, node_load, link_load, model)
        function c = getNetworkCost(this, node_load, link_load)
            if nargin <=1 || isempty(node_load)
                node_load = this.getDataCenterField('Load');
            end
            if nargin <= 2 || isempty(link_load)
                link_load = this.getLinkField('Load');
            end
            
            c = this.getTotalNodeCost(node_load) + this.getTotalLinkCost(link_load);
        end
        
        function theta = utilizationRatio(this, node_load, link_load)
            if nargin == 1
                node_load = this.getDataCenterField('Load');
                link_load = this.getLinkField('Load');
            end
            theta_v = sum(node_load)/this.totalNodeCapacity;
            theta_l = sum(link_load)/this.totalLinkCapacity;
            theta = 0.5*(theta_v + theta_l);
        end
        
        %% statistics of the output
        % type_index is a scalar.
        function [p,r] = statSlice(this, type_index, profit)
            s_index = this.findSlice(type_index);
            if isempty(s_index)
                p = [0, 0, 0, 0];
                if nargout >= 2
                    r = [0, 0, 0, 0];
                end
            else
                %%%
                % Only the statistics of the admitted slices are counted.
                p = [mean(profit(s_index)), max(profit(s_index)), min(profit(s_index)), ...
                    std(profit(s_index))];
                if nargout >= 2
                    rate = zeros(this.NumberFlows,1);
                    i = 0;
                    for s = s_index     % s_index is a row vector
                        num_flow = this.slices{s}.NumberFlows;
                        rate(i + (1:num_flow)) = this.slices{s}.FlowTable.Rate;
                        i = i +num_flow;
                    end
                    rate = rate(1:i);
                    r = [mean(rate), max(rate), min(rate), std(rate)];
                end
            end
        end 

        %%% compute link cost. Subclass may override this to provide cost.
        function link_uc = getLinkCost(this, link_id)
            if nargin == 1
                link_uc = this.getLinkField('UnitCost');
            else
                link_uc = this.getLinkField('UnitCost', link_id);
            end
        end

        %%% compute node cost. Subclass may override this to provide cost.
        % * *dc_id*: data center index (not the node index of the substrate physical node).
        function node_uc = getNodeCost(this, dc_id)
            if nargin == 1
                node_uc = this.getDataCenterField('UnitCost');
            else
                node_uc = this.getDataCenterField('UnitCost', dc_id);
            end
        end
    end
    
    methods (Access=protected) 
        function graph = residualgraph(this, slice_opt)
            if isfield(slice_opt, 'method') && strcmp(slice_opt.method, 'static-slicing')
                % If a link's residual capacity is zero, then this link should be removed
                % from the grpah.
                % If a node's residual capacity is zero, then this node and the adjacent
                % links should be removed from the graph.
                link_capacity = this.getLinkField('ResidualCapacity');
                node_capacity = this.getDataCenterField('ResidualCapacity');
                link_capacity(link_capacity<1) = 0;
                node_capacity(node_capacity<1) = 0;
                A = spalloc(this.NumberNodes, this.NumberNodes, this.NumberLinks);
                C = spalloc(this.NumberNodes, this.NumberNodes, this.NumberLinks);
                for i = 1:this.NumberLinks
                    if link_capacity(i) <= 0
                        continue;
                    end
                    h = this.graph.Head(i);
                    dc_h = this.Topology.Nodes.DataCenter(h);
                    if dc_h && node_capacity(dc_h) <= 0
                        continue;
                    end
                    t = this.graph.Tail(i);
                    dc_t = this.Topology.Nodes.DataCenter(t);
                    if dc_t && node_capacity(dc_t) <= 0
                        continue;
                    end
                    A(h, t) = this.graph.Adjacent(h, t); %#ok<SPRIX>
                    C(h, t) = link_capacity(i); %#ok<SPRIX>
                end
                graph = DirectedGraph(A, C);
            else
                graph = this.graph;
            end
        end
        
        % Implement the abstract function.
        function slice_opt = preAddingSlice(this, slice_opt) %#ok<INUSL>
            if ~isfield(slice_opt,'Weight') || isempty(slice_opt.Weight) ...
                    || slice_opt.Weight == 0
                error('error: invalid slice weight.'); %     slice_opt.Weight = 1;
            end
            if ~isfield(slice_opt, 'method') || isempty(slice_opt.method)
                slice_opt.method = 'dynamic-slicing';
            end
        end
        
        function sl = createslice(this, slice_opt, varargin)
            this.slices{end+1} = Slice(slice_opt);
            sl = this.slices{end};
        end
        %%%
        % calculate the profit ratio of slices and network
        % 
        function [b, profit_gap] = checkProfitRatio(this, node_price, link_price, options)
            slice_profit_ratio = zeros(this.NumberSlices,1);
            if nargin <= 3 || ~isfield(options, 'PricingPolicy')
                warning('pricing policy not set.');
                options.PricingPolicy = '';
            end
            if ~isfield(options, 'Threshold')
                warning('profit threshold not set.');
                options.Threshold = 'min';
            end
            for s = 1:this.NumberSlices
                sl = this.slices{s};
                revenue = sl.getRevenue;
                options.LinkPrice = link_price(sl.VirtualLinks.PhysicalLink);
                options.NodePrice = node_price(sl.getDCPI);
                % Prices have been announced to each slice.
                slice_profit_ratio(s) = Slice.fcnProfit([], sl, options)/revenue;
            end
            clear revenue;
            
            [sp_profit, sp_revenue] = this.getSliceProviderProfit(...
                node_price, link_price, options.PricingPolicy);
            network_profit_ratio = sp_profit/sp_revenue;
            % a = 1;        % {0.5|0.75|1}
            switch options.Threshold
                case 'min'
                    profit_threshold = min(slice_profit_ratio);
                case 'average'
                    profit_threshold = mean(slice_profit_ratio);
                case 'max'
                    profit_threshold = max(slice_profit_ratio);
                otherwise
                    error('error: invalid option (Threshold = %s)', options.Threshold);
            end
            if nargin == 4  && isfield(options, 'epsilon')
                if abs(network_profit_ratio - profit_threshold) < options.epsilon
                    b = true;
                else
                    b = false;
                end
            else
                if network_profit_ratio >= profit_threshold
                    b = true;
                else
                    b = false;
                end
            end
            if nargout == 2
                profit_gap = network_profit_ratio - profit_threshold;
            end
            if nargin== 4 && isfield(options, 'Display') && ...
                    ~isempty(strfind(options.Display, 'iter'))
                disp('Profit ratio {slices|network}:');
                disp([slice_profit_ratio; network_profit_ratio]);
            end
        end
           
        %%%
        % * *getSliceProviderProfit*
        % If |node_price| and |link_price| are not provided, the stored price are used.
        function [profit, revenue] = ...
                getSliceProviderProfit(this, node_price, link_price, pricing_policy)
            revenue = 0;
            [node_load, link_load] = this.getNetworkLoad;
            if isempty(node_price) || isempty(link_price)
                node_price = this.getDataCenterField('Price');
                link_price = this.getLinkField('Price');
            end                
            if nargin == 4
                if strcmpi(pricing_policy, 'quadratic-price')
                    for s = 1:this.NumberSlices
                        sl = this.slices{s};
                        link_id = sl.VirtualLinks.PhysicalLink;
                        dc_id = sl.getDCPI;
                        revenue = revenue + ...
                            sl.fcnLinkPricing(link_price(link_id), sl.getLinkLoad(sl.x_path)) + ...
                            sl.fcnNodePricing(node_price(dc_id), sl.getNodeLoad(sl.z_npf));
                    end
                end
            elseif nargin == 3
                revenue = dot(node_load, node_price) + dot(link_load, link_price);
            else
                error('error: input parameters are not enough.');
            end
            profit = revenue - this.getNetworkCost(node_load, link_load);
        end
        
        %%%
        % * *Finalize substrate network*
        %
        % # Record the resource allocation variables, flow rate, virtual node/link load of
        %   each slice.
        % # Virtual Nodes/Links' capacity is derived from node/link load;
        % # Calculate and announce the resource prices to each slice.
        % # Record the substrate network's node/link load, price.
        %
        % Usually, this function should be provided with 3 arguments, except that it is
        % called by
        % <file:///E:/workspace/MATLAB/Projects/Documents/CloudNetworking/singleSliceOptimization.html singleSliceOptimization>.
        % NOTE: the price here might be only prcing parameters (for varing pricing
        % policy). To calculate the payment, using _fcnLinkPricing_ and _fcnNodePricing_
        % function.
        function finalize(this, node_price, link_price)
            for i = 1:this.NumberSlices
                sl = this.slices{i};
                sl.Variables.x = sl.x_path;
                sl.Variables.z = sl.z_npf;
                sl.VirtualDataCenters.Load = sl.getNodeLoad;
                sl.VirtualLinks.Load = sl.getLinkLoad;
                sl.VirtualDataCenters.Capacity = sl.VirtualDataCenters.Load;
                sl.VirtualLinks.Capacity = sl.VirtualLinks.Load;
                sl.FlowTable.Rate = sl.getFlowRate;
                sl.setPathBandwidth;
                if nargin >= 2 && ~isempty(node_price)
                    sl.VirtualDataCenters.Price = node_price(sl.getDCPI);
                end
                if nargin >= 3 && ~isempty(link_price)
                    sl.VirtualLinks.Price = link_price(sl.VirtualLinks.PhysicalLink);
                end
            end
            [node_load, link_load] = this.getNetworkLoad;
            this.setLinkField('Load', link_load);
            this.setDataCenterField('Load', node_load);
            if nargin == 3
                this.setLinkField('Price', link_price);
                this.setDataCenterField('Price', node_price);
            end
        end

        function options = updatePathConstraints(this, slice_opt)
            options = this.updatePathConstraints@PhysicalNetwork(slice_opt);
            if this.NumberDataCenters < this.NumberNodes
                % if only part of the forwarding nodes is VNF-capable, we should make sure that the
                % path at least transit one VNF-capable node.
                % no matter when, the DataCenters is the middle nodes. However, if the MiddleNodes
                % option is not provided, the route calculation will be performed in a different way.
                options.MiddleNodes = this.DataCenters.NodeIndex;
            end
        end
        
        function slice_data = updateSliceData(this, slice_data, options)
            if strcmp(options.Method, 'normal')
                slice_data.VNFList = 1:this.NumberVNFs;
                %         b_vnf(this.slices{s}.VNFList) = true;
                %     slice_data.VNFList = find(b_vnf);
            end
        end
        function output = calculateOptimalOutput(this, ss, options, ~)
            output.welfare_optimal = sum(...
                ss.FlowTable.Weight.*fcnUtility(ss.FlowTable.Rate)) ...
                - this.getNetworkCost(ss.VirtualDataCenters.Load, ss.VirtualLinks.Load);
            if ~strcmp(options.Display, 'off') && ~strcmp(options.Display, 'none')
                fprintf('\tThe optimal net social welfare of the network: %G.\n', ...
                    output.welfare_optimal);
            end
            output.z_npf = reshape(full(ss.Variables.z), ...
                this.NumberDataCenters, this.NumberPaths, this.NumberVNFs);
        end
        argout = calculateOutput(this, argin, options);    
        [tb, stbs] = saveStatTable(PN, output, rt, slice_types, method);
        
        % Now the same as <PhysicalNetwork>, subclasses might override it.
        %         function [flow_table, phy_adjacent, flag] = ...
        %                 generateFlowTable(this, graph, slice_opt)
        %             [flow_table, phy_adjacent, flag] = ...
        %                 this.generateFlowTable@PhysicalNetwork(graph, slice_opt);
        %         end
    end
    
    methods (Access = private)
        [node_price, link_price, runtime] = pricingFactorAdjustment(this, options);
        function runtime = priceIteration(this, node_price, link_price, options)
            if nargout == 1
                slice_runtime = 0;
                runtime.Serial = 0;
            end
            for s = 1:this.NumberSlices
                sl = this.slices{s};
                link_id = sl.VirtualLinks.PhysicalLink;
                dc_id = sl.getDCPI;
                sl.VirtualLinks.Price = link_price(link_id);
                sl.VirtualDataCenters.Price = node_price(dc_id);
                %%%
                % optimize each slice with price and resource constraints.
                if nargout == 1
                    tic;
                end
                sl.optimalFlowRate(options);
                if nargout == 1
                    t = toc;
                    slice_runtime = max(slice_runtime, t);
                    runtime.Serial = runtime.Serial + t;
                end
            end
            if nargout == 1
                runtime.Parallel = slice_runtime;
            end
        end
    end
    
    methods (Static, Access=protected)
        function flag = assert_path_list(end_points, path_list,slice_opt)
            if isempty(path_list)
                if strcmp(slice_opt.method, 'static-slicing')
                    % two choice: reject the slice or reject the flow.
                    if isfield(slice_opt, 'admit_ploicy') && ...
                            strcmp(slice_opt.admit_ploicy, 'reject-slice')
                        warning('Reject the slice request.');
                        flag = 2;
                    else
                        % slice_opt.admit_ploicy = 'reject-slice', the actual number
                        % of generated flow may less than |number_flow|
                        warning('Reject the flow (%d,%d) in the slice request.',...
                            end_points(1), end_points(2));
                        flag = 1;
                    end
                else
                    flag = -1;
                end
            else
                flag = 0;
            end
        end
        
        [stat, slice_stat] = createStatTable(num_point, num_type, type);
    end
    
    methods
        [output, runtime] = optimizeResourcePrice(this, init_price, options);
        [output, runtime] = optimizeResourcePriceNew(this, init_price, options);
        [output, runtime] = singleSliceOptimization(this , options);
        output = StaticSlicing(this, slice, options);
    end
end

