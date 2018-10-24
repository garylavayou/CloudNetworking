%% CloudNetworkEx 
% Including experimental methods.   
classdef CloudNetworkEx < CloudNetwork
    properties (Dependent)
        phis_n;
        phis_l;
    end
    
    properties (SetAccess = private, GetAccess = {?Slice})
        delta = 0.5;  % the proportion of link when computing the network ultilization.
        lambda; % temporarily store the dual variables for the dual decomposition method
        eta;    % proportion of the embedded profit in the unit resource cost, the value is in the range of [0,1)
        static_factor = 1; % If |static_factor = 0|, do not count the static cost.
        enable_constant_profit = false;
    end
    
    methods 
        %%
        %   CloudNetworkEx(node_opt, link_opt, VNF_opt, options)
        % * *options*: 
        %       _ProfitType_ can be |AccuratePrice| or |ApproximatePrice|, used for
        %       calculating profit of Slices and Network;  
        %       _WelfareType_ can be |Accurate| or |Approximate|, used for calculating net
        %       social welfare of Network. 
        function this = CloudNetworkEx(varargin)
            this@CloudNetwork(varargin{:});
            if length(varargin)>=4
                this.options = structmerge(this.options, ...
                    getstructfields(varargin{4}, {'ProfitType', 'WelfareType'}));
            end
            if length(varargin) >= 4
                options = varargin{4};
                if isfield(options, 'Delta')
                    this.delta = options.Delta;
                    options = rmfield(options, 'Delta');
                end
                %%%
                % *Embedded Profit in the Resource Cost*:
                % when the load of network is low, the price the substrate network offer
                % to the slice cannot make profit. Thus this feature can be enabled to
                % guarantee a minimum profit of the substrate network.
                %
                % This feature is enabled by provide the |embed_profit| option in the
                % argument |options|.
                if isfield(options, 'EnableConstantProfit')
                    this.enable_constant_profit = options.EnableConstantProfit;
                    options = rmfield(options, 'EnableConstantProfit');
                end
                if isfield(options, 'EmbedProfit')
                    if ~isempty(options.EmbedProfit)
                        this.eta = options.EmbedProfit;
                        link_uc = this.readLink('UnitCost');
                        this.writeLink('UnitCost', link_uc/(1-this.eta));
                        node_uc = this.readDataCenter('UnitCost');
                        this.writeDataCenter('UnitCost', node_uc/(1-this.eta));
                    end
                    options = rmfield(options, 'EmbedProfit');
                end
                this.options = structmerge(this.options, options);
            end
            % Add slice data
            %             if nargin >=5
            %                 this.writeLink('Beta', beta{1});
            %                 this.writeDataCenter('Beta', beta{2});
            %             else
            %                 this.writeLink('Beta', 1*this.readLink('Capacity'));
            %                 this.writeDataCenter('Beta', 1*this.readDataCenter('Capacity'));
            %             end
            %                 options.beta.node = 1;          % NOTE: this is not used temporarily
            %                 options.beta.link = 1;            
        end
    end
    
    methods
        function p = get.phis_n(this)
            if this.static_factor == 0
                p = 0;
                return;
            end
            epsilon = this.unitStaticNodeCost;
            p = epsilon*this.delta*(this.NumberNodes-this.static_factor)/this.totalNodeCapacity;
        end
        
        function p = get.phis_l(this)
            if this.static_factor == 0
                p = 0;
                return;
            end
            epsilon = this.unitStaticNodeCost;
            p = epsilon*(1-this.delta)*(this.NumberNodes-this.static_factor)/this.totalLinkCapacity;
        end
    end
    
    methods
        function sl = AddSlice(this, slice_opt, varargin)
            if ~this.enable_constant_profit && isfield(slice_opt, 'ConstantProfit')
                slice_opt = rmfield(slice_opt, 'ConstantProfit');
            end
            sl = AddSlice@CloudNetwork(this, slice_opt, varargin{:});
        end
        
        function theta = utilizationRatio(this, load)
            if nargin == 1
                load.Node = this.readDataCenter('Load');
                load.Link = this.readLink('Load');
            end
            theta_v = sum(load.Node)/this.totalNodeCapacity;
            theta_l = sum(load.Link)/this.totalLinkCapacity;
            theta = this.delta*theta_v + (1-this.delta)*theta_l;
        end  
        
        function c = unitStaticNodeCost(this)
            c = mean(this.readDataCenter('StaticCost'));
        end
        
        function c = getStaticCost(this, node_load, link_load)
            if nargin <= 2
                link_load = this.readLink('Load');
            end
            if nargin <=1
                node_load = this.readDataCenter('Load');
            end
            if this.static_factor == 0 
                c = 0;
                return;
            end
            epsilon = this.unitStaticNodeCost;
            theta = this.utilizationRatio(node_load, link_load);
            c = epsilon*((this.NumberNodes-this.static_factor)*theta+this.static_factor);
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
        function c = getNetworkCost(this, node_load, link_load, model)
            if nargin <=1 || isempty(node_load)
                node_load = this.readDataCenter('Load');
            end
            if nargin <= 2 || isempty(link_load)
                link_load = this.readLink('Load');
            end
            if nargin <= 3
                warning('model is set as Approximate.');
                model = 'Approximate';
            end
            
            c = getNetworkCost@CloudNetwork(this, node_load, link_load);
            if strcmp(model, 'Approximate')
                c = c + this.getStaticCost(node_load, link_load);
            elseif strcmp(model, 'Accurate')
                if this.static_factor ~= 0
                    b_deployed = node_load > 0;
                    c = c + sum(this.readDataCenter('StaticCost', b_deployed));
                end
            else
                error('error: unknown model (%s).', model);
            end
        end

        %%% compute link cost. Subclass may override this to provide cost.
        function link_uc = getLinkCost(this)
            link_uc = this.readLink('UnitCost') + this.phis_l;
        end
        
        %%% compute node cost. Subclass may override this to provide cost.
        function node_uc = getNodeCost(this)
            node_uc = this.readDataCenter('UnitCost') + this.phis_n;
        end
    end

    methods (Access = protected)
        function initializeVNF(this, VNF_opt)
            initializeVNF@CloudNetwork(this, VNF_opt);
            switch VNF_opt.Model
                case VNFIntegrateModel.AllInOne
                    if ~isfield(VNF_opt, 'NodeStaticCostOption')
                        VNF_opt.StaticCostOption = NodeStaticCostOption.None;
                    end
                    if VNF_opt.StaticCostOption == NodeStaticCostOption.Random
                        % VNF_opt.static_cost_range = [0.1 0.3];
                        % VNF_opt.ProcessEfficiency is not set, using random value.
                        if ~isfield(VNF_opt, 'StaticCostRange') ||...
                                isempty(VNF_opt.StaticCostRange)
                            scr = [0.1 0.3];
                        else
                            scr = VNF_opt.StaticCostRange;
                        end
                        rng(VNF_opt.RandomSeed(2));
                        node_capacity = this.DataCenters.Capacity;
                        node_uc = this.DataCenters.UnitCost;
                        avg_static_cost = mean(node_capacity.*node_uc)...
                            .*(rand([this.NumberDataCenters, 1])*(scr(2)-scr(1))+scr(1));
                        this.writeDataCenter('StaticCost', round(avg_static_cost));
                    elseif VNF_opt.StaticCostOption == NodeStaticCostOption.NetworkSpecified
                        % TODO
                    elseif VNF_opt.StaticCostOption == NodeStaticCostOption.None
                        this.DataCenters.StaticCost = zeros(this.NumberDataCenters,1);
                        this.static_factor = 0;
                    end
                case VNFIntegrateModel.SameTypeInOne
                    % TODO
                case VNFIntegrateModel.Separated
                    % TODO
            end
        end
                
        % Save the computation results to individual slices and network.
        % TODO: rename as <save_results>.
        function saveStates(this, node_price, link_price, lambda)
            warning('this function should be redesigned.');
            for i = 1:this.NumberSlices
                sl = this.slices{i};
                sl.Variables.x = sl.temp_vars.x;
                sl.Variables.z = sl.temp_vars.z;
                sl.VirtualDataCenters.Load = sl.getNodeLoad;
                sl.VirtualLinks.Load = sl.getLinkLoad;
                sl.FlowTable.Rate = sl.getFlowRate;
                sl.setPathBandwidth;
                if nargin >= 3 && ~isempty(node_price)
                    sl.VirtualDataCenters.Price = node_price(sl.getDCPI);
                end
                if nargin >= 4 && ~isempty(link_price)
                    sl.VirtualLinks.Price = link_price(sl.VirtualLinks.PhysicalLink);
                end
            end
            if nargin >= 4 && ~isempty(lambda)
                this.lambda = lambda;
            end
        end
        function clearStates(this)
            for i = 1:this.NumberSlices
                sl = this.slices{i};
                sl.VirtualDataCenters.Load = zeros(sl.NumberDataCenters,1);
                sl.VirtualDataCenters.Price = zeros(sl.NumberDataCenters,1);
                sl.VirtualLinks.Load = zeros(sl.NumberVirtualLinks,1);
                sl.VirtualLinks.Price = zeros(sl.NumberVirtualLinks,1);
                sl.FlowTable.Rate = zeros(sl.NumberFlows,1);
                sl.setPathBandwidth(zeros(sl.NumberPaths,1));
                sl.Variables.x = [];
                sl.Variables.z = [];
                this.lambda = [];
            end
        end
        %%%
        % used for single slices.
        function slice_data = updateSliceData(this, slice_data)
            slice_data.ConstantProfit = 0;
            if this.options.SlicingMethod == SlicingMethod.SingleFunction
                slice_data.Alpha_f = zeros(NS, 1);
                slice_data.VNFList = 1;
            end
            for s = 1:this.NumberSlices
                sl = this.slices{s};
                slice_data.ConstantProfit = slice_data.ConstantProfit + sl.constant_profit;
                if this.options.SlicingMethod == SlicingMethod.SingleFunction
                    slice_data.Alpha_f(s) = sum(this.VNFTable{sl.VNFList,{'ProcessEfficiency'}});
                end
            end
        end
        
        function argout = calculateOptimalOutput(this, ss)
            tempout = calculateOptimalOutput@CloudNetwork(this, ss);
            if cellstrfind(welfare_type, 'Accurate')
                argout.WelfareAccurateOptimal = tempout.WelfareOptimal...
                    + ss.constant_profit;
            end
            if cellstrfind(welfare_type, 'Approximate')
                argout.WelfareApproxOptimal = ss.getProfit(struct('PricingPolicy', 'linear'));   % net profit
            end
            % if ~isempty(this.eta)
            %     embed_profit_approx = this.eta * ...
            %         this.getNetworkCost(ss.VirtualNodes.Load, ss.VirtualLinks.Load, 'Approximiate');
            %     embed_profit_accurate = this.eta * ...
            %         this.getNetworkCost(ss.VirtualNodes.Load, ss.VirtualLinks.Load, 'Accurate');
            % else
            %     embed_profit_approx = 0;
            %     embed_profit_accurate = 0;
            % end
            % argout.welfare_approx_optimal = argin.welfare_approx_optimal + embed_profit_approx;
            % argout.welfare_accurate_optimal = argin.welfare_accurate_optimal + embed_profit_accurate;
        end
    end
    
    methods
        [output, runtime] = partitionResourcePricing(this, init_price);
        [output, runtime] = optimizeResourcePrice0(this, init_price);
        [output, runtime] = optimizeResourcePrice2(this, init_price);
        [output, runtime] = resourcePartitionOptimization(this, slice_weight);
        [output] = optimizeNetSocialWelfare1( this );
        welfare = optimizeNetSocialWelfare2( this );
        welfare = optimizeNetSocialWelfare2a( this );
        welfare = optimizeNetSocialWelfare3( this );
        welfare = optimizeNetSocialWelfare4( this );
        % SingleSliceOptimization
        % the profit type with |Percent| has been deprecated.
        % if cellstrfind(options.ProfitType, 'Percent')
        %     warning('the profit type with Percent has been deprecated.');
        % end
    end
end

