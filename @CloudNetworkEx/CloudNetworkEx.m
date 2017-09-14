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
        %   FullCloudNetwork(node_opt, link_opt, VNF_opt, options)
        function this = CloudNetworkEx(varargin)
            this@CloudNetwork(varargin{:});
            if length(varargin) >= 4
                options = varargin{4};
                if isfield(options, 'delta')
                    this.delta = options.delta;
                    options = rmfield(options, 'delta');
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
                if isfield(options, 'embed_profit')
                    if ~isempty(options.embed_profit)
                        this.eta = options.embed_profit;
                        link_uc = this.getLinkField('UnitCost');
                        this.setLinkField('UnitCost', link_uc/(1-this.eta));
                        node_uc = this.getDataCenterField('UnitCost');
                        this.setDataCenterField('UnitCost', node_uc/(1-this.eta));
                    end
                    options = rmfield(options, 'embed_profit');
                end
                this.options = options;
            end
            % Add slice data
            %             if nargin >=5
            %                 this.setLinkField('Beta', beta{1});
            %                 this.setDataCenterField('Beta', beta{2});
            %             else
            %                 this.setLinkField('Beta', 1*this.getLinkField('Capacity'));
            %                 this.setDataCenterField('Beta', 1*this.getDataCenterField('Capacity'));
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
        
        function theta = utilizationRatio(this, node_load, link_load)
            if nargin == 1
                node_load = this.getDataCenterField('Load');
                link_load = this.getLinkField('Load');
            end
            theta_v = sum(node_load)/this.totalNodeCapacity;
            theta_l = sum(link_load)/this.totalLinkCapacity;
            theta = this.delta*theta_v + (1-this.delta)*theta_l;
        end  
        
        function c = unitStaticNodeCost(this)
            c = mean(this.getDataCenterField('StaticCost'));
        end
        
        function c = getStaticCost(this, node_load, link_load)
            if nargin <= 2
                link_load = this.getLinkField('Load');
            end
            if nargin <=1
                node_load = this.getDataCenterField('Load');
            end
            if this.static_factor == 0 
                c = 0;
                return;
            end
            epsilon = this.unitStaticNodeCost;
            theta = this.utilizationRatio(node_load, link_load);
            c = epsilon*((this.NumberNodes-this.static_factor)*theta+this.static_factor);
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
                        if ~isfield(VNF_opt, 'static_cost_range') ||...
                                isempty(VNF_opt.static_cost_range)
                            scr = [0.1 0.3];
                        else
                            scr = VNF_opt.static_cost_range;
                        end
                        rng(VNF_opt.RandomSeed(2));
                        node_capacity = this.DataCenters.Capacity;
                        node_uc = this.DataCenters.UnitCost;
                        avg_static_cost = mean(node_capacity.*node_uc)...
                            .*(rand([this.NumberDataCenters, 1])*(scr(2)-scr(1))+scr(1));
                        this.setDataCenterField('StaticCost', round(avg_static_cost));
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
        
        [node_price, link_price, runtime] = pricingFactorAdjustment(this, options);

        % Save the computation results to individual slices and network.
        % TODO: rename as <save_results>.
        function saveStates(this, node_price, link_price, lambda)
            warning('this function should be redesigned.');
            for i = 1:this.NumberSlices
                sl = this.slices{i};
                sl.Variables.x = sl.x_path;
                sl.Variables.z = sl.z_npf;
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
    end
    
    methods
        [output, runtime] = partitionResourcePricing(this, init_price, options);
        [output, runtime] = optimizeResourcePrice0(this, init_price, options);
        [output, runtime] = optimizeResourcePrice2(this, init_price, options);
        [output, runtime] = resourcePartitionOptimization(this, slice_weight, options);
        [output] = optimizeNetSocialWelfare1( this, options );
        welfare = optimizeNetSocialWelfare2( this, options );
        welfare = optimizeNetSocialWelfare2a( this );
        welfare = optimizeNetSocialWelfare3( this );
        welfare = optimizeNetSocialWelfare4( this );
    end
end

