classdef Slice < VirtualNetwork
	properties
		Type;           % Type from slice template
		weight;         % weight for the slice user's utility.
		
		FlowTable;      % flow information in the slice
		VNFList;        % List of virtual network functions in the slice

    b_final logical = false;
  end
	
	properties (Dependent)
    LinkCost;
    NodeCost;
		NumberFlows;					% Number of flows in the slice
		NumberPaths;					% Number of paths in the slice
		NumberVNFs;						% Number of virtual network functions in the slice
	end
	
	properties (SetAccess = {?Slice, ?SliceOptimizer})
		path_owner; % Associated flow of path, to provide fast inquiry of associated flow than using |I_flow_path|.
	end
	
  methods (Abstract)
    profit = getProfit(slice, options);
    r = getRevenue(this);
		vc = getVNFCapacity(this, ~);
		sc = getCost(this, load, model);
	end
    	
  %% Constructor and Destructor
  methods
    function this = Slice(slice_data)
			if nargin == 0
				args = cell(0);
			else
				args = {slice_data};
			end
			this@VirtualNetwork(args{:});
			if nargin == 0
				return;
			end
			
      if isfield(slice_data, 'Type')
        this.Type = slice_data.Type;
      end
      if isfield(slice_data,'Weight')
        this.weight = slice_data.Weight;
      end
      %%%
      % Link price
      if isfield(slice_data, 'LinkPrice')
        this.Links.Price = slice_data.LinkPrice;
      else
        this.Links.Price  = zeros(this.NumberLinks,1);
      end
      %%%
      % Data center node price
      if isfield(slice_data, 'NodePrice')
        this.ServiceNodes.Price = slice_data.NodePrice;
      else
        this.ServiceNodes.Price = zeros(this.NumberServiceNodes,1);
      end
      
      % Flow Table
      % convert node index to virtual node index
      phy_nodes = this.PhysicalNodeMap;
      this.FlowTable = slice_data.FlowTable;
      for k = 1:height(this.FlowTable)
        path_list = this.FlowTable{k, 'Paths'};
        for p = 1:path_list.Width
          path_list{p}.node_list = phy_nodes(path_list{p}.node_list);
        end
      end
      
      this.VNFList = slice_data.VNFList;
      
      % Flow options
			this.options = setdefault(this.options, struct(...
				'SlicingMethod', SlicingMethod.AdjustPricing,...
				'PricingPolicy', 'linear'...
				));
      this.options = structupdate(this.options, slice_data, ...
				{'SlicingMethod', 'PricingPolicy', 'FlowPattern', 'DelayConstraint', 'NumberPaths'});
    end
    
    function delete(this)
			for f = 1:this.NumberFlows
				delete(this.FlowTable{f,'Paths'});
			end
    end
    
  end
  
  %% Property Get Methods
	methods
		function c = get.LinkCost(this)
			% the virtual links's unit cost
			c = this.Parent.readLink('UnitCost', this.Links.PhysicalLink);
		end
		
		function c = get.NodeCost(this)
			% the virtual data center nodes's unit cost
			c = this.Parent.readDataCenter('UnitCost', this.getDCPI());
    end
    
    function f = get.NumberFlows(this)
			f = height(this.FlowTable);
    end
    
		function p = get.NumberPaths(this)
			p = 0;
			for i=1:this.NumberFlows
				p = p + this.FlowTable{i, 'Paths'}.Width;
			end
    end	
    
		function n = get.NumberVNFs(this)
			n = length(this.VNFList);
    end
    
  end
	
  %% Public Methods
	methods		
    % function eventhandler(this, source, eventData) [Abstract]
    
    %% Reset the resource allocation state of a slice
    % see also <Finalize>, <isFinal>.
    function initialize(this, prices)
      if nargin >= 2
        this.Links.Price = prices.Link;
        this.ServiceNodes.Price = prices.Node;
      else
        this.ServiceNodes{:,{'Price'}} = 0;
        this.Links{:,{'Price'}} = 0;
      end
      this.ServiceNodes{:,{'Load', 'Capacity'}} = 0;
      this.Links{:,{'Load', 'Capacity'}} = 0;
      this.FlowTable{:, 'Rate'} = 0;
      this.Optimizer.clear();
      this.b_final = false;
    end
    
    %% finalize
    % Called by the resource optimization method in the network.
    function finalize(this, prices)
			assert(this.b_final==false, 'Slice state has not been initialized (consider re-initialize it after optimization).');
      if this.NumberFlows > 0
        % subclass may override _postProcessing_, it is recommended that
        % the subclass maintains the default behavior of this function.
        this.Optimizer.postProcessing();
        this.ServiceNodes.Load = this.getNodeLoad();
        this.Links.Load = this.getLinkLoad();
        this.ServiceNodes.Capacity = this.ServiceNodes.Load;
        this.Links.Capacity = this.Links.Load;
        this.FlowTable.Rate = this.getFlowRate;
        %this.setPathBandwidth;
      end
      
      if nargin >= 2
        this.ServiceNodes.Price = prices.Node(this.getDCPI);
        this.Links.Price = prices.Link(this.Links.PhysicalLink);
      end
      
      this.b_final = true;
    end
    
    function tf = isFinal(this)
      tf = this.b_final;
		end
    
		function ye = getLinkLoad(this, varargin)
			ye = this.op.getLinkLoad(varargin{:});
		end
		
		function v_n = getNodeLoad(this, varargin)
			v_n = this.op.getNodeLoad(varargin{:});
		end
		
    %% Capacity and Load
    % Public interface for network to inquire the resource occupation of the slices.
    % In class <Slice>, <getLinkCapacity> and <getNodeCapacity> are equal to the
    % protected methods <getLinkLoad> and <getNodeLoad> respectively. But in
    % subclasses of <Slice>, the slice load might be less than its capacity, so that
    % the two group of methods return different results.
		function c = getNodeCapacity(this, isfinal)
			if nargin == 1 
				isfinal = true;
			end
			c = this.op.getNodeLoad(isfinal);
		end

		function c = getLinkCapacity(this, isfinal)
			if nargin == 1 
				isfinal = true;
			end
			c = this.op.getLinkLoad(isfinal);
		end
		
		function r = getFlowRate(this, varargin)
			r = this.op.getFlowRate(varargin{:});
		end
		%%
    % 		function pid = getLocalPathId(this, path)
    % 			%%%
    % 			% {need override}
    % 			cprintf('comment', '%s%s%s', ...
    % 				'getLocalPathId should be overrided,', ...
    % 				'if subclasses dynamically manage flows. ', ...
    % 				'Instead, using path.local_id, which should be dynamically maintained');
    % 			pid = path.id - this.FlowTable.Paths(1).paths{1}.id + 1;
    % 		end
		
		%         function pl = getPathLength(this)
		%             pl = zeros(this.NumberPaths,1);
		%             pid = 1;
		%             for i = 1:this.NumberFlows
		%                 pathlist = this.FlowTable{i,'Paths'}.paths;
		%                 for l = 1:length(pathlist)
		%                     pl(pid) = pathlist{l}.Length;
		%                 end
		%             end
		%         end
  end
	
  %% Protected Methods
	methods (Access = protected)
		function newobj = copyElement(this)
      if this.isShallowCopyable
        newobj = copyElement@VirtualNetwork(this);
      else
        newobj = this;
      end
			% Deep Copy
			% [FlowTable.Paths]
			for f = 1:this.NumberFlows
				% 'Paths' is a handle object, is should be copyed to the new table.
				newobj.FlowTable{f,'Paths'} = this.FlowTable{f,'Paths'}.copy;
			end
		end
    
	end
	
	methods (Access = {?Slice, ?SliceOptimizer})
%% Resource Cost of a Slice
    % Both the node and link have cost linear with its load. Subclasses may override this
    % method to provide mode realistic cost model, such as including static cost or
    % non-linear cost model.
    %
    % *Static Cost*: In the scenario of resource sharing, both static cost model and how
    % to share the static cost between slices are not determined. One option for sharing
    % static cost is shareing by usage, i.e., if a slice use a resource, then it shares a
    % part of the static cost. Another option is sharing by load, i.e., the static cost
    % depends on the portion of used resources at the network. 
    %
    % When calculate network cost as a single slice, this method equals to
    % <PhysicalNetwork.totalCost>.
    function rc = getResourceCost(this, load)
      if nargin <= 1 
        load.Node = this.ServiceNodes.Capacity;
        load.Link = this.Links.Capacity;
      end
      
      pn = this.Parent; % A slice should be have a parent network (so that |this.Parent| is valid).
      link_uc = pn.getLinkCost(this.Links.PhysicalLink);
      node_uc = pn.getNodeCost(this.getDCPI);
      rc = dot(link_uc, load.Link) + dot(node_uc, load.Node);
		end
	end
	
  %% Static Methods
	methods(Static)
		%% TODO: move to Network and split it according to type
		slice_template = loadSliceTemplate(index);
	end

	methods
		[payment, grad, pseudo_hess] = fcnLinkPricing(this, link_price, link_load);
		[payment, grad, pseudo_hess] = fcnNodePricing(this, node_price, node_load);
	end
	
end

