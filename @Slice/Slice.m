classdef Slice < VirtualNetwork
	properties
		Type;           % Type from slice template
		weight;         % weight for the slice user's utility.
		
		FlowTable;      % flow information in the slice
		VNFList;        % List of virtual network functions in the slice
	end
	
	properties (Access = protected)
    op;
  end
  
  properties
    b_final logical = false;
  end
	
	properties (Dependent)
    LinkCost;
    NodeCost;
		NumberFlows;					% Number of flows in the slice
		NumberPaths;					% Number of paths in the slice
		NumberVNFs;						% Number of virtual network functions in the slice
    Optimizer;
  end
	
  methods (Abstract)
    profit = getProfit(slice, options);
  end
  
  methods (Abstract, Access=protected)
    sc = getCost(this, load, model);
    ye = getLinkLoad(this, ~);
    v_n = getNodeLoad(this, ~);
    c = getLinkCapacity(this, isfinal);
    c = getNodeCapacity(this, isfinal);
    vc = getVNFCapacity(this, ~);
    r = getRevenue(this);
    setPathBandwidth(this, ~)
  end
  
	methods (Static, Abstract)
    [profit, grad] = fcnProfit(vars, slice, options);		% Objective function and gradient
    [profit, grad] = fcnSocialWelfare(x_vars, S, options);
    hs = fcnHessian(var_x, ~, slice, options);
  end
	
  %% Constructor and Destructor
  methods
    function this = Slice(slice_data)
      if nargin == 0
        return;
      end
      this@VirtualNetwork(slice_data);
      
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
        this.ServiceNodes.Price = zeros(this.NumberSeviceNodes,1);
      end
      
      % Flow Table
      % convert node index to virtual node index
      phy_nodes = this.PhysicalNodeMap;
      this.FlowTable = slice_data.FlowTable;
      for k = 1:height(this.FlowTable)
        path_list = this.FlowTable{k,{'Paths'}};
        for p = 1:path_list.Width
          path = path_list.paths{p};
          path.node_list = phy_nodes(path.node_list);
        end
      end
      
      this.VNFList = vnet_data.VNFList;
      
      % Flow options
      this.options = getstructfields(vnet_data, ...
        {'FlowPattern', 'DelayConstraint', 'NumberPaths'});
      this.options = structmerge(this.options, ...
        getstructfields(slice_data, 'SlicingMethod', 'default', SlicingMethod.AdjustPricing),...
        getstructfields(slice_data, 'PricingPolicy', 'ignore'));
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
			c = this.Parent.readNode('UnitCost', this.getSNPI);
    end
    
    function f = get.NumberFlows(this)
			f = height(this.FlowTable);
    end
    
		function p = get.NumberPaths(this)
			p = 0;
			for i=1:this.NumberFlows
				p = p + this.FlowTable.Paths(i).Width;
			end
    end	
    
		function n = get.NumberVNFs(this)
			n = length(this.VNFList);
    end
    
    function op = get.Optimizer(this)
      op = this.op;
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
      this.FlowTable.Rate = 0;
      this.op.clear();
      this.b_final = false;
    end
    
    %% finalize
    % Called by the resource optimization method in the network.
    function finalize(this, prices)
      if this.NumberFlows > 0
        % subclass may override _postProcessing_, it is recommended that
        % the subclass maintains the default behavior of this function.
        this.postProcessing();
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
		function initializeState(this)
			Nsn = this.NumberServiceNodes;
			Np = this.NumberPaths;
			Nve = this.NumberLinks;
			Nf = this.NumberFlows;
			
			this.I_dc_path = sparse(Nsn, Np);
			this.I_edge_path = sparse(Nve, Np);
			this.I_flow_path = sparse(Nf, Np);
			this.path_owner = zeros(Np,1);
			% this.local_path_id = zeros(this.NumberPaths, 1);
			pid = 0;
			for fid=1:Nf
				path_list = this.FlowTable{fid,{'Paths'}};
				for j = 1:path_list.Width
					pid = pid + 1;
					this.I_flow_path(fid,pid) = 1;
					this.path_owner(pid) = fid;
					path = path_list.paths{j};
					path.local_id = pid;    % record the local path in the slice.
					for k = 1:(path.Length-1)
						e = path.Link(k);
						eid = this.Topology.IndexEdge(e(1),e(2));
						this.I_edge_path(eid, pid) = 1;
						dc_index = this.Nodes{e(1),'ServiceNode'};
						if dc_index~=0
							this.I_dc_path(dc_index, pid) = 1;
						end
					end
					dc_index = this.Nodes{e(2),'ServiceNode'}; % last node
					if dc_index~=0
						this.I_dc_path(dc_index, pid) = 1;
					end
				end
			end
		end
		[payment, grad, pseudo_hess] = fcnLinkPricing(this, link_price, link_load);
		[payment, grad, pseudo_hess] = fcnNodePricing(this, node_price, node_load);
	end
	
	methods (Static)
		slice_template = loadSliceTemplate(index);
	end
	
end

