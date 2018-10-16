%% Virtual Network
% define resource information in a virtual network.
classdef VirtualNetwork < matlab.mixin.Copyable & matlab.mixin.Heterogeneous
	% Specify the properties that can only be modified by Physcial Network directly
	properties
		Parent;
		Identifier;     %
		
		Topology;       % <DirectedGraph> topology information of the slice
		Links;					% information of virtual links in the slice
		Nodes;					% information of virtual nodes
		DataCenters;		% data centers in the slice
		options;
	end
	
	properties (Dependent)
		NumberNodes;					% Number of nodes in the slice
		NumberDataCenters;    % Number of virtual data centers.
		NumberLinks;					% Number of edges in the slice
		PhysicalLinkMap;			%
		PhysicalNodeMap;
	end
	
	methods
		function this = VirtualNetwork(vnet_data)
			if nargin == 0
				return;
			elseif isa(vnet_data, 'VirtualNetwork')
				this = vnet_data.copy;
				return;
			end
			this.Parent = vnet_data.Parent;		% 'Parent' is mandatory
			if isfield(vnet_data, 'Identifier')
				this.Identifier = vnet_data.Identifier;
			end
			
			this.Topology = DirectedGraph(vnet_data.Adjacent);
			% Virtual Links
			this.Links = array2table(vnet_data.LinkMapS2P, 'VariableNames', {'PhysicalLink'});
			this.Links.Capacity = zeros(height(this.Links),1);	% Link capacity
			this.Links.Load = zeros(height(this.Links), 1); 	  % Link load
			% Virtual Nodes
			this.Nodes = array2table(vnet_data.NodeMapS2P, 'VariableNames', {'PhysicalNode'});
			% Virtual Data Center Nodes
			% Select the data center nodes from all the virtual nodes of this slice.
			dc_vn_index = ...
				find(this.Parent.getNodeField('Capacity', this.Nodes.PhysicalNode) > 0);
			this.DataCenters = array2table(dc_vn_index, 'VariableNames', {'VirtualNode'});
			this.Nodes.DataCenter = zeros(this.NumberNodes,1);
			this.Nodes{dc_vn_index, 'DataCenter'} = (1:height(this.DataCenters))';
			this.DataCenters.Capacity = zeros(height(this.DataCenters),1);	% Data center node capacity
			this.DataCenters.Load = zeros(height(this.DataCenters), 1); 		% Data center node load
		end
		
		function delete(this)
			delete(this.Topology);
		end
		
	end
	
	methods (Access = protected)
		function newobj = copyElement(this)
			% Make a shallow copy of all properties
			newobj = copyElement@matlab.mixin.Copyable(this);
			% Deep Copy
			% [Topology]
			newobj.Topology = this.Topology.copy;
			% [Parent]: should be updated by caller of <copy>.
		end
	end
	
	methods
		function n = get.NumberNodes(this)
			n = height(this.Nodes);% n = this.Topology.NumberNodes;
		end
		
		function n = get.NumberDataCenters(this)
			n = height(this.VirtualDataCenters);
		end
		
		function m = get.NumberLinks(this)
			m = height(this.Links); % m = this.Topology.NumberEdges;
		end
		
		function map = get.PhysicalLinkMap(this)
			map = zeros(this.Parent.NumberLinks,1);
			map(this.Links.PhysicalLink) = 1:this.NumberLinks;
		end
		
		function map = get.PhysicalNodeMap(this)
			map = zeros(this.Parent.NumberNodes,1);
			map(this.Nodes.PhysicalNode) = 1:this.NumberNodes;
		end
	end
	
	methods
		%%%
		% Get physical node index of data center
		function dc_node_id = getDCNI(this, dc_index)
			if nargin == 1
				dc_index = 1:this.NumberDataCenters;
			end
			vn_id = this.DataCenters.VirtualNode(dc_index);
			dc_node_id = this.Nodes.PhysicalNode(vn_id);
		end
		
		%%%
		% Get the physical index of data center.
		function dc_phy_index =getDCPI(this,dc_index)
			if nargin == 1
				dc_node_id = this.getDCNI;
			else
				dc_node_id = this.getDCNI(dc_index);
			end
			dc_phy_index = this.Parent.getNodeField('DataCenter', dc_node_id);
		end
		
	end
	
	methods
		%% Resource Utilization
		% Calculate the average and overall resource utilization ratio. 
		% Also, calculate the resource utilization of links and nodes separately. 
		%
		% NOTE: Exclude the resource with no capacity.
		%
		% Subclasses might override this method to provide different measure of
		% resource utilization.
		function [theta, t_link, t_node] = utilizationRatio(this)
			capacities = [this.DataCenters.Capacity; this.Links.Capacity];
			load = [this.DataCenters.Load; this.Links.Load];
			idx = capacities>eps;
			if isempty(idx)
				error('error: this virtual network has no capacity.');
			end

			theta.Mean = mean(load(idx)./capacities(idx));
			theta.Overall = sum(load)/sum(capacities);
			if nargout >= 2
				eidx = this.Links.Capacity>eps;
				if isempty(eidx)
					error('error: no link capacity.');
				end
				t_link.Mean = mean(this.Links.Load(eidx)./this.Links.Capacity(eidx));
				t_link.Overall = sum(this.Links.Load(eidx))/sum(this.Links.Capacity(eidx));
			end
			if nargout >= 3
				nidx = this.DataCenters.Capacity>eps;
				if isempty(nidx)
					t_link = struct([]);
					warning('no node capacity.');
				end
				t_node.Mean = mean(this.Nodes.Load(nidx)./this.Nodes.Capacity(nidx));
				t_node.Overall = sum(this.Nodes.Load(nidx))/sum(this.Nodes.Capacity(nidx));
			end
		end
		
	end
end

%% Properties
% * *Topology*: Normally, the network slice will not run shrtest path algorithm, so the
% absolute value of the adjacent matrix of Topology does not matter. On the other hand,
% the link and node capacity of the slice is also not determined until the substrate
% network allocate the resource to the slice.
% * *Links* : fields include _PhysicalLink_, _Price_, _Load_.
% * *Nodes* : fields include _PhysicalNode_, _Price_, _Load_.
%
% * *NumberNodes* |get|
%
% * *NumberLinks* |get|
%
