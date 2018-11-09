%% Virtual Network
% define resource information in a virtual network.
classdef VirtualNetwork < INetwork
	% Specify the properties that can only be modified by Physcial Network directly
	properties
		Parent;
		ServiceNodes;		% data centers in the slice
	end
	
	properties (Dependent)
		NumberServiceNodes;   % Number of virtual data centers.
		PhysicalLinkMap;			%
		PhysicalNodeMap;
	end
	
	methods
		function this = VirtualNetwork(vnet_data)
			if nargin == 0
				args = cell(0);
			else
				args = {vnet_data};
			end
			this@INetwork(args{:});
			if nargin == 0
				return;
			end
			
			this.Parent = vnet_data.Parent;		% 'Parent' fields is mandatory.
			% Virtual Links
			this.Links = array2table(vnet_data.LinkMapS2P, 'VariableNames', {'PhysicalLink'});
			this.Links.Capacity = zeros(height(this.Links),1);	% Link capacity
			this.Links.Load = zeros(height(this.Links), 1); 	  % Link load
			% Virtual Nodes
			this.Nodes = array2table(vnet_data.NodeMapS2P, 'VariableNames', {'PhysicalNode'});
			% Virtual Data Center Nodes
			% Select the data center nodes from all the virtual nodes of this slice.
			%
			%		intersect(this.Parent.DataCenters.Node, this.Nodes.PhysicalNode)

			dc_vn_index = find(this.Parent.Nodes{this.Nodes.PhysicalNode, 'DataCenter'});
			% 			dc_vn_index = ...
			% 				find(this.Parent.readNode('Capacity', this.Nodes.PhysicalNode) > 0);
			this.ServiceNodes = array2table(dc_vn_index, 'VariableNames', {'VirtualNode'});
			this.Nodes.ServiceNode = zeros(this.NumberNodes,1);
			this.Nodes{dc_vn_index, 'ServiceNode'} = (1:height(this.ServiceNodes))';
			this.ServiceNodes.Capacity = zeros(height(this.ServiceNodes),1);	% Data center node capacity
			this.ServiceNodes.Load = zeros(height(this.ServiceNodes), 1); 		% Data center node load
		end
		
		function delete(this)
			delete(this.graph);
		end
		
	end
	
	methods (Access = protected)
		% function newobj = copyElement(this)
		% [Parent]: should be updated by caller of <copy>.
		function [cap, load] = readCapacityLoad(this)
			cap.node = this.ServiceNodes{:,'Capacity'};
			cap.link = this.Links{:,'Capacity'};
			load.node = this.ServiceNodes{:,'Load'};
			load.link = this.Links{:,'Load'};
		end
	end
	
	methods
		function n = get.NumberServiceNodes(this)
			n = height(this.ServiceNodes);
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
		% Get physical node index of service nodes
		function sn_id = getSNPI(this, sn_index)
			if nargin == 1
				sn_index = 1:this.NumberServiceNodes;
			end
			vn_id = this.ServiceNodes{sn_index, 'VirtualNode'};
			sn_id = this.Nodes{vn_id, 'PhysicalNode'};
		end
		
		%%%
		% Get the data center index of the service node.
		function dc_index =getDCPI(this,sn_index)
			if nargin == 1
				sn_id = this.getSNPI;
			else
				sn_id = this.getSNPI(sn_index);
			end
			dc_index = this.Parent.readNode('DataCenter', sn_id);
		end
		
	end
	
end
