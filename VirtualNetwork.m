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
			
			this.Parent = vnet_data.Parent;		% 'Parent' is mandatory			
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
			this.ServiceNodes = array2table(dc_vn_index, 'VariableNames', {'VirtualNode'});
			this.Nodes.DataCenter = zeros(this.NumberNodes,1);
			this.Nodes{dc_vn_index, 'DataCenter'} = (1:height(this.ServiceNodes))';
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
		% Get physical node index of data center
		function dc_node_id = getDCNI(this, dc_index)
			if nargin == 1
				dc_index = 1:this.NumberDataCenters;
			end
			vn_id = this.ServiceNodes.VirtualNode(dc_index);
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
			capacities = [this.ServiceNodes.Capacity; this.Links.Capacity];
			load = [this.ServiceNodes.Load; this.Links.Load];
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
				nidx = this.ServiceNodes.Capacity>eps;
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
