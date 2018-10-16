%% CloudAccessNetwork
% Incorporate access network into cloud network.
% NOTE: CloudAccessNetwork can be derived from <CloudNetwork> and <AccessNetwork> (TODO),
% which are both derived from <PhysicalNetwork>.
classdef CloudAccessNetwork < CloudNetwork
	properties (GetAccess = public, SetAccess = protected)
		base_stations;      % node index of base stations
	end
	properties(Dependent)
		NumberBaseStations;
	end
	
	methods
		function n = get.NumberBaseStations(this)
			n = length(this.base_stations);
		end
		
		function this = CloudAccessNetwork(node_opt, link_opt, VNF_opt, options)
			%% Initialize as SD-RAN
			[node_opt, link_opt] = CloudAccessNetwork.loadNetworkData(node_opt, link_opt);
			this@CloudNetwork(node_opt, link_opt, VNF_opt, options);
			this.base_stations = [node_opt.bs_index; node_opt.dc_index.Gateway]; % the node should be naturally sorted
			this.Topology.Nodes.BaseStation = zeros(this.NumberNodes,1);
			this.Topology.Nodes{this.base_stations, 'BaseStation'} = ...
				(1:length(this.base_stations))';
			%             this.BaseStations = table([net.index.NBS'; dc.gateway'], ...
			%                 'VariableName', {'Node'});
		end
	end
	
	methods
		function value = getNodeField(this, name, node_id)
			if strcmp(name, 'BaseStation')
				value = this.Topology.Nodes{node_id, 'BaseStation'};
			else
				if nargin == 2
					value = getNodeField@CloudNetwork(this, name);
				else
					value = getNodeField@CloudNetwork(this, name, node_id);
				end
			end
		end
	end
	
	methods(Access = protected)
		function info = updateDemandInfo(this, slice_opt)
			switch slice_opt.FlowPattern
				case FlowPattern.RandomInterBaseStation
					info.NumberFlows = this.NumberBaseStations*(this.NumberBaseStations-1);
					if isfield(slice_opt, 'NodeSet')
						info.NodeSet = slice_opt.NodeSet;
					else
						info.NodeSet = this.base_stations;
					end
					info.NumberNodes = length(info.NodeSet);
				case FlowPattern.RandomDataCenter2BaseStation
					info.NumberFlows = this.NumberBaseStations*this.NumberDataCenters;
					if isfield(slice_opt, 'BSNodeSet')
						info.BSNodeSet = slice_opt.BSNodeSet;
					else
						info.BSNodeSet = this.base_stations;
					end
					if isfield(slice_opt, 'DCNodeSet')
						info.DCNodeSet = slice_opt.DCNodeSet;
					else
						info.DCNodeSet = this.DataCenters.Node;
					end
					info.NumberDCs = length(info.DCNodeSet);
				otherwise
					info = updateDemandInfo@CloudNetwork(this, slice_opt);
			end
		end
		
		function end_points = generateEndPoints(this, info, slice_opt)
			switch slice_opt.FlowPattern
				case FlowPattern.RandomInterBaseStation
					id = unique_randi(info.NumberNodes, 2, 'stable');
					end_points = info.NodeSet(id);
				case FlowPattern.RandomDataCenter2BaseStation
					id = randi(info.NumberDCs, 1);
					end_points(1) = info.DCNodeSet(id);
					bs = info.BSNodeSet;
					bs(bs==end_points(1)) = [];     % avoid the same node as BS and DC.
					end_points(2) = bs(randi(length(bs),1));
				otherwise
					end_points = generateEndPoints@CloudNetwork(this, info, slice_opt);
			end
		end
		%         [flow_table, phy_adjacent, flag] = generateFlowTable( this, graph, slice_opt )
		%         function slice_opt = preAddingSlice(this, slice_opt)
		%             warning('TODO');
		%         end
	end
	
	methods (Static)
		%% Override
		function [node_opt, link_opt] = loadNetworkData(node_opt, link_opt)
			% node is indexed by the order: router, GBS NBS.
			net = TwoLayerNetwork('SD-RAN');
			%             [this.index.Router; this.index.GBS]
			node_opt.dc_index.Core = [2; 3; 6; 7; 9; 12];
			node_opt.dc_index.Gateway = net.index.GBS';
			node_opt.bs_index = net.index.NBS';
			%% Set link capacity
			link_opt.Capacity = zeros(net.Size);
			link_opt.Cost = zeros(net.Size);
			link_capacity_map = [1000, 1000, 300, 100];
			link_cost_map = [1, 2, 4, 6];
			link_opt.CostModel = LinkCostOption.NetworkSpecified;
			for eid = 1:length(net.link)
				link_type = net.link(eid, 3);
				if link_type <= 3
					link_capacity = link_capacity_map(1+link_type);
					link_cost = link_cost_map(1+link_type);
				else
					link_capacity = 0;
					link_cost = inf;
				end
				link_opt.Capacity(net.link(eid,1),net.link(eid,2)) = link_capacity;
				link_opt.Capacity(net.link(eid,2),net.link(eid,1)) = link_capacity;
				link_opt.Cost(net.link(eid,1),net.link(eid,2)) = link_cost;
				link_opt.Cost(net.link(eid,2),net.link(eid,1)) = link_cost;
			end
			if isfield(link_opt, 'CapacityFactor') % 'CapacityFactor' used in <CloudNetwork>
				warning('CapacityFactor (%.2f) changes the link capacity configuration.',...
					link_opt.CapacityFactor);
			end
			[~, ~, link_opt.Cost] = find(link_opt.Cost);
			%             if isfield(link_opt, 'CostUnit')
			%                 link_opt.Cost = link_opt.Cost * link_opt.CostUnit;
			%             end
			%% Set node capacity
			node_opt.Model = NetworkModel.SD_RAN;
			node_opt.Location = net.node(:,1:2);
			node_opt.Capacity = zeros(net.Size,1);
			node_opt.Capacity(node_opt.dc_index.Core) = 5000;
			node_opt.Capacity(node_opt.dc_index.Gateway) = 1000;
			if isfield(node_opt, 'CapacityFactor')
				warning('CapacityFactor (%.2f) changes the node capacity configuration.',...
					node_opt.CapacityFactor);
			end
			node_opt.CapacityModel = NodeCapacityOption.NetworkSpecified;
			node_opt.CostModel = NodeCostOption.NetworkSpecified;
			node_opt.Cost = inf * ones(net.Size,1);
			node_opt.Cost(node_opt.dc_index.Core) = 1;
			node_opt.Cost(node_opt.dc_index.Gateway) = 3;
			%             if isfield(node_opt, 'CostUnit')
			%                 node_opt.Cost = node_opt.Cost * node_opt.CostUnit;
			%             end
		end
	end
	
end
