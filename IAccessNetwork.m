classdef IAccessNetwork < handle	
	properties (Abstract)
		Nodes;
		DataCenters;
	end
	properties (Abstract, Dependent)
		NumberNodes;
		NumberDataCenters;
	end
	properties (SetAccess = protected)
		base_stations;      % node index of base stations
	end
	properties(Dependent)
		NumberBaseStations;
	end
	
	%% Public Methods
	methods
		function this = IAccessNetwork(node_opt)
			this.base_stations = [node_opt.bs_index; node_opt.dc_index.Gateway]; % the node should be naturally sorted
			this.Nodes.BaseStation = zeros(this.NumberNodes,1);
			this.Nodes{this.base_stations, 'BaseStation'} = (1:length(this.base_stations))';
			% 			this.base_stations= table([net.index.NBS'; dc.gateway'], ...
			% 				'VariableName', {'Node'});
		end
	end

	%% Property Get Methods
	methods
		function n = get.NumberBaseStations(this)
			n = length(this.base_stations);
		end
	end
	
	methods
		function value = readNode(this, name, node_id)
			if nargin <= 2
				node_id = 1:this.NumberNodes;
			end
			switch name
				case 'BaseStation'
					value = this.Topology.Nodes{node_id, 'BaseStation'};
				otherwise
					value = double.empty();
			end
		end
	end
	
	%% Public Static Methods
	methods (Static)
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
			node_opt.Model = NetworkModel.SDRAN;
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
	
	%% Protected methods
	methods (Access = protected)
		function slice_opt = preAddingSlice(this, slice_opt)
			switch slice_opt.FlowPattern
				case FlowPattern.RandomInterBaseStation
					if ~slice_opt.DuplicateFlow
						slice_opt.NumberFlows = min(slice_opt.NumberFlows, this.NumberBaseStations*(this.NumberBaseStations-1));
					end
					if isfield(slice_opt, 'NodeSet')
						assert(isempty(setdiff(slice_opt.NodeSet, this.base_stations)));
						slice_opt.NodeSet = slice_opt.NodeSet;
					else
						slice_opt.NodeSet = this.base_stations;
					end
				case FlowPattern.RandomDataCenter2BaseStation
					if ~slice_opt.DuplicateFlow
						slice_opt.NumberFlows = min(slice_opt.NumberFlows, this.NumberBaseStations*this.NumberDataCenters);
					end
					if isfield(slice_opt, 'BSNodeSet')
						assert(isempty(setdiff(slice_opt.BSNodeSet, this.base_stations)));
					else
						slice_opt.BSNodeSet = this.base_stations;
					end
					if isfield(slice_opt, 'DCNodeSet')
						assert(isempty(setdiff(slice_opt.DCNodeSet, this.DataCenters.Node)));
					else
						slice_opt.DCNodeSet = this.DataCenters.Node;
					end
				otherwise
			end
		end
		
		function end_points = generateEndPoints(this, slice_opt) %#ok<INUSL>
      end_points = zeros(1,2);
			switch slice_opt.FlowPattern
				case FlowPattern.RandomInterBaseStation
					id = unique_randi(length(slice_opt.NodeSet), 2, 'stable');
					end_points(:) = slice_opt.NodeSet(id);
				case FlowPattern.RandomDataCenter2BaseStation
					id = randi(length(slice_opt.DCNodeSet), 1);
					end_points(1) = slice_opt.DCNodeSet(id);
					bs = slice_opt.BSNodeSet;
					bs(bs==end_points(1)) = [];     % avoid the same node as BS and DC.
					end_points(2) = bs(randi(length(bs),1));
				otherwise
					end_points = double.empty();
			end
		end
		%         [flow_table, phy_adjacent, flag] = generateFlowTable( this, graph, slice_opt )

	end
end

