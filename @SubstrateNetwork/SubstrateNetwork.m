%% Substrate Network
% The <SubstateNetwork> interface.
% The functionalities of <SubstateNetwork> partly rely on <CloudNetwork>.
%
% In <CloudNetwork>, VNF-capable nodes are co-located with forwarding
% nodes, while in <SubstateNetwork>, VNF-capable nodes are connected with
% the proximal forwarding node with high capacity low latency link. We
% assume that the VNF-capable nodes in <SubstrateNetwork> are located on 
% pseudo-forwarding nodes. The links between VNF-capable nodes and
% forwarding nodes are assumed to have infinite capacity, and the
% processing capacity of forwarding nodes is zero.
%
classdef (Abstract) SubstrateNetwork < handle
	
	properties (Access=protected)
		VNFCapableNodes;    % The information table of VNF-capable nodes;
		ForwardingNodes;
		% Augment the network topology with data center nodes connected to the
		% forwarding node;
		SubstrateTopology;
		substrate_graph;
	end
	
	properties (Access = private)
		init_gamma_k double;
		init_q_k double;
	end
	
	properties (Dependent)
		NumberLinksAugmented;
		NumberNodesAugmented;
		AugmentedLinks;
		%AugmentedNodes;
	end
	
	%% Preassumed Properties
	% 	DataCenters
	% 	Topology
	% 	NumberNodes
	% 	NumberDataCenters
	
	methods
		[output, runtime] = optimizeResourcePrice(this, slices);
		[output, runtime] = optimizeResourcePrice2(this, slices, options);
		[output, runtime] = optimizeResourcePriceScaling(this, slices, options);
	end
	
	methods (Access = protected)
		results = SolveSCP(this, slices, node_price, link_price, options);
		results = SolveSCPDD(this, slices, node_price, link_price, options);		% Dual Decomposition
		results = SolveSCPPP(this, slices, node_price, link_price, options);
		[sp_profit, b_violate, violates] = SolveSCPCC(this, slices, node_price, link_price, options)
	end
	
	methods
		%% NOTE: the augmented edges are not used in slices.
		% The only usage is to plot the augmented topolgy.
		% We only use the mapping between data centers and the augmented
		% nodes (co-located forwarding nodes).
		function this = SubstrateNetwork()
			%% Augment the toplogy
			% The orignial forwarding nodes and the links are keep unchanged.
			this.SubstrateTopology = this.Topology;
			% 'this.DataCenters.Node' maps to the original colocated node;
			% 'this.DataCenters.AugmentedNode' maps to the nodes in the augmented
			% graph. 
			this.DataCenters.AugmentedNode = this.NumberNodes + (1:this.NumberDataCenters)'; %#ok<MCNPR>
			new_heads = this.DataCenters.AugmentedNode;
			new_tails = this.DataCenters.Node;
			%% add position to data center nodes
			% TODO: find a better place to place the data center nodes.
			delta_pos = [0, 150];
			new_poss = zeros(this.NumberDataCenters,2);
			for i = 1:this.NumberDataCenters
				new_poss(i,:) = this.Topology.Nodes{new_tails(i), 'Location'} + delta_pos;
			end
			
			%% Properties of New Edges
			% Delay: 0 (default);
			% UnitCost: 0 (deafult);
			% Capacity: Infinity;
			% Weight: 0 (default);
			%
			% NOTE: addedge maintains the edge indexing order, i.e.,
			% according to the ascending order of the <head, tail> value.
			% This help us to indexing the edge when we creat the
			% <DirectedGraph> object. On the other hand, we must specify the
			% property value for new edges. Otherwise, it will needs extra
			% effort to locate those edges after the adding operation (find it
			% according to the node index, "G.findnode(G.Edges.EndNodes(:,1))>this.NumberNodes").
			%
			edge_props = this.Topology.Edges(1, 2:end);
			edge_props{1:2*length(new_heads),:} = 0;
			edge_props{:,'Capacity'} = inf;
			%% Properties of New Nodes
			% DataCenter Index: the index of Data Center that connect to the
			%			forwarding node should be changed; Now, the augmented data center
			%			node maps to the Data Center.
			% Cost: Data Center's cost (the cost of the colocated forwarding
			%       node is set to 0)
			this.SubstrateTopology = this.SubstrateTopology.addedge(...
				[new_heads; new_tails], [new_tails; new_heads], edge_props);
			this.SubstrateTopology.Nodes{:, 'DataCenter'} = 0;
			this.SubstrateTopology.Nodes{(this.NumberNodes+1):end, 'DataCenter'} = ...
				(1:height(this.DataCenters))';
			this.SubstrateTopology.Nodes{(this.NumberNodes+1):end, 'Location'} = new_poss;
			%%
			this.substrate_graph = DirectedGraph(this.SubstrateTopology);
			% Upate the edge indexing mapping, see <PhysicalNetwork>.
			[s,t] = this.SubstrateTopology.findedge;
			idx = this.substrate_graph.IndexEdge(s,t);
			this.SubstrateTopology.Edges.Index = idx;
			this.SubstrateTopology.Edges.Properties.VariableDescriptions{5} = ...
				'Index the edges by column in the adjacent matrix.';
		end
		
		function n = get.NumberLinksAugmented(this)
			n = height(this.SubstrateTopology.Edges);
		end
		function n = get.NumberNodesAugmented(this)
			n = height(this.SubstrateTopology.Nodes);
		end
		
		function idx = get.AugmentedLinks(this)
			end_nodes = this.SubstrateTopology.Edges.EndNodes;
			idx = (end_nodes(:,1)>this.NumberNodes) | (end_nodes(:,2)>this.NumberNodes);
			idx = find(idx);
		end
		
		function plot(this, b_undirect)
			hold on;
			fig = gcf;
			ax = fig.Children(1);
			aug_edge_idx = this.AugmentedLinks;
			aug_node_idx = unique(this.SubstrateTopology.Edges{aug_edge_idx, 'EndNodes'});
			aug_part = digraph(this.SubstrateTopology.Edges(aug_edge_idx, :));
			sub_aug_part = subgraph(aug_part, aug_node_idx);
			locations = this.SubstrateTopology.Nodes{aug_node_idx, 'Location'};
			nodelabels = cell(1,sub_aug_part.numnodes);
			num_new_nodes = this.NumberNodesAugmented - this.NumberNodes;
			nodelabels(1:(sub_aug_part.numnodes - num_new_nodes)) = {' '};
			for i = (sub_aug_part.numnodes-num_new_nodes+1):sub_aug_part.numnodes
				nodelabels{i} = ['DC', num2str(i-(sub_aug_part.numnodes - num_new_nodes))];
			end
			axes(ax);
			args = {'XData', locations(:,1), ...
				'YData', locations(:,2), ...
				'NodeLabel', nodelabels, ...
				'Marker', 's',...
				'MarkerSize', 7,...
				'LineStyle', ':', ...
				'LineWidth', 2, ...
				'NodeColor', PhysicalNetwork.NodeColor(2,:), ...
				'EdgeColor', PhysicalNetwork.EdgeColor(2,:)};
			if b_undirect
				g = graph(sub_aug_part.adjacency+sub_aug_part.adjacency');
				g.plot(args{:});
			else
				sub_aug_part.plot(args{:});
			end
			hold off;
		end
		

		
		
	end
	
	methods (Access = protected)
		function sl = createslice(this, slice_opt, varargin)
			this.slices{end+1} = FlowEdgeSlice(slice_opt);
			sl = this.slices{end};
		end
	end
end

