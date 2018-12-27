%% Normal Cloud Network
% VNF-capable nodes are connected with the proximal forwarding node with high capacity low
% latency link. We assume that the VNF-capable nodes are located on pseudo-forwarding
% nodes. The links between VNF-capable nodes and forwarding nodes are assumed to have
% infinite capacity, and the processing capacity of forwarding nodes is zero.
classdef NormalCloudNetwork < PhysicalNetwork
	properties (Access=protected)
		VNFCapableNodes;    % The information table of VNF-capable nodes;
		ForwardingNodes;
	end
	
	properties (SetAccess = protected)
		% Augment the network topology with data center nodes connected to the
		% forwarding node;
		ext_links;
		ext_nodes;
	end
	
	%% Constructor
	methods
		function this = NormalCloudNetwork(varargin)
			this@PhysicalNetwork(varargin{:});
			%% NOTE: the augmented edges are not used in slices.
			% The only usage is to plot the augmented topolgy.
			% We only use the mapping between data centers and the augmented
			% nodes (co-located forwarding nodes).
			%% Augment the toplogy
			% The orignial forwarding nodes and the links are keep unchanged.
			this.ext_links = this.Links;
			this.ext_nodes = this.Nodes;
			% 'this.DataCenters.Node' maps to the original colocated node;
			% 'this.DataCenters.AugmentedNode' maps to the nodes in the augmented
			% graph. 
			this.DataCenters.AugmentedNode = this.NumberNodes + (1:this.NumberDataCenters)';
			new_heads = this.DataCenters.AugmentedNode;
			new_tails = this.DataCenters.Node;
			%% add position to data center nodes
			% TODO: find a better place to place the data center nodes.
			delta_pos = [0, 150];
			new_poss = zeros(this.NumberDataCenters,2);
			for i = 1:this.NumberDataCenters
				new_poss(i,:) = this.Nodes{new_tails(i), 'Location'} + delta_pos;
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
			%% Properties of New Nodes
			% DataCenter Index: the index of Data Center that connect to the
			%			forwarding node should be changed; Now, the augmented data center
			%			node maps to the Data Center.
			% Cost: Data Center's cost (the cost of the colocated forwarding
			%       node is set to 0)
			warning off
			this.ext_links{this.NumberLinks+(1:2*length(new_heads)), 'EndNodes'} ...
				= [[new_heads; new_tails], [new_tails; new_heads]];
			this.ext_nodes{this.NumberNodes+(1:numel(new_heads)), 'Location'} = new_poss;
			warning on
			this.ext_links{(this.NumberLinks+1):end, 'Capacity'} = inf;
			this.ext_nodes{:, 'DataCenter'} = 0;
			this.ext_nodes{(this.NumberNodes+1):end, 'DataCenter'} = ...
				(1:height(this.DataCenters))';
		end
	end
	
	%% Public Methods
	methods
		[output, runtime] = optimizeResourcePrice1(this, slices);
		[output, runtime] = optimizeResourcePrice2(this, slices, options);
		[output, runtime] = optimizeResourcePriceScaling(this, slices, options);

		function op = getOptimizer(this, options)
			if nargin == 1
				this.op = NormalNetworkOptimizer(this);
			else
				this.op = NormalNetworkOptimizer(this, options);
			end
			op = this.op;
		end
		
		% [output, runtime] = optimizeResourcePrice(this, sub_slices, options);

		function plot(this, b_undirect)
			hold on;
			fig = gcf;
			ax = fig.Children(1);
			aug_edge_idx = (this.NumberLinks+1):height(this.ext_links);
			aug_node_idx = (this.NumberNodes+1):height(this.ext_nodes);
			aug_part = digraph(this.ext_links(aug_edge_idx, :));
			sub_aug_part = subgraph(aug_part, aug_node_idx);
			locations = this.ext_nodes{aug_node_idx, 'Location'};
			nodelabels = cell(1,sub_aug_part.numnodes);
			num_new_nodes = height(this.ext_nodes) - this.NumberNodes;
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
		results = SolveSCP(this, slices, node_price, link_price, options);
		results = SolveSCPDD(this, slices, node_price, link_price, options);		% Dual Decomposition
		results = SolveSCPPP(this, slices, node_price, link_price, options);
	
		function sl = createslice(this, slice_opt, varargin)
			this.slices(end+1) = NormalSlice(slice_opt);
			sl = this.slices(end);
			sl.getOptimizer(slice_opt);
		end
		
	end
end

