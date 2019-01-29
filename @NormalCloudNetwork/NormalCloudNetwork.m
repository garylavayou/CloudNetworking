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
		augraph DirectedGraph;
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
			num_aug_links = 2*length(new_heads);
			ext_link_idx = this.NumberLinks+(1:num_aug_links);
			ext_node_idx = this.NumberNodes+(1:numel(new_heads));
			warning off
			this.ext_links{ext_link_idx, 'EndNodes'} ...
				= [[new_heads; new_tails], [new_tails; new_heads]];
			this.ext_nodes{ext_node_idx, 'Location'} = new_poss;
			warning on
			this.ext_links{ext_link_idx, 'Delay'} = eps*ones(num_aug_links,1);
			this.ext_links{ext_link_idx, 'Capacity'} = inf;
			this.ext_nodes{:, 'DataCenter'} = 0;
			this.ext_nodes{ext_node_idx, 'DataCenter'} = (1:height(this.DataCenters))';
			this.augraph = DirectedGraph(this.ext_links, this.ext_nodes);
		end
	end
	
	%% Public Methods
	methods
		varargout = optimizeResourcePriceDual(this, slices, options);
		varargout = optimizeResourcePriceDual2(this, slices, options);
		varargout = optimizeResourcePriceScaling(this, slices, options);
		varargout = optimizeResourcePriceScaling2(this, slices, options);
		varargout = staticSlicing(this, slices, options);
		varargout = fixResourcePricing(this, slices, options)
		varargout = pricingResourceDual(this, slices, options);
		varargout = optimizeResourcePrice3(this, slices);
		
		function op = getOptimizer(this, varargin)
			this.op = NormalNetworkOptimizer(this, varargin{:});
			op = this.op;
		end
    
		% plot(this, b_undirect, b_extend)
		function plot(this, varargin)
			plot@PhysicalNetwork(this, varargin{:})
			if nargin >= 3
				b_extend = varargin{3};
			else
				b_extend = false;
			end
			if b_extend
				b_undirect = varargin{2};
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
		
		% Compatible with <PhysicalNetwork.LinkId>
		function [argout_1, argout_2] = LinkId(this, argin_1, argin_2)
			switch nargin
				case 1
					argout_1 = this.ext_links.EndNodes(:,1);
					argout_2 = this.ext_links.EndNodes(:,2);
				case 2
					argout_1 = this.ext_links.EndNodes(argin_1,1);
					argout_2 = this.ext_links.EndNodes(argin_1,2);
				otherwise
					% Link index: link is indexed by [head, tail] pairs. The appended augmented
					% links dose not change the indices of original links.
					argout_1 = this.augraph.IndexEdge(argin_1,argin_2);
			end
		end
		
	end
	
	methods (Access = protected)
		[sp_profit, b_violate, output] = SolveSCPDA(this, slices, prices, options);
		[sp_profit, b_violate, output] = SolveSCPDD(this, slices, prices, options);		% Dual Decomposition
		[sp_profit, b_violate, output] = SolveSCPPP(this, slices, prices, options);
	
		function sl = createslice(this, slice_opt, varargin)
			this.slices(end+1) = NormalSlice(slice_opt);
			sl = this.slices(end);
			sl.getOptimizer(slice_opt);
		end
		
	end
end

