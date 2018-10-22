%% Virtual Network
% Define resource information in a network.
classdef (Abstract) INetwork < matlab.mixin.Copyable & matlab.mixin.Heterogeneous & EventReceiver
	properties
		Identifier uint64;
		
		graph;					% <DirectedGraph> topology information
		Links;					% <table> information of links
		Nodes;					% <table> information of nodes
		options;
	end
	
	properties (Dependent)
		NumberNodes;					% Number of nodes in the network
		NumberLinks;					% Number of edges in the network
	end
	
	methods
		%
		% netdata: <struct>
		function this = Network(netdata)
			if nargin == 0	% default constructor
				return;
			end
			
			if nargin >= 1
				this.graph = DirectedGraph(netdata);
			end
			if isfied(netdata, 'Identifier')
				this.Identifier = netdata.Identifier;
			end
		end
		
	end
	
	methods (Access = protected)
		function newobj = copyElement(this)
			% Make a shallow copy of all properties
			if this.isShallowCopyable
				newobj = copyElement@matlab.mixin.Copyable(this);
			else
				newobj = this;
			end
			% Deep Copy
			% [graph]
			if newobj.graph == this.graph
				newobj.graph = this.graph.copy();
			end
		end
	end
	
	methods
		function n = get.NumberNodes(this)
			n = height(this.Nodes);% n = this.graph.NumberNodes;
		end
		
		function m = get.NumberLinks(this)
			m = height(this.Links); % m = this.graph.NumberEdges;
		end
	end
	
	methods

		function value = readNode(this, name, node_id)
			value = this.Nodes{node_id, {name}};
		end
		
		function value = readLink(this, name, link_id)
			value = this.Links{link_id, {name}};
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
			[cap,load] = this.readCapacityLoad();
			capacities = [cap.node; cap.link];
			load = [load.node; load.link];
			idx = capacities>eps;
			if isempty(idx)
				error('error:[%s] this network has no capacity.', calledby);
			end
			
			theta.Mean = mean(load(idx)./capacities(idx));
			theta.Overall = sum(load)/sum(capacities);
			if nargout >= 2
				eidx = cap.link>eps;
				if isempty(eidx)
					error('error: no link capacity.');
				end
				t_link.Mean = mean(load.link(eidx)./cap.link(eidx));
				t_link.Overall = sum(load.link(eidx))/sum(cap.link(eidx));
			end
			if nargout >= 3
				nidx = caop.node>eps;
				if isempty(nidx)
					t_link = struct([]);
					warning('no node capacity.');
				end
				t_node.Mean = mean(load.node(nidx)./cap.node(nidx));
				t_node.Overall = sum(load.node(nidx))/sum(cap.node(nidx));
			end
		end
	end
	
	methods (Abstract,Access = protected)
		[cap, load] = readCapacityLoad(this);
	end	
	
end

%% Properties
% * *graph*: Normally, the network slice will not run shrtest path algorithm, so the
% absolute value of the adjacent matrix of graph does not matter. On the other hand,
% the link and node capacity of the slice is also not determined until the substrate
% network allocate the resource to the slice.
% * *Links* : fields include _PhysicalLink_, _Price_, _Load_.
% * *Nodes* : fields include _PhysicalNode_, _Price_, _Load_.
%
% * *NumberNodes* |get|
%
% * *NumberLinks* |get|
%
