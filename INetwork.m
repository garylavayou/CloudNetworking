%% Virtual Network
% Define resource information in a network.
classdef (Abstract) INetwork < BaseCopyable & HeteroObject
	properties
		Identifier uint64;
		
		graph DirectedGraph;					% <DirectedGraph> topology information
		Links = table;					% <table> information of links
		Nodes = table;					% <table> information of nodes
		options Dictionary;
	end
	%%%
	% |Nodes|: the fields in node table include _Name_, _Location_.
	%
	% |Edges|: the fields in the edge table include _EndNodes_, _Weight_, _Capacity_,
	% _Index_, _Load_, _Price_.
	
	properties (Dependent)
    Optimizer;
		NumberNodes;					% Number of nodes in the network
		NumberLinks;					% Number of edges in the network
	end
	
	properties (Access = protected)
    op;
	end
		
	methods (Abstract)
		op = getOptimizer(this, options);
	end
	
	methods
		%
		% netdata: <struct>
		function this = INetwork(netdata)
			if nargin == 0	% default constructor
				return;
			end
			
			if nargin >= 1
				this.graph = DirectedGraph(netdata);
			end
			if isfield(netdata, 'Identifier')
				this.Identifier = netdata.Identifier;
			end		
			
			this.options = Dictionary();
		end
		
		function delete(this)
			if ~isempty(this.graph) && isvalid(this.graph)
				delete(this.graph);
			end
		end
	end
	
	methods (Access = protected)
		function newobj = copyElement(this)
			% Make a shallow copy of all properties
			if this.isShallowCopyable
				newobj = copyElement@BaseCopyable(this);
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
	
	%% Property Get Methods
	methods
		function n = get.NumberNodes(this)
			n = height(this.Nodes);% n = this.graph.NumberNodes;
		end
		
		function m = get.NumberLinks(this)
			m = height(this.Links); % m = this.graph.NumberEdges;
		end
		
		function op = get.Optimizer(this)
      op = this.op;
    end	
	end
	
	%% Public Methods
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
			loads = [load.node; load.link];
			idx = capacities>eps;
			if isempty(idx)
				error('error:[%s] this network has no capacity.', calledby);
			end
			
			theta.Mean = mean(loads(idx)./capacities(idx));
			theta.Overall = sum(loads)/sum(capacities);
			if nargout >= 2
				eidx = cap.link>eps;
				if isempty(eidx)
					error('error: no link capacity.');
				end
				t_link.Mean = mean(load.link(eidx)./cap.link(eidx));
				t_link.Overall = sum(load.link(eidx))/sum(cap.link(eidx));
			end
			if nargout >= 3
				nidx = cap.node>eps;
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
		% Get the link and node's capacity and load.
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
