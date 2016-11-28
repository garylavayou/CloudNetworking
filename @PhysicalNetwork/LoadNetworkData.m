%% Load network data
% load network data with speicifed network model, link options and node options.
%
%   graph_data = LoadNetworkData( model, link_opt, node_opt) 
% 
% |model|: specifies how to creat network data. See also NetworkModel.
%
% |link_opt|: link options.
%
% * *delay* |LinkDelayOption|: how to determine the links delay. See also LinkDelayOption.
%
% * *cost* |LinkCostOption|: how to decide the unit cost of the links. If
% |cost=LinkCostOption.NetworkSpecified|, then |link_opt| must contain a field of
% |link_cost|. See also LinkCostOption. 
%
% * *link_cost*: user specified link cost.
%
% * *delay2cost* |double|: when |link_opt.cost=LinkCostOption.LengthDependent|, this
% option specifies the coefficient to convert the delay to unit cost.
% 
% * *link_cost* |double|: unit cost of each link.
%
% |node_opt|: node options.
%
% * *capacity*: how to specified the network nodes' capacity data. See also
% NodeCapacityOption. 
%
% 
function graph_data = LoadNetworkData( model, link_opt, node_opt)
if nargin < 2 || ~isfield(link_opt, 'delay') || isempty(link_opt)
    link_opt.delay = LinkDelayOption.BandwidthPropotion;
end
if nargin < 2 || ~isfield(link_opt, 'cost') || isempty(link_opt)
    link_opt.cost = LinkCostOption.LengthDependent;
end
if nargin < 3 || ~isfield(node_opt, 'capacity') || isempty(node_opt)
    node_opt.capacity = NodeCapacityOption.Default;
end
if nargin < 3 || ~isfield(node_opt, 'cost') || isempty(node_opt.cost)
    node_opt.cost = NodeCostOption.Uniform;
end

switch model
    case NetworkModel.Sample1
        link_capacity =  1000*[      % Unit: Mbps
            0  6  8  0  0 13;
            7  0  9 14  0  0;
            9 10  0 10  0  7;
            0 15 11  0 13  0;
            0  0  0 14  0  8;
            14 0  8  0  9  0];
        node_name = {'Beijing'; 'Shanghai'; 'Wuhan'; 'Guangzhou'; 'Chengdu'; 'Xian'};
        location = [  40,116;    31,121;    30,114;   23,112;     30,104;   34,108];
        
        % node capacity
        switch node_opt.capacity
            case NodeCapacityOption.Default
                node_capacity = 1000*[50; 50; 100; 50; 40; 60];
            case NodeCapacityOption.UserSpecified
                node_capacity = node_opt.node_capacity;
            case NodeCapacityOption.BandwidthProportion
                node_capacity = PhysicalNetwork.AF0 * ...
                    (sum(link_capacity,2)+sum(link_capacity,1));
        end
        [tail, head, link_capacity]=find(link_capacity');
        
        % link delay
        if link_opt.delay==LinkDelayOption.BandwidthPropotion
            link_delay = PhysicalNetwork.LinkDelay(link_opt.delay, link_capacity);
            description = 'The link latency is proportial to the bandwidth, i.e. the longer path has larger bandiwdth';
        elseif link_opt.delay==LinkDelayOption.BandwidthInverse
            link_delay = PhysicalNetwork.LinkDelay(link_opt.delay, link_capacity);
            description = 'The link latency is inverse proportial to the bandwidth, i.e. the short path has larger bandiwdth';
        elseif link_opt.delay==LinkDelayOption.Constant
            link_delay = 10*ones(size(link_capacity));
            description = 'Each link has the same latency';
        elseif link_opt.delay==LinkDelayOption.Random
            rng(link_opt.RandomSeed);
            link_delay = 10*rand(size(link_capacity));
            description = 'Each link''s latency is random within a range';
        end
        
        graph_data = digraph(head, tail, link_delay, node_name);
        
        % Fields of Edge Table: EndNodes, Weight, Capacity, UnitCost
        graph_data.Edges.Properties.VariableDescriptions{2}=description;
        graph_data.Edges.Properties.UserData = {link_opt};
        graph_data.Edges.Capacity = link_capacity;
        switch link_opt.cost
            case LinkCostOption.Uniform
                link_opt.link_cost = ones(graph_data.numedges,1);
            case LinkCostOption.LengthDependent
                if ~isfield(link_opt, 'delay2cost') || isempty(link_opt.delay2cost)
                    link_opt.delay2cost = 0.1;
                end
                link_opt.link_cost = link_opt.delay2cost .* graph_data.Edges.Weight;
            case LinkCostOption.NetworkSpecified
                if ~isfield(link_opt, 'link_cost') || isempty(link_opt.link_cost)
                    error('Link cost information should be provided.');
                end
            otherwise
                error('invalid link cost option (%d).', link_opt.cost);
        end
        graph_data.Edges.UnitCost = link_opt.link_cost;
        graph_data.Edges.Properties.VariableUnits(2:3)={'ms','Gbps'};
        
        % Fields of Node Table: Name, Location, Capacity, StaticCost, UnitCost
        graph_data.Nodes.Location = location;
        graph_data.Nodes.Capacity = node_capacity;
        graph_data.Nodes.Properties.UserData = {node_opt};
        switch node_opt.cost
            case NodeCostOption.Uniform
                node_opt.node_cost = 3*ones(graph_data.numnodes,1)*mean(link_opt.link_cost);
            case NodeCostOption.NetworkSpecified
            otherwise
                error('invalid node cost option (%d).', node_opt.cost);
        end
        graph_data.Nodes.UnitCost = node_opt.node_cost;

    otherwise
        error('cannot process this network model.')
end

end
