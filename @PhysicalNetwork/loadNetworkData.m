%% Load network data
% load network data with speicifed network model, link options and node options.
%
%   graph_data = loadNetworkData( model, link_opt, node_opt)
%
% |model|: specifies how to creat network data. See also NetworkModel.
%
% |link_opt|: link options.
%
% * *delay* |LinkDelayOption|: how to determine the links delay. See also LinkDelayOption.
% * *cost* |LinkCostOption|: how to decide the unit cost of the links. If
% |cost=LinkCostOption.NetworkSpecified|, then |link_opt| must contain a field of
% |link_cost|. See also LinkCostOption.
% * *link_cost*: user specified link cost.
% * *delay2cost* |double|: when |link_opt.cost = LinkCostOption.LengthDependent|, this
% option specifies the coefficient to convert the delay to unit cost.
% * *link_cost* |double|: unit cost of each link.
%
% |node_opt|: node options.
%
% * *CapacityOption* |NodeCapacityOption|: how to specify the network nodes' capacity
% data. See also <NodeCapacityOption>.
% * *CapacityFactor* |double|: this parameter is used to map the link capacity to node
% capacity. Since we assume that the flow processing load is proportion to the flow rate,
% the node capacity is proportion to the adjacent link capacity. A unit of flow pass
% thorough a node will consume a unit of processing load, and two unit of link bandwidth
% (on the flow-in and flowout link), which is taken into consideration when mapping the
% node capacity.
% * *CostOption* |NodeCostOption|: how to specify the network node's unit cost. See also
% <NodeCostOption>.
% * *alpha* |double|: this parameter is used to specified the ratio of unit node cost to
% average unit link cost. With different |CostOption|, this parameter has different
% effect. When |CostOption=CapacityInverse|, this parameter is approximately equal to the
% ratio of average node unit cost to average link unit cost.
% * *cost_weight* |double|: different cost mapping factor that refects how the node
% resource is the scarce. 
function graph_data = loadNetworkData(node_opt, link_opt)
if nargin < 2 || ~isfield(link_opt, 'delay')
    warning('lacking link delay option, set to Random.'); 
    link_opt.delay = LinkDelayOption.Random;
end
if nargin < 2 || ~isfield(link_opt, 'cost')
    error('error: lacking link cost option.'); % link_opt.cost = LinkCostOption.LengthDependent;
end
if nargin < 1 || ~isfield(node_opt, 'capacity')
    error('error: lacking node capacity option.'); % node_opt.capacity = NodeCapacityOption.Default;
end
if nargin < 1 || ~isfield(node_opt, 'cost') 
    error('error: lacking node cost option.'); % node_opt.cost = NodeCostOption.Uniform;
end

switch node_opt.model
    case NetworkModel.Sample1
        capacity =  [      % Unit: Mbps
            0  6  8  0  0 13;
            7  0  9 14  0  0;
            9 10  0 10  0  7;
            0 15 11  0 13  0;
            0  0  0 14  0  8;
            14 0  8  0  9  0];
        node_name = {'Beijing'; 'Shanghai'; 'Wuhan'; 'Guangzhou'; 'Chengdu'; 'Xian'};
        location = [  116,40;    121,31;    114,30;   112,23;     104,30;   108,34];
        %% link information
        % For the convenience of assigning link capacity and other properties to the edge
        % table, we find the edges row-by-row.
        [tail, head, link_capacity]=find(capacity');
    case NetworkModel.Sample2
        capacity = [   % source: traffic engineering in software defined network
            0 1 1 1 0 0 0 0 0 0 0 0 0 0 0;    %1
            1 0 1 0 1 0 0 0 0 0 1 0 0 0 0;    %2
            1 1 0 1 0 1 1 0 0 0 0 0 0 0 0;    %3
            1 0 1 0 0 0 0 1 1 0 0 0 0 0 0;    %4
            0 1 0 0 0 0 0 0 0 0 0 1 0 0 0;    %5
            0 0 1 0 0 0 1 0 0 1 1 0 0 0 0;    %6
            0 0 1 0 0 1 0 0 1 1 0 0 0 0 0;    %7
            0 0 0 1 0 0 0 0 1 0 0 0 0 0 0;    %8
            0 0 0 1 0 0 1 1 0 1 0 0 0 0 1;    %9
            0 0 0 0 0 1 1 0 1 0 1 0 1 1 0;    %10
            0 1 0 0 0 1 0 0 0 1 0 1 1 0 0;    %11
            0 0 0 0 1 0 0 0 0 0 1 0 1 0 0;    %12
            0 0 0 0 0 0 0 0 0 1 1 1 0 1 0;    %13
            0 0 0 0 0 0 0 0 0 1 0 0 1 0 1;    %14
            0 0 0 0 0 0 0 0 1 0 0 0 0 1 0];   %15
        location = [96,-211; 221,-76; 173,-178; 138,-290; 327,-51; 278,-147; 234,-250;
            211,-341; 315,-327; 348,-256; 359,-119; 474,-84; 425,-184; 468,-259; 422,-327]; 
        [tail, head, link_capacity]=find(capacity');
    case NetworkModel.SD_RAN
        capacity = link_opt.Capacity;
        [tail, head, link_capacity]=find(capacity');
        location = node_opt.location;
    otherwise
        error('error: cannot process this network model.')
end
if isfield(link_opt, 'CapacityFactor')
    link_capacity = link_capacity * link_opt.CapacityFactor;
    capacity = capacity * link_opt.CapacityFactor;
end

%% link delay
switch link_opt.delay
    case LinkDelayOption.BandwidthPropotion
        link_delay = PhysicalNetwork.LinkDelay(link_opt.delay, link_capacity);
        description = 'The link latency is proportion to the bandwidth, i.e. the longer path has larger bandiwdth';
    case LinkDelayOption.BandwidthInverse
        link_delay = PhysicalNetwork.LinkDelay(link_opt.delay, link_capacity);
        description = 'The link latency is inverse proportion to the bandwidth, i.e. the short path has larger bandiwdth';
    case LinkDelayOption.Constant
        link_delay = 10*ones(size(link_capacity));
        description = 'Each link has the same latency';
    otherwise   % LinkDelayOption.Random
        if isfield(link_opt, 'RandomSeed')
            rng(link_opt.RandomSeed);
        else
            warning('random number seed is not specified (set as %d).', floor(now));
            rng(floor(now));
        end
        link_delay = 10*rand(size(link_capacity));
        description = 'Each link''s latency is random within a range';
end

if exist('node_name','var')
    graph_data = digraph(head, tail, link_delay, node_name);
else
    graph_data = digraph(head, tail, link_delay);
end
graph_data.Edges.Properties.VariableDescriptions{2}=description;

switch link_opt.cost
    case LinkCostOption.Uniform
        link_opt.link_cost = ones(graph_data.numedges,1);
    case LinkCostOption.LengthDependent
        if ~isfield(link_opt, 'delay2cost') || isempty(link_opt.delay2cost)
            link_opt.delay2cost = 1;
        end
        % the _Weight_ is equal to Latency (Length) of the link.
        link_opt.link_cost = link_opt.delay2cost .* graph_data.Edges.Weight;
    case LinkCostOption.NetworkSpecified
        if ~isfield(link_opt, 'link_cost') || isempty(link_opt.link_cost)
            error('Link cost information should be provided.');
        end
    case LinkCostOption.CapacityInverse
        link_opt.link_cost = 1./link_capacity;
    otherwise
        error('invalid link cost option (%d).', link_opt.cost);
end
if ~isfield(link_opt, 'CostUnit') || isempty(link_opt.CostUnit)
    link_opt.CostUnit = 1;
end
link_opt.link_cost = link_opt.CostUnit * link_opt.link_cost;

%% node capacity and cost
switch node_opt.capacity
    case NodeCapacityOption.NetworkSpecified
        node_capacity = reshape(node_opt.node_capacity,length(node_opt.node_capacity),1);
    case NodeCapacityOption.BandwidthProportion
        node_capacity = 1/2 * (sum(capacity,2)+sum(capacity,1)');
    otherwise
        error('error: Invalid option of NodeCapacityOption.');
end
if isfield(node_opt, 'CapacityFactor') && ~isempty(node_opt.CapacityFactor)
    node_capacity = node_capacity * node_opt.CapacityFactor;
end

switch node_opt.cost
    case NodeCostOption.Uniform
        if ~isfield(node_opt, 'alpha') || isempty(node_opt.alpha)
            node_opt.alpha = 1;
        end
        node_opt.node_cost = node_opt.alpha*mean(link_opt.link_cost)...
            *ones(graph_data.numnodes,1);
    case NodeCostOption.Weighted
        if ~isfield(node_opt, 'alpha') || isempty(node_opt.alpha)
            node_opt.alpha = 1;
        end
        if ~isfield(node_opt, 'cost_weight') || isempty(node_opt.cost_weight)
            node_opt.cost_weight = ones(graph_data.numnodes,1);
        end
        % NOTE: the node cost may be mapped according to the adjacent link's cost.
        % Here we use the average link cost to perform the cost mapping.
        node_opt.node_cost = node_opt.alpha*mean(link_opt.link_cost)...
            *(ones(graph_data.numnodes,1).*node_opt.cost_weight);
    case NodeCostOption.CapacityInverse
        node_opt.node_cost = 1 ./ node_capacity;
    case NodeCostOption.NetworkSpecified
        if ~isfield(node_opt, 'node_cost') || isempty(node_opt.node_cost)
            error('Node cost information should be provided.');
        end
    case NodeCostOption.Random
        error('TODO');
    otherwise
        error('error: Invalid node cost option (%d).', node_opt.cost);
end
if ~isfield(node_opt, 'CostUnit') || isempty(node_opt.CostUnit)
    node_opt.CostUnit = 1;
end
node_opt.node_cost = node_opt.CostUnit * node_opt.node_cost;

% Fields of Edge Table: EndNodes, Weight, Capacity, UnitCost
graph_data.Edges.Capacity = link_capacity;
graph_data.Edges.UnitCost = link_opt.link_cost;
graph_data.Edges.Properties.VariableUnits(2:3)={'ms','Gbps'};
graph_data.Edges.Properties.UserData = {link_opt};
% Fields of Node Table: Name, Location, Capacity, StaticCost, UnitCost
if exist('location','var')
    graph_data.Nodes.Location = location;
end
graph_data.Nodes.Capacity = node_capacity;
graph_data.Nodes.UnitCost = node_opt.node_cost;
graph_data.Nodes.Properties.UserData = {node_opt};
end
