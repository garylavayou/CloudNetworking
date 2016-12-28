%% Optimization in Single Slice 
% Optimize the resource allocation in a single slice.
%   All service flows are placed in one slice.
%   The VNFs in a slice are treated as one function.
%
% * *NOTE*: with the same optimal objective value, the two methods may result in different
% solutions.
% * *TODO*: to remove the unused VNFs from the single slice, check if b_vnf functions
% right with |I_path_function| and the output.
% * *TODO*: the solution of 'single-function' is infeasible.
% 
%%
function [output, ss] = singleSliceOptimization( this, options )
if nargin <= 1
    options.Display = 'final';
    options.Method = 'normal';
else
    if ~isfield(options, 'Display')
        options.Display = 'final';
    end
    if ~isfield(options, 'Method')
        options.Method = 'normal';
    end
end
this.clearStates;

%% Merge slices into one single big slice
NL = this.NumberLinks;
NN = this.NumberNodes;
NS = this.NumberSlices;
NP = this.NumberPaths;
slice_data.adjacent = this.graph.Adjacent;
slice_data.link_map_S2P = (1:NL)';
slice_data.link_map_P2S = (1:NL)';
slice_data.link_capacity = this.getLinkField('Capacity');
slice_data.node_map_S2P = (1:NN)';
slice_data.node_map_P2S = (1:NN)';
slice_data.node_capacity = this.getNodeField('Capacity');
slice_data.flow_table = table([],[],[],[],[],[],[],'VariableNames',...
    {this.slices{1}.FlowTable.Properties.VariableNames{:,:},'Weight'});
NF = this.NumberFlows;
flow_owner = zeros(NF, 1);
nf = 0;
% b_vnf = false(this.NumberVNFs, 1);
if strcmp(options.Method, 'single-function')
    slice_data.alpha_f = zeros(NS, 1);
end
for s = 1:NS
    sl = this.slices{s};
    new_table = sl.FlowTable;
    % Map the virtual nodes to physical nodes.
    new_table.Source = sl.VirtualNodes{new_table.Source, {'PhysicalNode'}};
    new_table.Target = sl.VirtualNodes{new_table.Target, {'PhysicalNode'}};
    for f = 1:height(sl.FlowTable)
        % path_list is handle object, is should be copyed to the new table.
        path_list = PathList(sl.FlowTable{f,'Paths'});
        for p = 1:path_list.Width
            path = path_list.paths{p};
            path.node_list = sl.VirtualNodes{path.node_list,{'PhysicalNode'}};
        end
        new_table{f,'Paths'} = path_list;
    end
    new_table.Weight = sl.weight*ones(height(new_table),1);
    slice_data.flow_table = [slice_data.flow_table; new_table];
    flow_owner(nf+(1:sl.NumberFlows)) = s;
    nf = nf + sl.NumberFlows;
    if strcmp(options.Method, 'normal')
%         b_vnf(this.slices{s}.VNFList) = true;
    elseif strcmp(options.Method, 'single-function')
        slice_data.alpha_f(s) = sum(this.VNFTable{sl.VNFList,{'ProcessEfficiency'}});
    end
end
if strcmp(options.Method, 'normal')
    slice_data.VNFList = 1:this.NumberVNFs;
%     slice_data.VNFList = find(b_vnf);
else
    slice_data.VNFList = 1;
end
slice_data.parent = this;
NV = length(slice_data.VNFList);        % NV is not the true number of VNFs when the method is 'single-function'
% the flow id and path id has been allocated in each slice already, no need to reallocate.
ss = Slice(slice_data);
I_flow_function = zeros(NF, this.NumberVNFs);
for f = 1:NF
    I_flow_function(f, this.slices{flow_owner(f)}.VNFList) = 1;
end
ss.I_path_function = ss.I_flow_path'*I_flow_function;

net_profit = ss.optimalFlowRate(options);

%% Finalize substrate network
% # Record the resource allocation variables, flow rate, virtual node/link load of each
% slice.   
% # Calculate and announce the resource prices to each slice.
% # Record the substrate network's node/link load, price.
link_uc = this.getLinkField('UnitCost');
node_uc = this.getNodeField('UnitCost');
epsilon = this.unitStaticNodeCost;
phis_n = epsilon*this.delta*(NN-1)/this.totalNodeCapacity;
phis_l = epsilon*(1-this.delta)*(NN-1)/this.totalLinkCapacity;
link_price = (link_uc + phis_l) * (1 + options.PricingFactor);
node_price = (node_uc + phis_n) * (1 + options.PricingFactor);
pid_offset = 0;
if strcmp(options.Method, 'single-function')
    z_npf = reshape(ss.Variables.z, NN, NP, this.NumberVNFs);
elseif strcmp(options.Method, 'normal')
    z_npf = reshape(ss.Variables.z, NN, NP, NV);
end
options.Tolerance = 10^-2;
for s = 1:NS
    sl = this.slices{s};
    pid = 1:sl.NumberPaths;
    sl.Variables.x = ss.Variables.x(pid_offset+pid);
    nid = sl.VirtualNodes.PhysicalNode;
    eid = sl.VirtualLinks.PhysicalLink;
    vid = sl.VNFList;
    sl.Variables.z = ...
        reshape(z_npf(nid,pid+pid_offset,vid),sl.num_vars-sl.NumberPaths, 1);
    pid_offset = pid_offset + sl.NumberPaths;
    if ~sl.checkFeasible([], options)
        error('error: infeasible solution.');
    end
    sl.VirtualNodes.Load = sl.getNodeLoad;
    sl.VirtualLinks.Load = sl.getLinkLoad;
    sl.FlowTable.Rate = sl.getFlowRate;
    sl.setPathBandwidth;
    sl.VirtualLinks.Price = link_price(eid);
    sl.VirtualNodes.Price = node_price(nid);
end
this.setLinkField('Load', ss.VirtualLinks.Load);
this.setNodeField('Load', ss.VirtualNodes.Load);
this.setLinkField('Price', link_price);
this.setNodeField('Price', node_price);

%% Calculate the output
% The output variables includes,
%
% # Price and Load of physical nodes and links;
% # The approximation of net social welfare, and the accurate net social welfare;
% # The net profit of each slice and the substrate network, including four results from
% different methods, i.e. |ApproximatePercent|, |ApproximatePrice|, |AccuratePercent|,
% |AccuratePrice|. 
% # Flow rate of all flows in the network.
output.link_load = ss.VirtualLinks.Load;
output.node_load = ss.VirtualNodes.Load;
output.link_price = link_price;
output.node_price = node_price;
output.welfare_approx = net_profit;
output.flow_rate = ss.FlowTable.Rate;
if ~strcmp(options.Display, 'off') && ~strcmp(options.Display, 'none')
    fprintf('\tThe optimal net social welfare of the network: %G.\n', ...
        output.welfare_approx);
end
options.Model = 'Accurate';
output.welfare_accurate = ...
    sum(ss.FlowTable.Weight.*log(ss.FlowTable.Rate)) ...
    - this.getNetworkCost([], [], options.Model);
t = zeros(this.NumberSlices+1, 1);
output.profit = table(t,t,t,t, 'VariableNames', ...
    {'ApproximatePercent', 'ApproximatePrice', 'AccuratePercent', 'AccuratePrice'});
clear t;
for s = 1:NS
    sl = this.slices{s};
    nid = sl.VirtualNodes.PhysicalNode;
    
    var_x = [sl.Variables.x;sl.Variables.z];
    p = -Slice.fcnNetProfit(var_x, sl);
    output.profit.ApproximatePercent(s) = options.PercentFactor * p;
    p = -Slice.fcnNetProfit(var_x, sl, options);
    idx = sl.VirtualNodes.Load>0;
    p = p - dot(sl.VirtualNodes.Load(idx)./ss.VirtualNodes.Load(nid(idx)),...
        this.getNodeField('StaticCost', nid(idx)));
    output.profit.AccuratePercent(s) = options.PercentFactor * p;
    
    output.profit.ApproximatePrice(s) = -Slice.fcnProfit(var_x, sl);
end
f = (1-options.PercentFactor)/options.PercentFactor;
output.profit.ApproximatePercent(end) = ...
    sum(output.profit.ApproximatePercent(1:(end-1))*f);
output.profit.AccuratePercent(end) = ...
    sum(output.profit.AccuratePercent(1:(end-1))*f);
output.profit.ApproximatePrice(end) = ...
    output.welfare_approx - sum(output.profit.ApproximatePrice(1:(end-1)));
output.profit.AccuratePrice = output.profit.ApproximatePrice;
output.profit.AccuratePrice(end) = ...
    output.welfare_accurate - sum(output.profit.AccuratePrice(1:(end-1)));
end
% delta_lambda.n = single_slice.getNodeLoad-this.getNodeField('Capacity')
% delta_lambda.e = single_slice.getLinkLoad-this.getLinkField('Capacity')
