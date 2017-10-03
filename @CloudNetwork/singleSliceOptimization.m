%% Optimization in Single Slice 
% Optimize the resource allocation in a single slice.
%   All service flows are placed in one slice.
%   The VNFs in a slice are treated as one function.
%
% * *NOTE*: with the same optimal objective value, the two methods may result in different
% solutions.
% * *TODO*: to remove the unused VNFs from the single slice, check if b_vnf functions
% right with |I_path_function| and the output.
% * *TODO*: the solution of |'single-function'| is infeasible.
% 
%%
function [output, runtime] = singleSliceOptimization( this )
% this.clearStates;

%% Merge slices into one single big slice
NL = this.NumberLinks;
NN = this.NumberNodes;
NS = this.NumberSlices;
slice_data.Adjacent = this.graph.Adjacent;
slice_data.LinkMapS2P = (1:NL)';
slice_data.LinkMapP2S = (1:NL)';
slice_data.LinkCapacity = this.getLinkField('Capacity');
slice_data.NodeMapS2P = (1:NN)';
slice_data.NodeMapP2S = (1:NN)';
slice_data.NodeCapacity = this.getDataCenterField('Capacity');
slice_data.FlowTable = table([],[],[],[],[],[],[],'VariableNames',...
    {this.slices{1}.FlowTable.Properties.VariableNames{:,:},'Weight'});
NF = this.NumberFlows;
flow_owner = zeros(NF, 1);
nf = 0;
% b_vnf = false(this.NumberVNFs, 1);
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
    slice_data.FlowTable = [slice_data.FlowTable; new_table];
    flow_owner(nf+(1:sl.NumberFlows)) = s;
    nf = nf + sl.NumberFlows;
end
slice_data.FlowPattern = FlowPattern.Default;
slice_data.DelayConstraint = inf;

slice_data = this.updateSliceData(slice_data);      % override by subclasses
slice_data.Parent = this;
% the flow id and path id has been allocated in each slice already, no need to reallocate.
ss = Slice(slice_data);
I_flow_function = zeros(NF, this.NumberVNFs);
for f = 1:NF
    I_flow_function(f, this.slices{flow_owner(f)}.VNFList) = 1;
end
ss.I_path_function = ss.I_flow_path'*I_flow_function;

if nargout == 2
    tic;
end
options = getstructfields(this.options, 'Method');
[~] = ss.optimalFlowRate(options);
if nargout == 2
    runtime.Serial = toc;
    runtime.Parallel = runtime.Serial;
end
output = calculateOptimalOutput(this, ss, slice_data);

%% Partition the network resources according to the global optimization
pid_offset = 0;
z_npf = output.Znpf;
output = rmfield(output, 'Znpf');
% node_load = zeros(this.NumberNodes, 1);
% link_load = zeros(this.NumberLinks, 1);
for s = 1:NS
    sl = this.slices{s};
    pid = 1:sl.NumberPaths;
    sl.x_path = ss.Variables.x(pid_offset+pid);
    nid = sl.getDCPI;       % here is the DC index, not the node index.
    vid = sl.VNFList;
    sl.z_npf = ...
        reshape(z_npf(nid,pid+pid_offset,vid),sl.num_vars-sl.NumberPaths, 1);
    pid_offset = pid_offset + sl.NumberPaths;
    if ~sl.checkFeasible([sl.x_path; sl.z_npf])
        error('error: infeasible solution.');
    end
    sl.VirtualDataCenters.Capacity = sl.getNodeLoad(sl.z_npf);
    sl.VirtualLinks.Capacity = sl.getLinkLoad(sl.x_path);
    % DEBUG
%     eid = sl.VirtualLinks.PhysicalLink;
%     node_load(nid) = node_load(nid) + sl.VirtualNodes.Capacity;
%     link_load(eid) = link_load(eid) + sl.VirtualLinks.Capacity;
end
% disp(max(node_load-this.getNodeField('Capacity')));
% disp(max(link_load-this.getLinkField('Capacity')));

%% Compute the real resource demand with given prices
if nargout == 2
    [node_price, link_price, rt] = pricingFactorAdjustment(this);
    runtime.Serial = runtime.Serial + rt.Serial;
    runtime.Parallel = runtime.Parallel + rt.Parallel;
else
    [node_price, link_price] = pricingFactorAdjustment(this);
end
% Finalize substrate network
this.finalize(node_price, link_price);

%% Calculate the output
output.SingleSlice = ss;
output = this.calculateOutput(output);
end
