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
%% Single-Function Model
% In the same slice, all flows request the same set of VNFs, and the processing demand is
% proportion to flow rate, so we can treat all VNF on the function chain as one VNF
% type,and its resource demand is the sum of all consisting VNFs.  
% Since slices distinguish in VNF chains, the assembled VNF's processing demand also
% diverses. Even that, the assembled VNFs of individual slices can also be treated as one
% type, although the processing demand is different, which can be addressed by setting
% coefficients for flows from different slices.   
%
% Post-processing is required to find an solution, we can allocate VNF resources by order.
%
% Note: The 'single-function' model can only be applied when all NFV-capable nodes can
% instantiate all types of VNF. Otherwise, the assembled VNF cannot be located at any
% NFV-capable nodes.
% 
% TODO: Part of the functionalities are moved to <ClouldNetworkEx>.
function [output, runtime] = singleSliceOptimization( this, new_opts )
% this.clearStates;
if nargin < 2
    new_opts = struct;
end
options = getstructfields(this.options, {'SlicingMethod', 'PricingFactor'});
assert(options.SlicingMethod.IsSingle,...
	'error[%s]: unrecognized method (%s).', calledby(0), options.SlicingMethod.char);
assert(isfield(options, 'PricingFactor'), ...
	'error[%s]: PricingFactor not specified, %s, %s',  calledby(0), ...
	'considerng provide it when creating the network',...
	'or specify it before calling this method.');
options = structmerge(options, getstructfields(new_opts, 'bCompact', 'default', true));
options.PricingPolicy = 'linear';       % can be specified by the input argument.

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

slice_data = this.updateSliceData(slice_data, options);      % override by subclasses
slice_data.NumberPaths = [];             % avoid warning, no use here.
slice_data.SlicingMethod = [];
slice_data.Parent = this;
% the flow id and path id has been allocated in each slice already, no need to reallocate.
ss = SimpleSlice(slice_data);
if options.SlicingMethod == SlicingMethod.SingleFunction
    ss.getAs_res(flow_owner, slice_data.Alpha_f);
    options.Alpha_f = slice_data.Alpha_f;
elseif options.SlicingMethod == SlicingMethod.SingleNormal
    %% Coefficient for global optimization
    % When all slices are combined into one slice, a VNF might not be used by all paths
    % (_i.e._ all flows). If a VNF |f| is not used by a path |p|, there is no
    % processing-rate constraints on $f \times p$(|delete_items|). To form the constraint
    % coefficient matrix, the related items in |As| should be removed. See also <Slice
    % file://E:/workspace/MATLAB/Projects/Documents/CloudNetworking/Slice.html>.
    %
    % By the way, $z_{n,p,f}=0, \forall n$, if |p| does not use NFV |f|.
    I_flow_function = zeros(NF, ss.NumberVNFs);
    for f = 1:NF
        [~, vid] = ismember(this.slices{flow_owner(f)}.VNFList, ss.VNFList);
        I_flow_function(f, vid) = 1;
    end
    I_path_function = ss.I_flow_path'*I_flow_function;
    ss.As_res = ss.As_res(logical(I_path_function(:)),:);
end

if nargout == 2
    tic;
end
% Only return intermediate results, so no return value provided.
ss.optimalFlowRate(options);  
if nargout == 2
    runtime.Serial = toc;
    runtime.Parallel = runtime.Serial;
end
if options.SlicingMethod == SlicingMethod.SingleNormal
    nz = ss.NumberDataCenters*ss.NumberPaths;
    z_index = 1:nz;
    for v = 1:ss.NumberVNFs
        mask_npf = ss.I_dc_path.*I_path_function(:,v)'; % compatible arithmetic operation
        ss.temp_vars.z(z_index) = mask_npf(:).*ss.temp_vars.z(z_index);
        z_index = z_index + nz;
    end
end
output = calculateOptimalOutput(this, ss);

%% Partition the network resources according to the global optimization
pid_offset = 0;
z_npf = reshape(full(ss.temp_vars.z), ss.NumberDataCenters, ss.NumberPaths, ss.NumberVNFs);
% node_load = zeros(this.NumberNodes, 1);
% link_load = zeros(this.NumberLinks, 1);
fmincon_opt = optimoptions('fmincon');
for s = 1:NS
    sl = this.slices{s};
    pid = 1:sl.NumberPaths;
    sl.temp_vars.x = ss.temp_vars.x(pid_offset+pid);
    nid = sl.getDCPI;       % here is the DC index, not the node index.
    [~, vid] = ismember(sl.VNFList, ss.VNFList);
    sl.temp_vars.z = ...
        reshape(z_npf(nid,pid+pid_offset,vid),sl.num_vars-sl.NumberPaths, 1);
    pid_offset = pid_offset + sl.NumberPaths;
    assert(sl.checkFeasible([sl.temp_vars.x; sl.temp_vars.z], ...
        struct('ConstraintTolerance', fmincon_opt.ConstraintTolerance)), 'error: infeasible solution.');
    sl.VirtualDataCenters.Capacity = sl.getNodeCapacity(false);
    sl.VirtualLinks.Capacity = sl.getLinkCapacity(false);
    % DEBUG
%     eid = sl.VirtualLinks.PhysicalLink;
%     node_load(nid) = node_load(nid) + sl.VirtualNodes.Capacity;
%     link_load(eid) = link_load(eid) + sl.VirtualLinks.Capacity;
end
% disp(max(node_load-this.getNodeField('Capacity')));
% disp(max(link_load-this.getLinkField('Capacity')));

%% Compute the real resource demand with given prices
options = structmerge(options, getstructfields(new_opts, 'PricingPolicy', 'ignore'));
if nargout == 2
    [prices, rt] = pricingFactorAdjustment(this, options);
    runtime.Serial = runtime.Serial + rt.Serial;
    runtime.Parallel = runtime.Parallel + rt.Parallel;
else
    [prices] = pricingFactorAdjustment(this, options);
end
% Finalize substrate network
this.finalize(prices);

%% Calculate the output
output.SingleSlice = ss;
options.Slices = this.slices;
output = this.calculateOutput(output, options);
end
