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
function [output, prices, runtime] = singleSliceOptimization( this, new_opts )
% this.clearStates;
net = this.hn;
defaultopts = structmerge(...
	getstructfields(net.options, {'SlicingMethod', 'PricingPolicy', 'PricingFactor'}, 'error'), ...
	getstructfields(this.options, {'Form'}, 'error'));
if nargin < 2
    options = defaultopts;
else
	options = structmerge(defaultopts, new_opts);
end
assert(options.SlicingMethod.IsSingle,...
	'error[%s]: unrecognized method (%s).', calledby(0), options.SlicingMethod.char);
assert(isfield(options, 'PricingFactor'), ...
	'error[%s]: PricingFactor not specified, %s, %s',  calledby(0), ...
	'considerng provide it when creating the network',...
	'or specify it before calling this method.');

%% Merge slices into one single big slice
Nl = net.NumberLinks;
Nn = net.NumberNodes;
Ns = net.NumberSlices;
slice_data.Adjacent = net.graph.Adjacent;
slice_data.LinkMapS2P = (1:Nl)';
slice_data.LinkMapP2S = (1:Nl)';
slice_data.LinkCapacity = net.readLink('Capacity');
slice_data.NodeMapS2P = (1:Nn)';
slice_data.NodeMapP2S = (1:Nn)';
slice_data.NodeCapacity = net.readDataCenter('Capacity');
slice_data.FlowTable = table([],[],[],[],[],[],[],'VariableNames',...
    {net.slices{1}.FlowTable.Properties.VariableNames{:,:},'Weight'});
flow_owner = zeros(net.NumberFlows, 1);
nf = 0;
for s = 1:Ns
    sl = net.slices{s};
    new_table = sl.FlowTable;
    % Map the virtual nodes to physical nodes.
    new_table.Source = sl.Nodes{new_table.Source, {'PhysicalNode'}};
    new_table.Target = sl.Nodes{new_table.Target, {'PhysicalNode'}};
    for f = 1:height(sl.FlowTable)
        % path_list is handle object, is should be copyed to the new table.
        path_list = PathList(sl.FlowTable{f,'Paths'});
        for p = 1:path_list.Width
            path_list{p}.node_list = sl.Nodes{path_list{p}.node_list,{'PhysicalNode'}};
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
slice_data.Parent = net;
% the flow id and path id has been allocated in each slice already, no need to reallocate.
ss = SimpleSlice(slice_data);
op = ss.getOptimizer(slice_data);
slice_data.flow_owner = flow_owner;
pricing_policy = options.PricingPolicy;
options.PricingPolicy = 'linear';  % the first step use the cost as price, so the policy is linear
if nargout >= 3 
	runtime = op.optimalFlowRateSingleSlice(slice_data, options);
else
	op.optimalFlowRateSingleSlice(slice_data, options);
end
output = calculateOptimalOutput(this, ss);

%% Compute the real resource demand with given prices
options.PricingPolicy = pricing_policy;
if nargout == 3
    [prices, rt] = pricingFactorAdjustment(net, options);
    runtime.Serial = runtime.Serial + rt.Serial;
    runtime.Parallel = runtime.Parallel + rt.Parallel;
else
    [prices] = pricingFactorAdjustment(net, options);
end

%% Calculate the output
output.SingleSlice = ss;
output.options = options;
end
