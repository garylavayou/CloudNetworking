%% Load Slice Template
%loadSliceTemplate Slice Template Type used to generate slice request.
%   slice_template = loadSliceTemplate(index)

%% Slice Template
% 
function [ slice_template ] = loadSliceTemplate(type_index)
if nargin == 0
    type_index = 1:100;
end
%% T1.1
i = 1;
if find(type_index==1,1)
    slice_template(i).Type = 1;
    slice_template(i).Pattern = FlowPattern.RandomMultiFlow;
    slice_template(i).Weight = 8;
    slice_template(i).NumberPaths = 3;
    slice_template(i).NumberFlows = 30;
    slice_template(i).NumberVNFs = 2;
    slice_template(i).ServiceInterval = 100;
    slice_template(i).Probability = 0.2;
    i = i + 1;
end
% slcie_opt.DelayConstraint = [];
% slice_opt.VNFList = [];

%% T1.2
if find(type_index==2,1)
    slice_template(i).Type = 2;
    slice_template(i).Pattern = FlowPattern.RandomMultiFlow;
    slice_template(i).Weight = 30;
    slice_template(i).NumberPaths = 3;
    slice_template(i).NumberFlows = 10;
    slice_template(i).NumberVNFs = 3;
    slice_template(i).ServiceInterval = 30;
    slice_template(i).Probability = 0.3;
    i = i + 1;
end

%% T1.3
if find(type_index==3,1)
    slice_template(i).Type = 3;
    slice_template(i).Pattern = FlowPattern.RandomMultiFlow;
    slice_template(i).Weight = 400;
    slice_template(i).NumberPaths = 1;
    slice_template(i).NumberFlows = 3;
    slice_template(i).NumberVNFs = 5;
    slice_template(i).ServiceInterval = 10;
    slice_template(i).Probability = 0.5;
    i = i + 1;
end

%% T2.1
if find(type_index==4,1)
    slice_template(i).Type = 4;
    slice_template(i).Pattern = FlowPattern.RandomMultiFlow;
    slice_template(i).Weight = 8;
    slice_template(i).NumberPaths = 3;
    slice_template(i).NumberFlows = 120;
    slice_template(i).NumberVNFs = 3;
    slice_template(i).ServiceInterval = 100;
    slice_template(i).Probability = 0.2;    
    i = i + 1;
end
%% T2.2
if find(type_index==5,1)
    slice_template(i).Type = 5;
    slice_template(i).Pattern = FlowPattern.RandomMultiFlow;
    slice_template(i).Weight = 40;
    slice_template(i).NumberPaths = 2;
    slice_template(i).NumberFlows = 40;
    slice_template(i).NumberVNFs = 4;
    slice_template(i).ServiceInterval = 30;
    slice_template(i).Probability = 0.3;    
    i = i + 1;
end
%% T2.3
if find(type_index==6,1)
    slice_template(i).Type = 6;
    slice_template(i).Pattern = FlowPattern.RandomMultiFlow;
    slice_template(i).Weight = 400;
    slice_template(i).NumberPaths = 1;
    slice_template(i).NumberFlows = 4;
    slice_template(i).NumberVNFs = 3;
    slice_template(i).ServiceInterval = 10;
    slice_template(i).Probability = 0.5; 
    i = i + 1;
end
end

%% Adjust Parameters
% If the network is averagely partitioned between network slices, the low weight slice
% will cannot fully utilize the allocated network resources, and the high weight slice
% will run out of resources to achieve the maximal profit, which leads to lower profit 
% of the network.
%
% * *Weight* and *Cost* are in proportion. With one parameter fixed, adjusting another
% parameter can control the optimal serving rate of each type of slice.
% * *Network Capacity*: network capacity determine how many slices can run
% simultaneousely,  which is related to the arriving process and serving interval of each
% slice. For example, if there are at most 10 slices are running at the same time, we can
% set the network capacity so that the network utilization is approaches 1 when there is
% 10 slices.
