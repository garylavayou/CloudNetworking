%% Load Slice Template
%loadSliceTemplate Slice Template Type used to generate slice request.
%   slice_template = loadSliceTemplate(index)

%% Slice Template
% 
function [ slice_template ] = loadSliceTemplate(type_index)
type_index = unique(type_index);
slice_template = struct;
for i = 1:length(type_index)
    switch type_index(i)
        case 11      % Type 1 for experiment 1
            slice_template(i).Type = 11;
            slice_template(i).FlowPattern = FlowPattern.RandomMultiFlow;
            slice_template(i).Weight = 15;
            slice_template(i).NumberPaths = 1;
            slice_template(i).NumberFlows = 30;
            slice_template(i).NumberVNFs = 2;
            slice_template(i).ServiceInterval = 10000;
            slice_template(i).Probability = [];
            slice_template(i).ConstantProfit = 100;
            % slcie_opt.DelayConstraint = [];
            % slice_opt.VNFList = [];
        case 12     % Type 2 for experiment 1
            slice_template(i).Type = 12;
            slice_template(i).FlowPattern = FlowPattern.RandomMultiFlow;
            slice_template(i).Weight = 40;
            slice_template(i).NumberPaths = 2;
            slice_template(i).NumberFlows = 10;
            slice_template(i).NumberVNFs = 3;
            slice_template(i).ServiceInterval = 12;
            slice_template(i).Probability = 1/4;
            slice_template(i).ConstantProfit = 100;
        case 13  % Type 3 for experiment 1
            slice_template(i).Type = 13;
            slice_template(i).FlowPattern = FlowPattern.RandomMultiFlow;
            slice_template(i).Weight = 200;
            slice_template(i).NumberPaths = 3;
            slice_template(i).NumberFlows = 3;
            slice_template(i).NumberVNFs = 4;
            slice_template(i).ServiceInterval = 4;
            slice_template(i).Probability = 3/4;
            slice_template(i).ConstantProfit = 100;
        case 21  % Type 1 for experiment 2
            slice_template(i).Type = 21;
            slice_template(i).FlowPattern = FlowPattern.RandomMultiFlow;
            slice_template(i).Weight =25;
            slice_template(i).NumberPaths = 1;
            slice_template(i).NumberFlows = 120;    % at most 15*14 flows
            slice_template(i).NumberVNFs = 3;
            slice_template(i).ServiceInterval = 1000;
            slice_template(i).Probability = [];
            slice_template(i).MinRate = 1;      % Mbps
            slice_template(i).ConstantProfit = 600;
        case 22  % Type 2 for experiment 2
            slice_template(i).Type = 22;
            slice_template(i).FlowPattern = FlowPattern.RandomMultiFlow;
            slice_template(i).Weight = 50;
            slice_template(i).NumberPaths = 2;
            slice_template(i).NumberFlows = 40;
            slice_template(i).NumberVNFs = 4;
            slice_template(i).ServiceInterval = 15;
            slice_template(i).Probability = 1/6;
            slice_template(i).MinRate = 5;      % Mbps
            slice_template(i).ConstantProfit = 600;
        case 23  % Type 3 for experiment 2
            slice_template(i).Type = 23;
            slice_template(i).FlowPattern = FlowPattern.RandomMultiFlow;
            slice_template(i).Weight = 300;
            slice_template(i).NumberPaths = 3;
            slice_template(i).NumberFlows = 4;
            slice_template(i).NumberVNFs = 3;
            slice_template(i).ServiceInterval = 6;
            slice_template(i).Probability = 5/6;
            slice_template(i).MinRate = 50;      % Mbps
            slice_template(i).ConstantProfit = 600;
        case 31     % Type 1 for experiment 3
            slice_template(i).Type = 31;
            slice_template(i).FlowPattern = FlowPattern.RandomInterDataCenter;
            slice_template(i).Weight = 40;
            slice_template(i).NumberPaths = 1;
            slice_template(i).NumberFlows = 25;
            slice_template(i).VNFList = [1,2,3];
            slice_template(i).ServiceInterval = 1000;
            slice_template(i).Probability = [];
        case 32      % Type 2 for experiment 3
            slice_template(i).Type = 32;
            slice_template(i).FlowPattern = FlowPattern.RandomInterBaseStation;
            slice_template(i).Weight = 10;
            slice_template(i).NumberPaths = 2;
            slice_template(i).NumberFlows = 100;
            slice_template(i).VNFList = [1,2,4,5];
            slice_template(i).ServiceInterval = 12;
            slice_template(i).Probability = 1/6;
        case 33      % Type 3 for experiment 3
            slice_template(i).Type = 33;
            slice_template(i).FlowPattern = FlowPattern.RandomDataCenter2BaseStation;
            slice_template(i).Weight = 100;
            slice_template(i).NumberPaths = 3;
            slice_template(i).NumberFlows = 20;
            slice_template(i).VNFList = [1,4,6];
            slice_template(i).ServiceInterval = 4;
            slice_template(i).Probability = 5/6;
        case 41     % Type 1 for experiment 4: dynamic flow arrival and departure.
            slice_template(i).Type = 41;
            slice_template(i).Weight = 40;
            slice_template(i).NumberPaths = 1;
            slice_template(i).VNFList = [1,2,3];
            slice_template(i).ArrivalRate = 0.05;   % 0.05 arrivals/hour  =>  20 hours/arrival
            slice_template(i).ServiceInterval = 100;
            %%%
            % Flow parameters for the slice.
            % Since the slice is dynamic, the number of flows is not fixed. |NumberFlows|
            % specify the average number of flows of this slice, add the initial number of
            % flows in the slice. |ArrivalRate| and |ServiceInterval| specifiy the flow's
            % dynamic attributes.
            slice_template(i).NumberFlows = 40;         % 25/45
            slice_template(i).Flow.ArrivalRate = 15; % 25 arrivals/hour  
            slice_template(i).Flow.ServiceInterval = 1; 
            slice_template(i).FlowPattern = FlowPattern.RandomInterDataCenter;
        case 42     % Type 2 for experiment 4
            slice_template(i).Type = 42;
            slice_template(i).Weight = 10;
            slice_template(i).NumberPaths = 2;
            slice_template(i).VNFList = [1,2,4,5];
            slice_template(i).ArrivalRate = 1;     % 1 arrivals/hour
            slice_template(i).ServiceInterval = 12;
            slice_template(i).NumberFlows = 100;
            slice_template(i).Flow.ArrivalRate = 1000; % 1000 arrivals/hour  
            slice_template(i).Flow.ServiceInterval = 0.1; 
            slice_template(i).FlowPattern = FlowPattern.RandomInterBaseStation;
        case 43     % Type 3 for experiment 4
            slice_template(i).Type = 43;
            slice_template(i).Weight = 100;
            slice_template(i).NumberPaths = 3;
            slice_template(i).VNFList = [1,4,6];
            slice_template(i).ArrivalRate = 4;   % 4 arrivals/hour
            slice_template(i).ServiceInterval = 4;
            slice_template(i).NumberFlows = 20;
            slice_template(i).Flow.ArrivalRate = 100; % 100 arrivals/hour  
            slice_template(i).Flow.ServiceInterval = 0.2; 
            slice_template(i).FlowPattern = FlowPattern.RandomDataCenter2BaseStation;
        case 44     % Type 1 for experiment 4-2: DynamicNetwork
            slice_template(i).Type = 44;
            slice_template(i).Weight = 25;
            slice_template(i).NumberPaths = 1;
            slice_template(i).VNFList = [1,2,3];
            slice_template(i).ArrivalRate = 0.01;   % 0.01 arrivals/min  =>  100 min/arrival
            slice_template(i).ServiceInterval = 300;
            slice_template(i).NumberFlows = 80;
            slice_template(i).Flow.ArrivalRate = 120; % 1/120 min/arrival
            slice_template(i).Flow.ServiceInterval = 0.5; 
            slice_template(i).FlowPattern = FlowPattern.RandomMultiFlow;
        case 45     % Type 2 for experiment 4-2: DynamicNetwork
            slice_template(i).Type = 45;
            slice_template(i).Weight = 40;
            slice_template(i).NumberPaths = 2;
            slice_template(i).VNFList = [1,2,4,5];
            slice_template(i).ArrivalRate = 1;     % 1 arrivals/hour
            slice_template(i).ServiceInterval = 12;
            slice_template(i).NumberFlows = 20;
            slice_template(i).Flow.ArrivalRate = 1000; % 1000 arrivals/hour  
            slice_template(i).Flow.ServiceInterval = 0.1; 
            slice_template(i).FlowPattern = FlowPattern.RandomMultiFlow;
        case 46     % Type 3 for experiment 4-2: DynamicNetwork
            slice_template(i).Type = 46;
            slice_template(i).Weight = 200;
            slice_template(i).NumberPaths = 3;
            slice_template(i).VNFList = [1,4,6];
            slice_template(i).ArrivalRate = 4;   % 4 arrivals/hour
            slice_template(i).ServiceInterval = 4;
            slice_template(i).NumberFlows = 5;
            slice_template(i).Flow.ArrivalRate = 100; % 100 arrivals/hour  
            slice_template(i).Flow.ServiceInterval = 0.2; 
            slice_template(i).FlowPattern = FlowPattern.RandomMultiFlow;
        otherwise
            error('error: unrecognized slice type.')
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
