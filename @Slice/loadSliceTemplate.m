%% Load Slice Template
%loadSliceTemplate Slice Template Type used to generate slice request.
%   slice_template = loadSliceTemplate(index)

%% Slice Template
% 
%% TODO
% time unit: second;
%% Meaning of digits from right to left:
% 1: experiment number;
% 2: slice type;
% 3: dimensioning trigger: 1 - Event-based, 2 - Time-based, 3 - Threshold-based;
% 4: Ad-hoc mode;
% 5: test mode.
function [ slice_template ] = loadSliceTemplate(type_index)
type_index = unique(type_index, 'stable');
slice_template = struct;
for i = 1:length(type_index)
    switch mod(type_index(i),100)
        case 11      % Type 1 for experiment 1
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
        case 21     % Type 2 for experiment 1
            slice_template(i).FlowPattern = FlowPattern.RandomMultiFlow;
            slice_template(i).Weight = 40;
            slice_template(i).NumberPaths = 2;
            slice_template(i).NumberFlows = 10;
            slice_template(i).NumberVNFs = 3;
            slice_template(i).ServiceInterval = 12;
            slice_template(i).Probability = 1/4;
            slice_template(i).ConstantProfit = 100;
        case 31  % Type 3 for experiment 1
            slice_template(i).FlowPattern = FlowPattern.RandomMultiFlow;
            slice_template(i).Weight = 200;
            slice_template(i).NumberPaths = 3;
            slice_template(i).NumberFlows = 3;
            slice_template(i).NumberVNFs = 4;
            slice_template(i).ServiceInterval = 4;
            slice_template(i).Probability = 3/4;
            slice_template(i).ConstantProfit = 100;
        case 12  % Type 1 for experiment 2
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
            slice_template(i).FlowPattern = FlowPattern.RandomMultiFlow;
            slice_template(i).Weight = 50;
            slice_template(i).NumberPaths = 2;
            slice_template(i).NumberFlows = 40;
            slice_template(i).NumberVNFs = 4;
            slice_template(i).ServiceInterval = 15;
            slice_template(i).Probability = 1/6;
            slice_template(i).MinRate = 5;      % Mbps
            slice_template(i).ConstantProfit = 600;
        case 32  % Type 3 for experiment 2
            slice_template(i).FlowPattern = FlowPattern.RandomMultiFlow;
            slice_template(i).Weight = 300;
            slice_template(i).NumberPaths = 3;
            slice_template(i).NumberFlows = 4;
            slice_template(i).NumberVNFs = 3;
            slice_template(i).ServiceInterval = 6;
            slice_template(i).Probability = 5/6;
            slice_template(i).MinRate = 50;      % Mbps
            slice_template(i).ConstantProfit = 600;
        case 13     % Type 1 for experiment 3
            slice_template(i).FlowPattern = FlowPattern.RandomInterDataCenter;
            slice_template(i).Weight = 40;
            slice_template(i).NumberPaths = 1;
            slice_template(i).NumberFlows = 25;
            slice_template(i).VNFList = [1,2,3];
            slice_template(i).ServiceInterval = 1000;
            slice_template(i).Probability = [];
        case 23      % Type 2 for experiment 3
            slice_template(i).FlowPattern = FlowPattern.RandomInterBaseStation;
            slice_template(i).Weight = 20;
            slice_template(i).NumberPaths = 2;
            slice_template(i).NumberFlows = 200;
            slice_template(i).VNFList = [1,2,4,5];
            slice_template(i).ServiceInterval = 12;
            slice_template(i).Probability = 1/6;
        case 33      % Type 3 for experiment 3
            slice_template(i).FlowPattern = FlowPattern.RandomDataCenter2BaseStation;
            slice_template(i).Weight = 200;
            slice_template(i).NumberPaths = 3;
            slice_template(i).NumberFlows = 10;
            slice_template(i).VNFList = [1,4,6];
            slice_template(i).ServiceInterval = 4;
            slice_template(i).Probability = 5/6;
        case 14     % Type 1 for experiment 4: dynamic flow arrival and departure.
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
        case 24     % Type 2 for experiment 4
            slice_template(i).Weight = 10;
            slice_template(i).NumberPaths = 2;
            slice_template(i).VNFList = [1,2,4,5];
            slice_template(i).ArrivalRate = 1;     % 1 arrivals/hour
            slice_template(i).ServiceInterval = 12;
            slice_template(i).NumberFlows = 100;
            slice_template(i).Flow.ArrivalRate = 1000; % 1000 arrivals/hour  
            slice_template(i).Flow.ServiceInterval = 0.1; 
            slice_template(i).FlowPattern = FlowPattern.RandomInterBaseStation;
        case 34     % Type 3 for experiment 4
            slice_template(i).Weight = 100;
            slice_template(i).NumberPaths = 3;
            slice_template(i).VNFList = [1,4,6];
            slice_template(i).ArrivalRate = 4;   % 4 arrivals/hour
            slice_template(i).ServiceInterval = 4;
            slice_template(i).NumberFlows = 20;
            slice_template(i).Flow.ArrivalRate = 100; % 100 arrivals/hour  
            slice_template(i).Flow.ServiceInterval = 0.2; 
            slice_template(i).FlowPattern = FlowPattern.RandomDataCenter2BaseStation;
        case 44   
            %% Type-1 for experiment 4-2: DynamicNetwork
            % Type-1 is a long term service (1 arrivals/day, 3 days lifetime). So the
            % average number of Type-1 slices is 3. This type can also be used as
            % persistent slice. In that case, the settings of |ArrivalRate| and
            % |SrviceInterval| do not take effect. 
            %
            % The target network <Sample-2> has 15 nodes, so we set the number of flows
            % to: 75 = 1/3*(15*15);
            %
            % Flow arrival rate is set to 4 sec/arrival, i.e., 900 arrivals/hour, for
            % emulating the dynamics of users. The service interval is set to 300s(5min),
            % so that the average number of users stays at 75.
            %
            % The L1 approximation of reconfiguration cost should be normalized, so it
            % will have the same magnitude as the original formulation. 'ReconfigScaler'
            % is an absolute value to serve as the normalizer. We can set this value by
            % first optimizing without normalization, and then computing the ratio between
            % the original reconfiguration cost and the L1 Approximation.
            slice_template(i).Weight = 25;
            slice_template(i).NumberPaths = 1;
            slice_template(i).VNFList = [1,2,3];
            slice_template(i).ArrivalRate = 1/(3600*24);   % 1 arrivals/24 hour => 3600 sec/arrival
            slice_template(i).ServiceInterval = 3600*24*7;
            slice_template(i).NumberFlows = 75;
            slice_template(i).Flow.ArrivalRate = 900/3600; 
            slice_template(i).Flow.ServiceInterval = 300; 
            slice_template(i).FlowPattern = FlowPattern.RandomMultiFlow;
            slice_template(i).ReconfigScaler = 2;       % beta
        case 54
            %% Type-2 for experiment 4-2: DynamicNetwork
            % Type-2 is a middle term service, e.g. inter-data centers communications (1
            % arrival/hour, 6 hours lifetime). So the average number of Type-2 slices is
            % 6. 
            %
            % Assuming that the target network <Sample-2> has 6 data centers, so we set
            % the number of flows to: 18 = 1/2*(6*6);
            %
            % Flow arrival rate is set to 60 sec/arrival, i.e., 60 arrivals/hour, while
            % the service interval is set to 1080s(18min), so that the average number of
            % flows stays at 18. 
            slice_template(i).Weight = 50;
            slice_template(i).NumberPaths = 2;
            slice_template(i).VNFList = [1,2,4,5];
            slice_template(i).ArrivalRate = 1/3600;  % 1 arrival/hour => 3600 sec/arrival
            slice_template(i).ServiceInterval = 3600*6; % 
            slice_template(i).NumberFlows = 18;
            slice_template(i).Flow.ArrivalRate = 60/3600;  % 60 arrivals/hour  
            slice_template(i).Flow.ServiceInterval = 1080; 
            slice_template(i).FlowPattern = FlowPattern.RandomMultiFlow;
            slice_template(i).ReconfigScaler = 2.4;
        case 64 
            %% Type-3 for experiment 4-2: DynamicNetwork
            % Type-3 is a middle term service with high QoS-demand (10 arrivals/hour, 1
            % hour lifetime). So the average number of Type-3 slices is 10. 
            %
            % While the service has stringent QoS requirement, the number of flows is
            % relatively fewer, thus set to 4.
            %
            % Those flows always arrive at the creation of the slice and depart when the
            % slice is released. Thus the arrival rate of flow is relatively lower. Flow
            % arrival rate is set to 450 sec/arrival, i.e., 8 arrivals/hour, while  
            % the service interval is set to 1800s(30min), so that the average number of
            % flows stays at 4. 
            slice_template(i).Weight = 300;
            slice_template(i).NumberPaths = 3;
            slice_template(i).VNFList = [1,4,6];
            slice_template(i).ArrivalRate = 10/3600;   % 4 arrivals/hour
            slice_template(i).ServiceInterval = 3600;
            slice_template(i).NumberFlows = 4;
            slice_template(i).Flow.ArrivalRate = 8/3600; % 8 arrivals/hour  
            slice_template(i).Flow.ServiceInterval = 1800; 
            slice_template(i).FlowPattern = FlowPattern.RandomMultiFlow;
            slice_template(i).ReconfigScaler = 0.5;
        otherwise
            error('error: unrecognized slice type (%d).', type_index(i));
    end
    slice_template(i).Type = type_index(i);

    switch getdigit(type_index(i),3)
        case 1
            slice_template(i).Trigger = 'EventBased';
            switch mod(type_index(i),100)
                case 44
                    slice_template(i).EventInterval = 10;       % {50}
                case 54
                    slice_template(i).EventInterval = 10;
                case 64
                    slice_template(i).EventInterval = 4;
                otherwise
                    error('error: unrecognized slice type (%d).', type_index(i));
            end
        case 2
            slice_template(i).Trigger = 'TimeBased';
            switch mod(type_index(i),100)
                case 44
                    slice_template(i).TimeInterval = 100;
                case 54
                    slice_template(i).TimeInterval = 300;
                case 64
                    slice_template(i).TimeInterval = 900;
                otherwise
                    error('error: unrecognized slice type (%d).', type_index(i));
            end
        case 3
            slice_template(i).Trigger = 'ThresholdBased';
            %         otherwise no trigger.
        case 0
        otherwise
            error('error: unrecognized slice type (%d).', type_index(i));
    end
    switch getdigit(type_index(i),4)   % Ad-hoc mode
        case 1
            %%
            % If Ad-hoc mode is set, flows may originate from uncovered nodes. Unused nodes and links
            % will be recycled.
            % Otherwise, flows only originated within the existing virtual nodes. Unused node and link
            % resources will be recycled, but the nodes and links themselves will be maintained.
            slice_template(i).Adhoc = true;
        case 0
        otherwise
            error('error: unrecognized slice type (%d).', type_index(i));
    end
    
    switch getdigit(type_index(i),5)
        case 1
            switch mod(type_index(i),100)
                case {44,54,64}
                    slice_template(i).ClassName = 'DynamicSliceTest';
            end
        case 0
        otherwise
            error('error: unrecognized slice type (%d).', type_index(i));
    end
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
% simultaneously,  which is related to the arriving process and serving interval of each
% slice. For example, if there are at most 10 slices are running at the same time, we can
% set the network capacity so that the network utilization is approaches 1 when there is
% 10 slices.
