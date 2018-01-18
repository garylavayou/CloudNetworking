global DEBUG g_results event_num; 
event_num = 0;
g_results = table;
DEBUG = true;

%% TODO
% Add event listener to SFED, PN, and slices in PN


%% Event handling
% Event dispatcher's listeners
SFED.AddListener(PN, {'SliceArrive', 'SliceDepart', 'FlowArrive', 'FlowDepart'}, ...
    @PN.eventhandler);
% Network's listeners
% There will be slices listen to the network, after they are added to the network.
PN.AddListener(SFED, {'AddSliceSucceed', 'AddSliceFailed', 'RemoveSliceSucceed', ...
    'RemoveSliceFailed', 'AddFlowSucceed', 'AddFlowFailed', 'RemoveFlowSucceed', ...
    'RemoveFlowFailed'}, @SFED.eventhandler);
for temp_sid = 1:PN.NumberSlices
    sl = PN.slices{temp_sid};
    if sl.isDynamicFlow
        PN.AddListener(sl, {'FlowArrive', 'FlowDepart'}, @sl.eventhandler);
        sl.AddListener(PN, {'AddFlowSucceed', 'AddFlowFailed', ...
            'RemoveFlowSucceed', 'RemoveFlowFailed',...
            'RequestDimensioning', 'DeferDimensioning'}, @PN.eventhandler);
        sl.setOption('ReconfigMethod', options.ReconfigMethod);
    end
end

%% Main loop
iter_num = 0;
while iter_num < NUM_EVENT
    %%%
    % *output*
    iter_num = iter_num + 1;
    total_iter_num = total_iter_num + 1;
    waitbar(total_iter_num/TOTAL_NUM, progress_bar, ...
            sprintf('Simulation Progress: %d/%d', total_iter_num, TOTAL_NUM));
    ev = SFED.nextEvent;
end