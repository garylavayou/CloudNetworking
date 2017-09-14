%% Cloud Networking
%% Dimensioning Network Slicing by Adjusting Resource Price
% The objective is to maximize the net social welfare, at the same time each network slice
% want to maximize their own net profit. 

%%%
% *Note*
%
% # Compared with the shortest path, longer path will be sold at a chiper price to
% encourage the virtual slice operator to use it, since the slice owner will pay for more
% links (and to guarantee the latency constraint, the physical owner should allocate more
% bandwidth to the slice owner.)
% # The pricing should distinguish different service type (e.g. different latency).
% Usually, stringent latency constraint correpsonds to high user utility. For example,
% considering file downloading, watching online videos and playing online games. Playing
% online games has most stringent delayconstraint, it has highest utility, while file
% downloading almost has no delay constraint, so it has lowest utility.
% # The cost function of underlying resources add a penalty. 