classdef ReconfigMethod < uint32
    
    enumeration
        Baseline(0); 
        Fastconfig(1);
        FastconfigReserve(2);
        Fastconfig1(3);
        Fastconfig2(4);
        Fastconfig1Reserve(5);
        Fastconfig2Reserve(6);
        % Periodically performing slice dimensioning
        % Reconfiguring slice with the aim to reduce reconfiguration cost.
        Dimconfig(11);
        % Periodically performing slice dimensioning
        % Reconfiguring slice with the aim to reduce reconfiguration cost. 
        % Compared with DIMCONFIG, this method further provides mechanisms, i.e., resource
        % reservation and partial reconfiguration to reduce the amount of recofigurations. 
        DimconfigReserve(12);
        DimconfigReserve0(13);  % implicit resource reservation.
        Dimconfig1(14);
        Dimconfig2(15);
        % Periodically performing slice dimensioning
        % Reconfiguring slice without considering reconfiguration cost. 
        DimBaseline(101);
    end

end

