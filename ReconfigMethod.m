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
	
	methods
		function str = FullName(this)
			str = this.char;
			switch this
				case ReconfigMethod.Baseline
					str = 'Baseline Reconfiguration';
				case ReconfigMethod.Fastconfig
					str = 'Fast Reconfiguration';
				case ReconfigMethod.FastconfigReserve
					str = 'Fast Reconfiguration with Resource Reservation';
				case ReconfigMethod.Fastconfig1
				case ReconfigMethod.Fastconfig1Reserve
				case ReconfigMethod.Fastconfig2
					str = 'Fast Reconfiguration 2';
				case ReconfigMethod.Fastconfig2Reserve
				case ReconfigMethod.Dimconfig
					str = 'Hybrid Slicing Scheme';
				case ReconfigMethod.DimconfigReserve
					str = 'Hybrid Slicing Scheme with Resource Reservation';
				case ReconfigMethod.DimconfigReserve0
				case ReconfigMethod.Dimconfig1
				case ReconfigMethod.Dimconfig2
					str = 'Hybrid Slicing Scheme 2';
				case ReconfigMethod.DimBaseline
					str = 'Baseline Reconfiguration With Dimensioning';
			end
		end
	end
	
	methods (Static)
		function obj = Enumerate(str)
			m = enumeration('ReconfigMethod');
			for i = 1:length(m)
				if m(i) == str
					obj = m(i);
					return;
				end
			end
			error('error: unidentified entry.');
		end
		
	end
	
end

