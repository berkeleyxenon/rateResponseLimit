function phi=getSolarNuFlux(mm)
	% flux from solar models in 2005.ApJ.621.L85.Bahcall
	switch mm
	case 1 % BP04(Yale)
		phi.pp = 5.94e10; %/cm^2/s
		phi.hep = 7.88e3; %/cm^2/s
		phi.Be7 = 4.86e9; %/cm^2/s
		phi.N13 = 5.71e8; %/cm^2/s
		phi.B8 = 5.79e6; %/cm^2/s
	case 2 % BP04(Garching)
		phi.pp = 5.94e10; %/cm^2/s
		phi.hep = 7.88e3; %/cm^2/s
		phi.Be7 = 4.84e9; %/cm^2/s
		phi.N13 = 5.70e8; %/cm^2/s
		phi.B8 = 5.74e6; %/cm^2/s
	case 3 % BS04
		phi.pp = 5.94e10; %/cm^2/s
		phi.hep = 7.86e3; %/cm^2/s
		phi.Be7 = 4.88e9; %/cm^2/s
		phi.N13 = 5.62e8; %/cm^2/s
		phi.B8 = 5.87e6; %/cm^2/s
	case 4 % BP05(14N)
		phi.pp = 5.99e10; %/cm^2/s
		phi.hep = 7.91e3; %/cm^2/s
		phi.Be7 = 4.89e9; %/cm^2/s
		phi.N13 = 3.11e8; %/cm^2/s
		phi.B8 = 5.83e6; %/cm^2/s
	case 5 % BP05(OP)
		phi.pp = 5.99e10; %/cm^2/s
		phi.hep = 7.93e3; %/cm^2/s
		phi.Be7 = 4.84e9; %/cm^2/s
		phi.N13 = 3.07e8; %/cm^2/s
		phi.B8 = 5.69e6; %/cm^2/s
	case 6 % BP05(AGS,OP)
		phi.pp = 6.06e10; %/cm^2/s
		phi.hep = 8.25e3; %/cm^2/s
		phi.Be7 = 4.34e9; %/cm^2/s
		phi.N13 = 2.01e8; %/cm^2/s
		phi.B8 = 4.51e6; %/cm^2/s
	case 7 % BP05(AGS,OPAL)
		phi.pp = 6.06e10; %/cm^2/s
		phi.hep = 8.23e3; %/cm^2/s
		phi.Be7 = 4.38e9; %/cm^2/s
		phi.N13 = 2.03e8; %/cm^2/s
		phi.B8 = 4.59e6; %/cm^2/s
	end
	
