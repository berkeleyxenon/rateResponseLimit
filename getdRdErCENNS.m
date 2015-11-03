function C = getdRdErCENNS(C);
% C = getdRdErCENNS(C)
% C is spec'd in defineExperiment.m
%
% dRdEr units cts/keV/kg/day
%
%
% 121212 pfs
% 130306 pfs - switch to standard cts/kev/kg/day
% 130910 pfs switch C.ma =>C.mp for consistency with DM code
% 150902 pfs dusting off for LZ neutrino Task Force

[rr,cc]=size(C.Er);
if rr>cc % flip dimensionality
	C.Er = C.Er';
end

%% Fundamental constants
	%C.G_F = 1.16639*1e-5; % in GeV^-2 % this is GF/(hbarc)^3
	C.mp = 0.938; % % proton Mass, GeV
	C.hbar = 6.5822*1e-25; % in GeV * s
	C.c = 2.9979e10; % cm/s
	C.N_A = 6.02e23;
	C.M_N = C.A*C.mp; % GeV
	C.N = round(C.A)-C.Z;

Emin = sqrt(C.M_N .* C.Er/2); % min neutrino energy to cause a recoil of energy Er
C = getHelmFF(C);

%% neutrino flux
%if ~isfield('C','source');C.source=1;end;
mm=1; % which solar model
phi=getSolarNuFlux(mm);

switch C.source
case 1 % 8B
	%load('1996.B8.8hep.solar.neutrinos.mat');
	if 1 % read in Bahcall et al. txt file
		fid = fopen('./dataFiles/b8spectrum.txt','r');
		headr = fscanf(fid,'%s',82);
		for ii=1:829
			B8.Enu(ii) = str2double(fscanf(fid,'%s',1));
			B8.lambda(ii) = str2double(fscanf(fid,'%s',1));
			B8.lambdap(ii) = str2double(fscanf(fid,'%s',1)); % +3sigma
			B8.lambdam(ii) = str2double(fscanf(fid,'%s',1)); % -3sigma
			ii=ii+1;
		end
		fclose(fid);	
		norm = sum(B8.lambda)*0.02; % = 1
	end
	%phi_8B = 5.79e6; % flux (cts/cm^2/sec) from model BP04(Yale) from ApJ 621 L85 (2005)
	%phi_8B = 7.4e6; % agrees w Strigari
	Enu = B8.Enu(2:end)';dEnu = Enu(2)-Enu(1);
	dN_dEnu = B8.lambda(2:end)' * phi.B8; % (cts/0.02 MeV) * 0.02 MeV * (cm^2/sec)^-1
	%dN_dEnu = B8.lambda(2:end)'*phi.B8;
case 2 % hep
	fid=fopen('./dataFiles/hep.txt');
	count=0;
	for z=1:500
		for i=1:2
			count=count+1;
			hep.x(count) = str2num(fscanf(fid,'%s',1));
			hep.y(count) = str2num(fscanf(fid,'%s',1));
		end
	end
	fclose(fid);
	Enu=hep.x';dEnu=Enu(2)-Enu(1);dN_dEnu=hep.y' * phi.hep;
case 3 % atmospheric nu
	switch 1
	case 0 % strigari rip
		rip = [13.2 117.4; 19.7 151.4; 26.8 168.3; 30.2 179.5; 36 180; 46 174; 58 134; 74 115; 103 86; 128 66; 170 47; 235 28.3; 378 12.3; 556 6.0; 983 1.67]; % ripped points are nu_mu flux
		x = rip(:,1);y = rip(:,2)*1e-4; % cts/cm^2/s/MeV
		dEnu=1;Enu = [x(1):dEnu:x(end)]';
	case 1 % original Astropart Phys 23 (2005) 526
		x = [0.013 0.015 0.017 0.019 0.021 0.024 0.027 0.030 0.033 0.038 0.042 0.047 0.053 0.060 0.067 0.075 0.084 0.094 0.106 0.119 0.133 0.150 0.168 0.188 0.211 0.237 0.266 0.299 0.335 0.376 0.422 0.473 0.531 0.596 0.668 0.750 0.841 0.944]*1e3; % MeV
		switch 1
		case 0 % Super Kamiokande site, nu_mu and nu_e
			ym = [0.114 0.124 0.138 0.146 0.155 0.159 0.164 0.181 0.174 0.179 0.178 0.176 0.153 0.131 0.123 0.114 0.107 0.0963 0.0842 0.0727 0.0635 0.0552 0.0477 0.0412 0.0344 0.0284 0.0236 0.0196 0.0158 0.0128 0.0103 0.0082 0.00649 0.00515 0.00398 0.00313 0.00241 0.00182]*1e6/1e7; % n/cm^2/s/MeV	
			ye = [0.0696 0.0746 0.0797 0.0874 0.0942 0.101 0.103 0.109 0.108 0.107 0.101 0.0885 0.0696 0.0644 0.0593 0.0543 0.0497 0.0451 0.0406 0.0358 0.0317 0.0273 0.0239 0.0204 0.0170 0.0145 0.0120 0.00996 0.00811 0.00662 0.00527 0.00423 0.00337 0.00266 0.00209 0.00162 0.00124 0.00095]*1e6/1e7; % n/cm^2/s/MeV	
		case 1	% Gran Sasso site, nu_mu and nu_e
			ym = [0.174 0.190 0.211 0.220 0.235 0.243 0.248 0.275 0.263 0.270 0.269 0.266 0.232 0.197 0.182 0.169 0.157 0.141 0.123 0.105 0.0912 0.0785 0.0672 0.0575 0.0476 0.0391 0.0325 0.0262 0.0212 0.0169 0.0134 0.0105 0.00825 0.00640 0.00490 0.00378 0.00287 0.00214]*1e6/1e7; % n/cm^2/s/MeV
			ye = [0.105 0.114 0.121 0.133 0.140 0.152 0.158 0.165 0.164 0.161 0.153 0.134 0.105 0.096 0.0885 0.0808 0.0735 0.0671 0.0595 0.0519 0.0457 0.0394 0.0340 0.0292 0.0241 0.0203 0.0168 0.0140 0.0110 0.00894 0.00706 0.00557 0.00445 0.00344 0.00267 0.00211 0.00153 0.00116]*1e6/1e7; % n/cm^2/s/MeV
		end
		y = (ym+ye); % sum since coherent scatter
		%y = y*0.67; % guess for SURF
	end
	dEnu=1;Enu = [x(1):dEnu:x(end)]';
	dN_dEnu = interp1(x,y,Enu,'pchip'); % cts/cm^2/s/MeV
	%figure(1);clf;loglog(x,y,'o');hold on;plot(Enu,dN_dEnu,'r-');
	%dN_dEnu = dN_dEnu * 1.57; % simple approx to add in nu_e contribution, to be fixed...
case 4 % DSN
	% rip from Horiuchi 2009 PRD
	dEnu=1;Enu=[2:dEnu:50]';
	% 8 MeV
	x = [2:10 12 15 20 30 40];
	y = [0.51 0.83 1.03 1.15 1.17 1.14 1.08 1.01 0.93 0.75 0.53 0.28 0.08 0.02];
	dN_dEnu8 = interp1(x,y,Enu,'pchip',0); % cts/cm^2/s/MeV
	dN_dEnu8 = dN_dEnu8;
	% 6 MeV
	x = [2 3 4 5 7 10 20 30 38];
	y = [1.3 1.8 2.05 2.06 1.75 1.15 0.21 0.037 0.010];
	dN_dEnu6 = interp1(x,y,Enu,'pchip',0); % cts/cm^2/s/MeV
	% 4 MeV
	x = [2 3 4 5 10 20 28];
	y = [4.08 4.67 4.38 3.69 1.11 0.082 0.010];
	dN_dEnu4 = interp1(x,y,Enu,'pchip',0); % cts/cm^2/s/MeV

	% add em up
	dN_dEnu = dN_dEnu4 + dN_dEnu6 + dN_dEnu8 * 4; % factor of 4, yo
	%dN_dEnu = dN_dEnu * 0.5; % lower family of curves
end
figure(11);plotlogstairs(Enu,dN_dEnu,'k-');hold on;box on;grid on;axis([0 100 1e-4 1e6]);

%% coherent cross section
Gf = 1.166e-5; % Fermi Coupling Constant, units GeV^-2
Qw = C.N - C.Z*(1-4*0.231); % assume sin^2(thetaW) = 0.231;

m = size(Enu,1);
n = size(C.Er,2);

dsigma_dEr = Gf.^2 /(4*pi) * Qw^2 * C.M_N ...
			* (1 - C.M_N.*repmat(C.Er,m,1)./(2*repmat(Enu,1,n).^2)) ... % GeV*keV/(MeV)^2
			.* repmat(C.F,m,1).^2 ... % form factor
			.* (C.hbar*C.c)^2 / 1e6 ... % throw in hbarc^2 / 1e6 to get units correct
			; % cm^2/keV
			

%% calculate expected event rate -- folowing NJP 11 105011 (2009)
dRdEr =  sum( ...
			repmat(dN_dEnu,1,n) * dEnu ... % neutrino flux, cts/MeV/cm^2/sec * MeV
			.* dsigma_dEr ... % neutrino cross section, cm^2/keV
			.* (repmat(Enu,1,n)>repmat(Emin,m,1)) ... % integral limits (integrate over MeV)
			); % cts/keV/nucleus/sec

%dRdEr = dRdEr ...
%	* 3600*24*365 ... % convert seconds => year
%	* C.N_A/C.A *1e6 ... % convert nuclei => tonne
%	; % now we have cts/keV/tonne/yr

dRdEr = dRdEr ...
	* 3600*24 ... % convert seconds => day
	* C.N_A/C.A *1e3 ... % convert nuclei => kg
	; % now we have cts/keV/kg/day

C.dRdEr = dRdEr;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END

if 0 % verify agreement with NJP 11 105011 (2009)
	figure(2);clf;
	loglog(C.Er,dRdEr ,'r-'); % 
	axis([0.1 10 1e-5 1e4]);
	set(gca,'ytick',10.^[-5:1:4]);
	grid on;
end

if 0 % verify agreement with neutrino spectrum
	figure(8);clf;
	loglog(B8.Enu,B8.lambda*5.79e6,'r-');
	axis([0.1 1e3 1e-4 1e6]);
	set(gca,'xtick',10.^[-1:1:3]);
	set(gca,'ytick',10.^[-4:1:6]);
end

