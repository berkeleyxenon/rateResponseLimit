% quick set up and plot
%
% 130110 pfs

%%%%%%%%%%%%%%%%%%%%%%%%% Fundamental constants
	C.mp = 0.938; % % proton Mass, GeV
	C.hbar = 6.5822*1e-25; % in GeV * s
	C.c = 2.9979e10; % cm/s
	C.N_0 = 6.022e23; % Avogadro number, atoms/mol
%%%%%%%%%%%%%%%%%%%%%%%%% Astrophysical constants
	C.rho_chi = 0.3; % GeV cm-3 c-2
	C.v_0 = 220/3e5; % velocity dispersion of isotropic MB distribution
	% rotational speed for LSR is similar to M-B velocity dispersion -- see getDMvelocityDistribution.m
	switch 100
	case 1
		C.v_esc = 550/3e5; % units of c -- following PRD 79 043513 (2009)
	case 100
		C.v_esc = 544/3e5; % units of c -- following PRD 79 043513 (2009)
	end
%%%%%%%%%%%%%%%%%%%%%%%%% Nuclear physics constants
	C.f_p = 1;
	C.f_n = 1;
%%%%%%%%%%%%%%%%%%%%%%%%% the energy points for the expected spectrum
	C.dEr = 0.01; % should be sufficient granularity...
	C.Er = C.dEr/2:C.dEr:100;

%%%%%%%%%%%%%%%%%%%%%%%%% time of year
	C.t = 2/3; % basically average value over the year
	%C.t = 0.91; % minimum flux
	C.t = 0.49; % LUX flux

%%%%%%%%%%%%%%%%%%%%%%%%%% particulars
	C.m_chi = 2; % GeV
	C.delta = 0; % keV
	C.sigma_n = 5e-40; % cm^2

%%%%%%%%%%%%%%%%%%%%%%%%%% experiment
	switch 2
	case 0
		C.A = 131.3;
		C.Z = 54;
        ax=[0 4 1e-2 1e3];
	case 1
		C.A = 72.6;
		C.Z = 32;	
        ax=[0 4 1e-2 1e3];
	case 2
		C.A = 28.0;
		C.Z = 14;
        ax=[0 10 1e-1 50];
    end
	C.liveDays = 88.3;
	C.kg = 118;


C.M_N = C.A*C.mp; % GeV

C = getHelmFF(C);
C = getBetaMin(C); % calculate beta_min
dR_dEr = getdRdErDM(C); % dru
cts = dR_dEr *  C.dEr * C.liveDays * C.kg;

if 1
figure(1);%clf;
	semilogy(C.Er,dR_dEr,'k-');
	hold on;
	axis(ax);
	set(gca,'xtick',[0:1:10],'ytick',10.^[-6:1:6]);
	xlabel('recoil energy [keV]');
	ylabel('dru');
	
	setplot;
end

% title('$\sigma=5\times 10^{-40}$ cm$^2$ on silicon')
% how many events?
%sum(dR_dEr(C.Er>3))*C.dEr*118*88.3
