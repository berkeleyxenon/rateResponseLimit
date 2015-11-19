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
	C.Er = C.dEr/2:C.dEr:20;

%%%%%%%%%%%%%%%%%%%%%%%%% time of year
	C.t = 2/3; % basically average value over the year
	%C.t = 0.91; % minimum flux
	%C.t = 0.49; % LUX flux

%%%%%%%%%%%%%%%%%%%%%%%%%% particulars
	m_chi = [1 2 5 10 50]; % GeV
	C.delta = 0; % keV
	C.sigma_n = 1e-44; % cm^2

hleg=[];leg=[];
for ii=1:4
	C.m_chi = m_chi(2);	
	%%%%%%%%%%%%%%%%%%%%%%%%%% experiment
		switch ii
		case 1
			C.A = 131.3;
			C.Z = 54;
		case 2
			C.A = 72.6;
			C.Z = 32;	
		case 3
			C.A = 28.0;
			C.Z = 14;
		case 4
			C.A = 16.0;
			C.Z = 10;
		end
		C.liveDays = 30;
		C.kg = 5600;


	C.M_N = C.A*C.mp; % GeV
	%sty={'r-' 'k--' 'b-.' 'g-' 'm-'};
	define_rainbow;fn=fieldnames(cols);
	C = getHelmFF(C);
	C = getBetaMin(C); % calculate beta_min
	dR_dEr = getdRdErDM(C); % dru
	%cts = dR_dEr *  C.dEr * C.liveDays * C.kg;

	if 1
	figure(1);%clf;
		h=plot(C.Er,dR_dEr,'-');set(h,'color',0.9*cols.(fn{ii}));hleg(end+1)=h(1);
		hold on;
		%Qy = 4; % broad brush :)
		%h=plot(C.Er,PoissonConvolution(dR_dEr,C.Er.*Qy,1),'--');set(h,'color',0.9*cols.(fn{ii}));
		ax=[0.1 15 1e-6 0.1];
		axis(ax);
		%set(gca,'xtick',[0:1:15],'ytick',10.^[-6:1:6]);
		xlabel('nuclear recoil energy / keV');
		ylabel('dru');
		set(gca,'xsc','log','ysc','log');
		setplot;
		leg{end+1}=dis('A~%1.0f',C.A);
		legend(hleg,leg,'location','ne');
	end
end

% title('$\sigma=5\times 10^{-40}$ cm$^2$ on silicon')
% how many events?
%sum(dR_dEr(C.Er>3))*C.dEr*118*88.3
