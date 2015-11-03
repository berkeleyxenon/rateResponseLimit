%
%
% 100615 pfs
% 100617 fixed S1 coincidence treatment
% 120706 pfs - checking Zep III result
% 120730 pfs - checking XENON100 result
% 121210 pfs - updated to include cnns
% 130710 pfs - updated for LUX

%typ=10 % leave this at 0.1 unless you know for sure that you want it otherwise :)
pathh = '/';

%%% simulation params
	switch typ
	case 0.1
		Er_sim_pts = [1:0.25:100]; % 40 is ~enough for up to S1=20 (studying CENNS from 8B)		
		load('./dataFiles/XENON10_single_elastic_AmBe_fiducial.mat');
		Nstats0 = hist(MCdata,Er_sim_pts); Nstats0(1) = Nstats0(2); Nstats0(end)=0;
	case 0.2 % flat spectrum
		%Er_sim_pts = [1:0.25:50]; % OK up to ~25 phe not higher
		Er_sim_pts = [1:0.1:60]; % OK up to ~30 phe not higher
		Nstats0 = 100*ones(size(Er_sim_pts));
	case 0.3 % flat spectrum, YBe like
		Er_sim_pts = [1:0.001:4.5];
		Nstats0 = 100*ones(size(Er_sim_pts));
	case 0.31 % flat spectrum, BiBe like
		Er_sim_pts = [1:0.00025:2.7];
		Nstats0 = 200*ones(size(Er_sim_pts));
	case 0.32 % flat spectrum, AlBe like
		Er_sim_pts = [1:0.1:3.8];
		Nstats0 = 100*ones(size(Er_sim_pts));
	case 0 % neutrino-nucleus scattering
		C.Z = 54;
		C.A = 131.3; % average		
		C.dEr = 0.001;
		%C.Er = [C.dEr:C.dEr:10];
		switch source
		case 1
			C.source=1;rate=12;
			C.Er = [1:C.dEr:10];
		case 2
			C.source=2;rate=13;
			C.Er = [1:C.dEr:20];
		case 3
			C.source=3;rate=14;
			C.Er = [1:C.dEr:50];
		case 4
			C.source=4;rate=14;
			C.Er = [1:C.dEr:50];
		end
		Er_sim_pts = C.Er;
		C = getdRdErCENNS(C); % cts/keV/kg/day
		dRdEr = C.dRdEr;
		if ~exist('rate','var');rate=1;end
		switch rate
		case 0
			targetMass = 118; % target mass, kg
			liveDays = 85.3; % yes, days
		case 1 % 1 LZ = 15.34 tonne years
			targetMass = 5600; % target mass, kg
			liveDays = 1000; % yes, days
		case 2 % 10 tonne years
			targetMass = 1000; % target mass, kg
			liveDays = 365*10; % yes, days
		case 10 % 10 LZ = 153.4 tonne years
			targetMass = 5600; % target mass, kg
			liveDays = 1000*10; % yes, days
		case 12 % 1e2 LZ = 15.34*1e2 tonne years
			targetMass = 5600; % target mass, kg
			liveDays = 1000*100; % yes, days
		case 13 % 1e4 LZ = 15.34*1e3 tonne years
			targetMass = 5600; % target mass, kg
			liveDays = 1000*1e3; % yes, days
		case 14 % 1e4 LZ = 15.34*1e4 tonne years
			targetMass = 5600; % target mass, kg
			liveDays = 1000*1e4; % yes, days
		end
		Nstats0 = round(dRdEr * C.dEr * targetMass * liveDays); % cts
		
		%figure(501);ERcut=ER>=1;plot(S1(ERcut),log10(S2(ERcut)./S1(ERcut)),'ro')
	case 1000
		Er_sim_pts = [2:0.25:20]; 
		addpath(genpath('../limitCode/')); % assumes you are working in bandSimulation dir
		plotDMspectrum;
		Nstats0 = round(cts);
		
	otherwise	%case {1 2 3 4 5 6 8 10 12 16 20 24 32} %
		Er_sim_pts = [typ]; 
		Nstats0 = 5000;
	end	
	fs = 16;

%%% detector params
for iii=1
	
	alpha1 = 0.075; 
	nPMTtop = 247;
	nPMTbot = 241;
	PMTsig = 0.5*ones(1,nPMTtop+nPMTbot);
	alpha2 = 50; % phe / e-
	alpha2Sigma = 1/3;
	nco = 3;
	minpheArea = 1/3;
	dt_window = 100; % ns
	eta_extraction = 0.95;
	
	alpha_t = ones(nPMTtop,1) * 0.20 / nPMTtop; 
	alpha_b = ones(nPMTbot,1) * 0.80 / nPMTbot;
	alpha = [alpha_t ; alpha_b];
	modifyAcceptance = 0;
	suffix='LZ';
	
	%% physical params
		tau1 = 27; % ns
		%tau1 = 45; % ns - elongate to account for shaping?
			
	%% define the model prediction
		Z = 54;
		A = 131.3;
		wq = 0.0138;
		epsilon = 11.5 * Er_sim_pts / Z^(7/3);
		%k = 0.133 * Z^(2/3)*A^(-1/2);
		g = 3*epsilon.^0.15 + 0.7*epsilon.^0.6 + epsilon;
	
		for ii=1 % use 1 (default)
			if ii==1 % model case A
				k = 0.110;
				NexOverNi = 1;
				alphaOveraSquaredv = 0.037;
			elseif ii==1.05
				k = 0.110;
				NexOverNi = ii;
				alphaOveraSquaredv = 0.042;
			elseif ii==1.09 % XENON10 S2-only paper
				k = 0.110;
				NexOverNi = ii;
				alphaOveraSquaredv = 0.032;
			elseif ii==2 % NEST prediction for XENON10 730 V/cm
				k = 0.110;
				NexOverNi = 0.773;
				alphaOveraSquaredv = 0.0378;			
			elseif ii==15 % 2015 LUX prelim approx
				k = 0.166;
				NexOverNi = 0.5;
				alphaOveraSquaredv = 0.04;				
			elseif ii==-1
				k = kay; % somewhere, k got redefined, and I am too tired to find that shit
			end			
			fn = k .* g ./ (1 + k.*g);
			Ntot = Er_sim_pts .* fn / wq;
			Ni = Ntot./(1+NexOverNi);
			%Nex = Ntot./(1+1./NexOverNi);
			xi = Ni/4 * alphaOveraSquaredv;
			Fe = 1./xi.*log(1+xi).*Ni./(Ntot);
			modelQy = 1./xi .* log(1+xi) .* (fn/wq) ./ (1+NexOverNi);
			switch 1
			case 0 % old skool
				modelLeff = (1-(1./xi .* log(1+xi).*Ni./Ntot)).*(fn/wq).*(alpha1/Ly).*(S(1)/S(2));
			case 1
				modelLeff = (1-(1./xi .* log(1+xi).*Ni./Ntot)).*(fn/wq).*0.015;
			end
		end
		Qy = modelQy;
		Leff = modelLeff;

	figure(500);clf;
	topPlot=axes('Position',[0.15    0.55    0.80    0.375]);% [left, bottom, width, height] %subplot(2,1,1);
		plotQypts;
		h=plot([3 5 10 20],[4.5 4.6 4.51 4.07],'m^');
		h=plot(Er_sim_pts,modelQy,'k--','linew',1.);set(h,'Color',[0.9 0.6 0.1]);
		h=plot(Er_sim_pts,Qy,'k-','linew',1.);
		
		axis([0.9 80 2.5 11]);ax=axis;
		xlabel('');ylabel('');
		h=text(ax(1)-0.3,(ax(4)+ax(3))/2,'$\mathcal{Q}_{y}$','rotation',90,'horizontalAlignment','center');
		set(gca,'xsc','log','ysc','lin');
		set(gca,'xtick',[2:1:10 20:10:100],'xticklabel',[],'ytick',[0:1:16],'yticklabel',[]);
        botPlot=axes('Position',[0.15    0.15    0.80    0.375]);%subplot(2,1,2);
		plotLeffpts;
		h=plot(Er_sim_pts,modelLeff,'k--','linew',1.);set(h,'Color',[0.9 0.6 0.1]);	
		h=plot(Er_sim_pts,Leff,'k-','linew',1.);

		plot([5 10 20],[0.0821 0.0996 0.1248],'m^');
		axis([ax(1:2) 0 0.22]); ax=axis;
		xlabel('');ylabel('');
		h=text(10,ax(3)-0.06,'E$_{nr}$ [keV]','horizontalAlignment','center');
		h=text(ax(1)-0.42,(ax(4)+ax(3))/2,'$\mathcal{L}_{eff}$','rotation',90,'horizontalAlignment','center');
		set(gca,'xsc','log','ysc','lin');
		set(gca,'xtick',[1:1:10 20:10:100],'xticklabel',[],'ytick',[0:0.025:0.30],'yticklabel',[]);		

		set(gca,'Color','white');
	drawnow;

		Ne = 	zeros(sum(Nstats0),1);
		Ng = 	zeros(sum(Nstats0),1);
		S1 = 	zeros(sum(Nstats0),1);
		S2 = 	zeros(sum(Nstats0),1);
		tmpNe = zeros(sum(Nstats0),1);
		ER =  	zeros(sum(Nstats0),1);
		S1all = zeros(sum(Nstats0),1);

	Nlast=1;
	for ii=1:length(Er_sim_pts)
		Er = Er_sim_pts(ii);
		Nstats = Nstats0(ii); %round(Nstats0*exp(-Er/20));
		if Nstats>0 %Er==5
			R = [Nlast:Nlast+Nstats-1]; Nlast = R(end);
			ER(R) = Er*ones(size(R,1),size(R,2));
			%% generate a statistical sample
			Nt0 = round(Ntot(ii)); % Ntot is a mean value from the model
		    Ne(R) = binornd(Nt0*ones(1,Nstats),Fe(ii))'; % Nt0 is split into electrons and photons
			S2(R) = alpha2* (sqrt(Ne(R)).*alpha2Sigma.*randn(Nstats,1)+Ne(R));
			% photons
			Ng(R) = (Nt0-Ne(R)); % number fluctuations already accounted for
			hitpatterns = binornd(repmat(Ng(R),1,size(alpha,1)),alpha1*repmat(alpha',Nstats,1)); % distribute Ng photons onto the PMTs, with specified probability to create a phe
			phe = (hitpatterns>0).*randn(size(hitpatterns,1),size(hitpatterns,2)).*repmat(PMTsig,Nstats,1)+hitpatterns; % phe distributions are Gaussian with spcified sigma/mu
			phe(phe<0) = 0; % physicality req.
			for k=1:size(phe,1)	% did we keep the event? -- check coincidence timing
				S1all(R(k)) = sum(phe(k,:).*(phe(k,:)>minpheArea),2);
				nphe = sum(phe(k,:)>minpheArea);
				if nphe >= nco
					tp = sort( exprnd(tau1*ones(nphe,1) ) );
					dt_diff = ( tp(nco:end) - tp(1:end-nco+1) );
					if min(dt_diff) <= dt_window 
						S1(R(k)) = sum(phe(k,:).*((phe(k,:)>minpheArea)),2);
					end
				end
            end
			% last, reduce Ne for extraction efficiency
			if exist('eta_extraction','var') % extraction efficiency
				Ne(R) = binornd(Ne(R),eta_extraction);
				S2(R) = alpha2* (sqrt(Ne(R)).*alpha2Sigma.*randn(Nstats,1)+Ne(R));
            end
        end
    end
end % iii
return

n=20;
p=0.5;
bins=[0:1:20];
figure(11);
logstep( binornd(n*ones(1e3,1),p),bins,'b');
hold on;
logstep(poissrnd(n*p*ones(1e3,1)),bins,'k');
axis([0 30 1 1e3]);
