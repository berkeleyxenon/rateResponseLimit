% calculate pp neutrino rate in LZ, mostly following 
%   bahcall 1987 RMP 
%   PDG 2013 
%   also consulting 1309.7024, 1307.5081
%   http://www.sns.ias.edu/~jnb/SNdata/sndata.html
%
% 150721 pfs
% 150806 pfs adding atomic constraints approximation
% 150831 pfs add 7Be and automate for various solar models


sigma0 = 88.083e-46; % cm^2
%sigma_enu = 11.6e-46; % cm^2 ... not used
Gf = 1.166e-5; % GeV^-2 (this is Gf/(hbarc)^3)
C.hbar = 6.5822*1e-25; % in GeV * s
C.c = 2.9979e10; % cm/s    sin2ThetaW = 0.231;
dq = 1; % keV
me = 511; % keV
F = 1; % "near-unity correction factor" due to Weinberg angle uncertainty, not implemented

for mm=1%:7
	phi=getSolarNuFlux(mm);
	% which interaction?
	for nn=[1 2 3 4]
		clear dN_dq dsigma_dT0 dsigma_dT0_mutau dNe_by_dT0 dNe_by_dT0_mutau;
		switch nn
		case 1
			disColor = 'b';
			Q = 420.3; % pp endpoint, keV
			q = dq:dq:Q;
			% pp neutrino differential energy spectrum
			% Phys. Rev. C. 56 3391(1997)
			% http://www.sns.ias.edu/~jnb/SNdata/ppspectrum.html
			fid=fopen('./dataFiles/ppenergytab.txt');
			A=fscanf(fid,'%s',8); % header
			count=0;
			for z=1:21
				for i=1:4
					count=count+1;
					b.x(count) = str2num(fscanf(fid,'%s',1));
					b.y(count) = str2num(fscanf(fid,'%s',1));
				end
			end
			fclose(fid);
	
			dN_dq = (Q+me-q) .* sqrt( (Q+me-q).^2-me^2) .* q.^2 * F;
			dN_dq = dN_dq/sum(dN_dq)*phi.pp; % normalize to total flux

			if 0 % verify that parameterization from 1309.7024 matches Bahcall data (it does)
			%figure(1);plot(q,dN_dq);xlabel('neutrino energy / keV');ylabel('flux');
			figure(1);clf;
				loglog(b.x,b.y*phi.pp,'ro'); hold on;
				loglog(q/1e3,dN_dq*1e3/dq);
			end
		case 2
			disColor = 'g';
			Q = 862; % 7Be line, keV
			q = dq:dq:Q;
			dN_dq = zeros(size(q));
			dN_dq(end) = 1;
			dN_dq = dN_dq/sum(dN_dq)*phi.Be7*0.897; % normalize to total flux
		case 3
			disColor = 'r';
			Q = 384; % 7Be line, keV
			q = dq:dq:Q;
			dN_dq = zeros(size(q));
			dN_dq(end) = 1;
			dN_dq = dN_dq/sum(dN_dq)*phi.Be7*0.103; % normalize to total flux
		case 4
			disColor = 'm';
			Q = 1199; % 13N endpoint, keV
			q = dq:dq:Q;
			% CNO neutrino differential energy spectrum
			% J. N. Bahcall and R. K. Ulrich, Rev. Mod. Phys. 60, 297 (1988)
			fid=fopen('./dataFiles/n13.dat.txt');
			count=0;
			for z=1:100
				for i=1:2
					count=count+1;
					b.x(count) = str2num(fscanf(fid,'%s',1));
					b.y(count) = str2num(fscanf(fid,'%s',1));
				end
			end
			fclose(fid);
			dN_dq = (Q+me-q) .* sqrt( (Q+me-q).^2-me^2) .* q.^2 * F;
			dN_dq = dN_dq/sum(dN_dq)*phi.N13; % normalize to total flux
			if 0 % verify that parameterization from 1309.7024 matches Bahcall data (it does)
			figure(1);clf;
				loglog(b.x,b.y*phi.N13,'ro'); hold on;
				loglog(q/1e3,dN_dq*1e3/dq);
			end
			
		end

		T = q; % for convenience in calculating
		Tmax = 2*(q/me).^2./(1+2*q/me) * me; % check
		qmin = (T/me+sqrt(T/me.*(T/me+2)))/2 * me; % enforced by Tmax so not needed explicitly
		bahcallsigma0 = 2*Gf^2*me^2/pi .* (C.hbar*C.c)^2 / (1e6)^2;


	
		sin2ThetaW = 0.231;
		gL = 0.5+sin2ThetaW;
		gR = sin2ThetaW;

		for jj=1:length(q)
			for ii = 1:length(q)
				% bahcall
				dsigma_dT0(ii,jj) = sigma0/me*(gL^2+gR^2*(1-T(ii)/q(jj))^2-gL*gR*me*T(ii)/q(jj)^2) ... % *me
					.* (T(ii)<Tmax(jj)) ...
					....* (q(jj)>qmin(ii)) ...
					;
				dsigma_dT0_mutau(ii,jj) = sigma0/me*((sin2ThetaW-0.5)^2+gR^2*(1-T(ii)/q(jj))^2-(sin2ThetaW-0.5)*gR*me*T(ii)/q(jj)^2) ... % *me
					.* (T(ii)<Tmax(jj)) ...
					....* (q(jj)>qmin(ii)) ...
					;
		
				if 0
				% 1309.7024
				dsigma_dT1(ii,jj) = 2*Gf^2*me/pi .* (C.hbar*C.c)^2 / (1e6)^2 ...
					.* (gL^2 + gR^2.*(1-T(ii)./q(jj)).^2 - gL*gR.*me./T(ii)./q(jj)) ...
					.* (T(ii)<Tmax(jj)) ...
					;
				% 1307.5081
				dsigma_dT2(ii,jj) = 2*Gf^2*me/pi .* (C.hbar*C.c)^2 / (1e6)^2 ...
					.* (gR^2 + gL^2.*(1-T(ii)./q(jj)).^2 - gL*gR.*me.*T(ii)./q(jj)^2) ...
					.* (T(ii)<Tmax(jj)) ...
					;
				end

			end
		end

		%figure(2);semilogy(T,dsigma_dT)


		t = 60*60*24*365; % convert s => year 
		switch 0
		case 0 % free electron approximation
			No = ones(size(dN_dq))*54*6.02e23/131.3*1e6; % electrons in target, assume 1 tonne target
		case 1 % simple atomic binding treatment -- WRONG !!!
			% NOTE: using DQ shell energies from LUXDB348, these are different from xdb.lbl.gov, should probably double check
			for ii=1:length(T)
			if T(ii)<=0.19
				No(ii) = (54-2-8-18-32)*6.02e23/131.3*1e6; % electrons in target, assume 1 tonne target
			elseif T(ii)<=1.1
				No(ii) = (54-2-8-18)*6.02e23/131.3*1e6; % electrons in target, assume 1 tonne target
			elseif T(ii)<=5.2
				No(ii) = (54-2-8)*6.02e23/131.3*1e6; % electrons in target, assume 1 tonne target
			elseif T(ii)<=33.2
				No(ii) = (54-2)*6.02e23/131.3*1e6; % electrons in target, assume 1 tonne target
			else
				No(ii) = 54*6.02e23/131.3*1e6; % electrons in target, assume 1 tonne target
			end
			end
		end

		%test=100;qmin(test)
		%[q' dsigma_dT(test,:)'*1e46]
		for ii=1:length(T)
			dNe_by_dT0(ii) = No(ii) * t * sum( dN_dq .* dsigma_dT0(ii,:) );
			dNe_by_dT0_mutau(ii) = No(ii) * t * sum( dN_dq .* dsigma_dT0_mutau(ii,:) );
			P(ii) =  sum( dN_dq .* dsigma_dT0(ii,:) );
	
			%dNe_by_dT1(ii) = No * t * sum( dN_dq .* dsigma_dT1(ii,:) );
			%dNe_by_dT2(ii) = No * t * sum( dN_dq .* dsigma_dT2(ii,:) );

		end

		% just to compare shape:
		P = P/sum(P)*1000; % basically agrees with Bahcall Table VIII and Fig. 5
		%[T' P']

		%% neutrino mixing from sun propagates as if in vacuum (see PDG review)
			% survival probability for nu_e
			sin2theta12 = 0.307; % +0.018-0.016
			sin2theta13 = 0.024-0.0025; % +-0.0025

			theta12 = asin(sqrt(sin2theta12)); % 33.6 degrees
			theta13 = asin(sqrt(sin2theta13)); % 8.9 degrees


			Pnue = sin(theta13)^4 + (1-0.5*sin(2*theta12)^2)*cos(theta13)^4;

		%% primary result
		dNe_by_dT0_total = dNe_by_dT0*Pnue + dNe_by_dT0_mutau*(1-Pnue);
		% interpolate to evade binning effects (just for headline pp rate):
		dTi=1e-3;
		Tinterp = dTi:dTi:Q;
		dppinterp = interp1(T,dNe_by_dT0_total,Tinterp,'pchip',0);

		%cts=sum(dNe_by_dT0_total(q>=1.5 & q<=6.5)) * dq * 2.73 * 5.6
		ctsinterp=sum(dppinterp(Tinterp>=1.5 & Tinterp<=6.5)) * dTi * 2.73 * 5.6;
		
		% the answer!
		flux(mm,nn) = ctsinterp;


		if 1 % plot fluxes
		figure(4); hold on;
			h=plotstairs(T,dNe_by_dT0*Pnue,'k--'); set(h,'color',disColor);
			h=plotstairs(T,dNe_by_dT0_mutau*(1-Pnue),'k:');set(h,'color',disColor);
			h=plotstairs(T,dNe_by_dT0_total,'k-.');set(h,'color',disColor);
			h=plot(Tinterp,dppinterp,'k-');set(h,'color',disColor);

			xlabel('electron energy / keV')
			ylabel('cts/keV/tonne/year');
	
			%legend('$\nu_e$ $~$','$\nu_{\mu},\nu_{\tau}$','all $\nu$','$\nu_e$ 1309.7024 Fig. 1$~$','location','ne');
			plot([1.5 1.5],[1e-2 10],'k-');
			plot([6.5 6.5],[1e-2 10],'k-');
			%text(1.6,1,dis('LZ search \nwindow'),'color','red');
			title(['$\phi.{pp} = $ ' dis('%3.2f',phi.pp/1e10) '$\times10^{10}$ / cm$^2$/s ']);
			set(gca,'xsc','lin','ysc','lin');
			ax=axis;axis([0 300 1e-2 10]);
			%set(gca,'xtick',[0:50:300]);
			setplot;
			box on;
			axis([0 10 0 3.3]);
		end
	end
end

flux
tflux = sum(flux,2)
mean(tflux)
std(tflux)

if 0 % compare with Baudis
figure(3); clf;
	plotstairs(T,dNe_by_dT0*Pnue,'k--','linew',1);
	hold on;
	plotstairs(T,dNe_by_dT0_mutau*(1-Pnue),'k:','linew',1);
	plotstairs(T,dNe_by_dT0_total,'k-','linew',1);
	plot(Tinterp,dppinterp,'g--');

	%plot(T,dNe_by_dT1,'r-');
	%plot(T,dNe_by_dT1*0.54,'r--');
	%plot(T,dNe_by_dT2*0.54,'b-');
	xlabel('electron energy / keV')
	ylabel('cts/keV/tonne/year');
	
	if 1 % plot 1309.7024
		b.x=[1 8.62 29.0 56.9 83.6 106 130 149 170 188 204 218 228 239 247 254 261 267 270 273];
		b.r=[2.675 2.567 2.330 1.981 1.675 1.453 1.175 0.933 0.719 0.542 0.395 0.279 0.192 0.121 0.077 0.048 0.027 0.0131 0.0071 0.0037];
		plot(b.x,b.r,'bo');
	end

	legend('$\nu_e$ $~$','$\nu_{\mu},\nu_{\tau}$','all $\nu$','$\nu_e$ 1309.7024 Fig. 1$~$','location','ne');
	plot([1.5 1.5],[1e-2 10],'r-');
	plot([6.5 6.5],[1e-2 10],'r-');
	text(1.6,1,dis('LZ search \nwindow'),'color','red');
	title(['$\phi.{pp} = $ ' dis('%3.2f',phi.pp/1e10) '$\times10^{10}$ / cm$^2$/s ']);
	set(gca,'xsc','lin','ysc','lin');
	ax=axis;axis([0 300 1e-2 10]);
	%set(gca,'xtick',[0:50:300]);
	setplot;
end

% check binning effects:
%dqs=reverse([0.2:0.1:1]);for qq=1:length(dqs);dq=dqs(qq);ppNu;ctss(qq)=cts;ctssinterp(qq)=ctsinterp;end;figure(99);plot(dqs,ctss,'bo');hold on;plot(dqs,ctssinterp,'go');xlabel('bin width / keV');ylabel('LZ pp cts');


