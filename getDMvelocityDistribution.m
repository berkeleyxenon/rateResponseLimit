function eta = getDMvelocityDistribution(C,t);
%
% returns halo DM velocity distribution based on REF1: PRD 74 043531 (2006)
%
% 090714 pfs (Bastille Day! ooo-la!)
%
% NOTE: XENON10 used
% v0 = 220 km/s
% v_esc = 650 km/s
% rho_d = 0.3 GeV/cm^3
% Helm FF with s = 1.0; c = 1.2 * A**(1/3.);rn = sqrt(c*c - 5 * (s*s));
%

v_esc = C.v_esc*3e5; % take this as input, 500 km/s default
v0 = C.v_0*3e5; % local standard of rest, taken to be 220 km/s
%%% define v_min (REF1 eq 3)
%v_min = sqrt(M.*Er./2./mu.^2); % minimum WIMP velocity that can result in a recoil energy Er
% I need to use beta_min
v_min = C.beta_min*3e5;

%%% Standard Halo Model
% define v_obs(t) from REF1 eq 18
% this involves assumption that v_earth << v_halo, i.e. 30 << 200, which seems somewhat OK
t_c = 0.415; % june 1
switch 1 % USE 1
case 0 % simplest assumption
	v_sun = v0+12;
	b = 0.51;
case 1 % following REF1
	v_sun = sqrt(10^2+(13+v0)^2+7^2); % km/s
	b = 0.49; % this is cos(gamma), giving the inclination of the ecliptic
case 2 % following 1999 Dehnen & Binney
	v_sun = sqrt(10^2+(5.23+v0)^2+7.17^2); % km/s
	b = 0.49; % this is cos(gamma), giving the inclination of the ecliptic
end
%figure(1);t=0:0.001:1;plot(t,cos(2*pi*(t-t_c)),'r-')

v_earth 	= 29.8; % km/s
v_obs 		= v_sun .* (1 + b.*v_earth./v_sun.*cos(2*pi*(t-t_c))); % REF1 eq 18

%v_e 		= v0*( 1.05-0.07.*cos(2*pi*(t-t_c)) ); %
v_e = v_obs;

% example modulation with May 1 and Aug 1 approximately indicated:
%figure(1);clf;t=0:0.001:1;plot(t,v_obs);xlabel('time of year [fraction]');ylabel('$v_0$ [km/s]');hold on;plot([1 1]*4/12,[210 250],'r--');plot([1 1]*7/12,[210 250],'r--')

%%% NOTE: up to this point I am working in km, but I need to use cm for the normalization since sigma is measured in cm^2 (in other words need a factor x1e5)
conversh = 1e5;

if ~isfield(C,'velocityIntegrationMethod'); 
	C.velocityIntegrationMethod = 'PRD_74_043531_2006'; 
end;
switch C.velocityIntegrationMethod
case 'JKG_w_escape' % from JKG p122 eq8.16
	col='k';
		TQ = sqrt(pi)/4*(v0/v_e)*( erf( (v_min+v_e)/v0 ) - erf( (v_min-v_e)/v0 ) ) ...
		-exp(-v_esc.^2./v0.^2) ...
		;
		eta = (TQ*2)/(conversh*sqrt(pi)*v0);
		eta(eta<0)=0; 
case 'explicit_JKG_w_escape' % integrate it mah sef
	col='b';
	dv = 1;
	v_inf=1000; % infinity, numerically ;p
	switch 2 % use 2
	case 0 % slow matlab loop
		for i=1:length(v_min)
			v = v_min(i):dv:v_inf; % NOTE: limits of integration are infinity.  excape velocity is accounted for separately.
			f1 = (1/sqrt(pi))*v./(v_e*v0) .* ( exp(-(v-v_e).^2./v0^2) - exp(-(v+v_e).^2./v0^2) );
			TQ(i) = sqrt(pi)/2 .* v0 .* sum(f1./v.*dv) - exp(-v_esc.^2./v0.^2);
		end
	case 1 % slow vector coding (memory allocation burns it)
		v=repmat([dv:dv:vinf]',1,size(v_min,2)); 
		v=v.*(v>repmat(v_min,size(v,1),1)); v(v==0)=eps;
		f1 = (1/sqrt(pi))*v./(v_e*v0) .* ( exp(-(v-v_e).^2./v0^2) - exp(-(v+v_e).^2./v0^2) );
		TQ = sqrt(pi)/2 .* v0 .* sum(f1./v.*dv) - exp(-v_esc.^2./v0.^2);
	case 2 % mex loop - word.
		if C.accountForInelasticNucleusScattering
			elasticFraction = C.elasticFraction;
			%keV = elasticFraction.eV/1e3; % initial particle energy assuming m_n = 1 GeV
			keV = elasticFraction.eV/1e3 * (1+C.m_chi/elasticFraction.(C.element).amu);  % initial particle energy assuming m_chi
			f = 1 - (elasticFraction.(C.element).Inelastic ./ elasticFraction.(C.element).Total)';
		else
			keV = -1;
			f = -1;
		end
		TQ = doDMvelocityIntegrationLoop(v_min,dv,v_inf,v_e,v_esc,v0,C.m_chi,keV,f);		
		%pause;
	end
	
	eta = (TQ*2)/(conversh*sqrt(pi)*v0);
	eta(eta<0)=0; 
case 'PRD_74_043531_2006'
	col='r';
	%%% define v0 == v_0_bar (REF1 eq 9)
	% here we assume that v0 == sigma_v, i.e. velocity dispersion == velocity(local standard of rest)
	% and v0 == "most probable speed" (this is just v_0 in PRD 72 063509 (2005) )
	% so basically people say "v_0" or "sigma_v" or "v0" interchangeably, but they are different 
	% by a rather important sqrt(2/3)
	% i.e.  v0 = 220.  Not 270, not 180.

	%%% define dimensionless variables (REF1 eq 11)
	x = v_min ./ v0;
	y = v_obs ./ v0;
	z = v_esc ./ v0;

	%%% define dimensionless normalization (REF1 eq 8)
	N_esc = erf(z) - 2.*z.*exp(-z.^2)./sqrt(pi);
	
	%%% define boolean for evaluating eta
		logi(1,:) = z<y & x<abs(y-z);
		logi(2,:) = z>y & x<abs(y-z);
		logi(3,:) = x>abs(y-z) & x<(y+z);
		logi(4,:) = x>(y+z);
	
	%%% define eta(E,t)
		eta(1,:) = 1./(v0*conversh.*y) .* ones(1,size(v_min,2));
		eta(2,:) = 1./( 2.*N_esc.*v0*conversh.*y) .* (erf(x+y) - erf(x-y) - 4/sqrt(pi).*y.*exp(-z.^2) );
		eta(3,:) = 1./( 2.*N_esc.*v0*conversh.*y) .* (erf(z) - erf(x-y) - 2/sqrt(pi).*(y+z-x).*exp(-z.^2) );
		eta(4,:) = 0 .* ones(1,size(v_min,2));

	if 0 % check kinetic energy ?  relevant to population of inelastic states
		keyboard		
		T = sqrt(0.5*(v_min/3e5).^2*C.m_chi*1e9);
		figure(2);clf;
			loglog(C.Er,T);grid on
			xlabel('recoil energy     [keV]');ylabel('DM kinetic energy     [keV]');myFigView;
	end
	
	eta = max(eta .* logi);
	eta(eta<0)=0; % don't know why it goes negative...

end

if 0 % plot
	v=dv:dv:1000;
	T = 0.5*C.m_chi*1e9*(v/3e5).^2 /1e3; % last 1e3 converts eV=>keV
	T1 = 0.5*50*1e9*(v/3e5).^2 /1e3; % last 1e3 converts eV=>keV
	T2 = 0.5*500*1e9*(v/3e5).^2 /1e3; % last 1e3 converts eV=>keV
	f1 = (1/sqrt(pi))*v./(v_e*v0) .* ( exp(-(v-v_e).^2./v0^2) - exp(-(v+v_e).^2./v0^2) );
	figure(2);clf;
		plot(v,f1./v,'k-');
		hold on;
		%plot(T,f1./v,'b-');
		plot(T1,f1./v,'b--');
		plot(T2,f1./v,'r-.');
		legend('km s^{-1}' ...
				,'keV, m_{\chi}=50 GeV' ... ()
				,'keV, m_{\chi}=500 GeV' ... ()				
				);
		xlabel('km s^{-1} or keV     (see legend)');
		ylabel('Velocity Distribution f_1(v)/v');
		myFigView(18);
		%title(['m_{\chi} = ' dis('%1.0f GeV',C.m_chi)]);
		if 1
			set(gca,'xsc','log');
			axis([20 5000 0 1.1e-5]);
			set(gca,'xtickl',[100 1000],'ytickl',[0:0.2:1]);
		end
	pause
	%keyboard
	%sav(2,'fofv_vs_km_and_E')
end

if 0 % debug
	figure(2); hold on;
		loglog(C.Er,eta,col);grid on
		xlabel('recoil energy     [keV]');ylabel('eta');%myFigView;
	%keyboard
end
