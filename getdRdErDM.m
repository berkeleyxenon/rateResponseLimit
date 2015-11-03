function dR_dEr = getdRdErDM(C);
% function dR_dEr = getdRdErDM(C);
%
% 090715 pfs
% 130306 pfs - output in cts/keV/kg/day (rather than sec)

%%
C.M_N = C.A*C.mp; % GeV
C = getHelmFF(C);
C = getBetaMin(C);%%% calculate beta_min

%% do the velocity integration
	C.eta = getDMvelocityDistribution(C,C.t);

%% calculate WIMP-nucleon reduced mass
	C.mu_n = (C.m_chi.*C.mp) ./ (C.m_chi+C.mp); % reduced mass of WIMP-nucleon system (GeV)


%% calculate expected event rate -- folowing PRD 79 043513 (2009)
dR_dEr = 0.5 * C.N_0/C.A*1e3 .* C.M_N ... kg^-1
	.* C.rho_chi./C.m_chi ...
	.* C.sigma_n./C.mu_n.^2 ...
	.* (C.f_p.*C.Z + C.f_n.*(C.A-C.Z)).^2 ./ C.f_n.^2 ... % A^2 coherence term
	.* C.F.^2 ... % form factor
	.* C.eta ... % integral of velocity distribution function
	.* C.c.^2 ... % adds in the missing factor c^2
	./ 1e6 ... % convert GeV => keV
	.* 86400 ... % sec => day
	;

