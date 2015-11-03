function C = getBetaMin(C);
% function C = getBetaMin(C);
%
% 090715 pfs

C.mu_wn = (C.m_chi.*C.M_N) ./ (C.m_chi+C.M_N); % reduced mass of WIMP-nucleus system (GeV)

C.beta_min = sqrt(1./(2*C.M_N.*C.Er*1e-6)) ...
		.* ((C.M_N.*C.Er*1e-6)./C.mu_wn ...
		+ C.delta*1e-6); % units of c		
