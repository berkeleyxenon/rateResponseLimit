function C = getHelmFF(C)
% function C = getHelmFF(C)
%
% 090715 pfs adapted from previous code


C.q = sqrt( 2 * C.M_N * C.Er ) ./ (C.hbar*C.c*1e3/100*1e15); % M_A in GeV, E_R in keV, 197.3 in MeV*fm
switch C.Z %experiment
%case {'XENON10' 'XENON100' 'DAMA' 'LUX' 'CRESST-II'} % Helm, following 1996 Lewin & Smith
case {14 18 23 32 53 54 74} % Helm, following 1996 Lewin & Smith
	s = 1; % 1 fm -- skin depth of nucleus
	%r = 1.2*C.M_N^(1/3);  % think this is wrong
	r = 1.2*C.A^(1/3);  % as in 2009 Chang et al. (??)
	r_0 = sqrt(r^2 - 5*s^2);
case 'NOT_IN_USE'	% is this the 2-parameter Fermi ?
	%r = 0.91*C.M_N^(1/3) + 0.3;  % don't use this, it is "less accurate"
	s = 0.9;
	a = 0.52;
	%r = 1.23*C.M_N^(1/3); - 0.60;
	r = 1.23*C.A^(1/3); - 0.60;
	r_0 = sqrt( r^2 + 7/3*pi^2*a^2 - 5*s^2 );  % 2-param Fermi model (aka Woods-Saxon?) from hep-ph/0608035
end

qR = C.q.*r_0;
C.qR = qR;
% ------------------
F = (3./qR) .* ( sin(qR)-qR.*cos(qR) )./qR.^2  .* exp(-(C.q.^2.*s.^2/2)); % factor 2 !
C.F = F;
if 0
	figure(111); clf;
		subplot(2,1,1)
		semilogy( C.q , C.F.^2 , ['ro']); hold on;
		xlabel('qr')
		ylabel('F^2');
		axis([0 1 1e-3 1]);
		grid on;
		subplot(2,1,2);
		semilogy( C.Er , C.F.^2 , ['bo']); hold on;
		xlabel('Er')
		ylabel('F^2');
		axis([0 100 1e-3 1]);
		grid on;
		pause
end


	