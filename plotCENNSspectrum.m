% quick set up and plot
%
% 130110 pfs


%% particulars
	C.A = 131.3;
	C.Z = 54;
	C.dEr = 0.01; % should be sufficient granularity...
	C.Er = C.dEr/2:C.dEr:100;


	
	
figure(1);clf;
figure(11);clf;
	
col = {'r-' 'b:' 'g--' 'k-.'};
for ii=4%1:4
	C.source = ii;
	C = getdRdErCENNS(C); % dru

	figure(1);%clf;
		plotlogstairs(C.Er,C.dRdEr*1e3*365,col{ii});
		hold on;
		set(gca,'xsc','log');
		%set(gca,'xtick',[0:1:10]);
		xlabel('recoil energy [keV]');
		ylabel('cts/keV/tonnne/year');
		axis([0.1 100 1e-5 1e4]);
		set(gca,'ytick',10.^[-6:1:4]);

	for rr=1:length(C.Er)
		C.integratedRate(rr) = sum(C.dRdEr(rr:end))*C.dEr;
	end
	figure(2);%clf;
		plotlogstairs(C.Er,C.integratedRate*1e3*365,col{ii});
		hold on;
		plotlogstairs(C.Er,C.integratedRate*5.6e3*1000,col{ii});
		set(gca,'xsc','log');
		%set(gca,'xtick',[0:1:10]);
		xlabel('recoil energy [keV]');
		%ylabel('dru');
		axis([0.1 5 0.1 1e3]);
		%set(gca,'ytick',10.^[-6:1:4]);

	% counts in LZ?
	sum( C.dRdEr(C.Er>6 & C.Er<30) ) * C.dEr * 5.6e3 * 1000	
end
	
% how many events?
% sum(C.dRdEr(C.Er>3))*C.dEr*118*88.3
% 
% sum(C.dRdEr(C.Er>2))*C.dEr*118*88.3
% ans =
%     0.1308
% 
% sum(C.dRdEr(C.Er>1))*C.dEr*118*88.3
% ans =
%     2.7030
% 
% sum(C.dRdEr(C.Er>0))*C.dEr*118*88.3
% ans =
%   27.6222

