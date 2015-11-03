% summarize CENNS vs coincidence
%
% 150903 pfs

nco = [1 2 3 4 5];
col = {'ro' 'bx' 'gv'};
%col = {'rb' 'bb' 'gb'};

% events below NR mean

% 8B hep atm DSN
N(1,nco) = [0 45 7.0 1.1 0.11];
N(2,nco) = [0 0.97 0.24 0.07 0.008];
N(3,nco) = [0 0.14 0.125 0.104 0.096];

U(1,nco) = [0 0.67 0.26 0.1 0.03];
U(2,nco) = [0 0.031 0.0015 0.008 0.003];
U(3,nco) = [0 0.0037 0.0035 0.0032 0.0031];

figure(8);clf;hold on; box on;
	for k=1:3
		h=erb(nco,N(k,:),U(k,:),col{k},0);
		%plot(nco,N(k,:),col{k});
		%plot(nco,[N(k,:)-U(k,:) ; N(k,:)+U(k,:)],col{k}(1));

	end
set(gca,'ysc','log');
axis([1.5 5.5 5e-3 1e2]);
grid off;
set(gca,'ytick',10.^[-2:1:2]);
xlabel('PMT coincidence');
ylabel('LZ counts 5600 tonne-day');
title('CENNS events below NR mean');
