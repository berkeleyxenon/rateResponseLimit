% code for plotting output of NRbandsim
%
% updated 120706 pfs

%% load or define external bands data
	load([pathh '/dataFiles/xebands_S1_090903.mat']);			
	
	load([pathh '/dataFiles/xenon100_centroids.mat']); % loads variable 'xenon100'
	xenon100.s1a = [1 2 3 4 5 6 7 8 9 10 11 12 13 14];
	xenon100.s1acc = [0.25 0.46 0.65 0.78 0.86 0.9 0.94 0.96 0.97 0.98 0.985 0.99 0.995 1.00];
    yy = log10(S2./S1);

	clear mu sig ncts;
	logBins = [1:0.05:4];
	binEdges=[1 2 3 4 5 6 7 8 9 10 12 14 16 18 20 25:5:50 60:10:100];
	minS2phe = 250; % 5 e-
	binCenters = (binEdges(2:end)+binEdges(1:end-1))/2;
	kut = ER>=0 & S2>minS2phe;
	kut = S1>0 & S2>minS2phe;
	
%% fit the NR band in the same manner that the data was treated
if typ==0.1 | typ==0.2
	disp('** fitting bands and storing band parameters **');
	mu=[]; sig=[]; Nl=zeros(1,size(Er_sim_pts,2));
	cuts.box = false(size(S1,1),size(S1,2)); % NR 50% box
	% loop over S1 bins
	for k=1:length(binEdges)-1 
		cuts.S1 = S1>=binEdges(k) & S1<binEdges(k+1) & kut;
		ncts(k) = sum(cuts.S1);
		N = hist(yy(cuts.S1),logBins);
		if sum(N)>20
			switch 1
			case 0
				[r,c]=max(N);
				fitLowerLim = logBins(c-2);			
				while N(c)>2%(r/50)
					c=c+1;
				end
				fitUpperLim = logBins(c);		
				try
					[A,Y,chi,J] = FH_Gauss(logBins,N,fitLowerLim,fitUpperLim,[]);
				catch
					A = [0 0 0];
				end
				mu(k) = A(2);
				sig(k) = A(3);
	
				if 1
				figure(5000);clf; axes('Position',[0.15 0.15 0.80 0.70]);
					erb(logBins,N,sqrt(N),'kb');
					hold on;
					plot(logBins,gaussian(A,logBins),'ro');
					plot([1:0.01:3.3],gaussian(A,[1:0.01:3.3]),'r-');
					plot(fitUpperLim*[1 1],[1 1e3],'r-');
					plot(fitLowerLim*[1 1],[1 1e3],'r-');
	
					set(gca,'ysc','log');
					set(gca,'xtick',[1:0.1:4],'xticklabel',{'1' '' '' '' '' '1.5' '' '' '' '' '2' '' '' '' '' '2.5' '' '' '' '' '3' '' '' '' '' '3.5' '' '' '' '' '4'});
					set(gca,'ytick',10.^[0:1:4]);
					axis([1 3.5 0.8 5e4]);
					%title(dis('%d-%d S1 photoelectrons',binEdges(k),binEdges(k+1)));
					grid on;
					drawnow;
					pause;
				% calculate lower limit of box from wmean (more inclusive)
				end
				
			case 1
				[mm ss] = wmean(N,logBins);
				mu(k) = mm;
				sig(k) = ss;				
			end
			cuts.box = cuts.box | (cuts.S1 & yy<mu(k));
	
		else
			%dis('insufficient stats, bailing on fit/mean calc...')
			mu(k) = 0;
			sig(k) = 0;
		end
		
		%dis('done with bins %d-%d (of %d)',binEdges(k),binEdges(k+1),binEdges(end));
	end
	wb = 1:13; % which bins
	%chisq_pieces = (mu(wb) - xebands.NR.mu(wb)).^2 ./ (sig(wb)./sqrt(ncts(wb))).^2;
	chisq_pieces = (mu(wb) - xebands.NR.mu(wb)).^2 ;
	chisq_mu = sum(chisq_pieces);
	
	% also add (a kludged) chisq_pieces for sigma: 
	tmp = sig(wb) - xebands.NR.sig(wb); tmp(tmp<0)=0; % can't have sigma_sim > sigma_data !!
	chisq_pieces_sig = tmp.^2 .* ncts(wb) * 100; % pulled factor x100 from thin air
	chisq_sig = sum(chisq_pieces_sig);
	 
	
	xi = [1:1:40];
	%xi = [1:1:75];
	plot(binCenters,mu,'go')
	yi_mu = interp1(binCenters,mu,xi,'pchip');
	yi_m3s = interp1(binCenters,mu-3*sig,xi,'pchip');
	
	NRsimBands.be = binEdges;
	NRsimBands.bc = (binEdges(1:end-1)+binEdges(2:end))/2;
	NRsimBands.m = mu;
	NRsimBands.s = sig;
	save('./dataFiles/NRsimBands.mat','NRsimBands');
else%if typ==1000 | typ==0
	load('./dataFiles/NRsimBands.mat');
	mu = interp1(NRsimBands.bc,NRsimBands.m,S1,'pchip');%plot(S1,mu,'r.')
	m3s = interp1(NRsimBands.bc,(NRsimBands.m-3*NRsimBands.s),S1,'pchip');%plot(S1,m3s,'r.')
end

%% plots
hleg=[];leg=[];
figure(501);clf;
hold on;
h=plot( S1(kut) , yy(kut) , 'bo','markers',2);
axis([0 25 1.5 3.5]);ax=axis;
h=plot([1:1:20],log10(minS2phe./[1:1:20]),'k--');
set(gca,'layer','top');
set(gca,'xminortick','on','yminortick','on');
xx=[0.1:0.1:100];
plot(xx,interp1(NRsimBands.bc,(NRsimBands.m+3*NRsimBands.s),xx,'pchip'),'r-','linew',0.5);
plot(xx,interp1(NRsimBands.bc,NRsimBands.m,xx,'pchip'),'r-','linew',2);
plot(xx,interp1(NRsimBands.bc,(NRsimBands.m-3*NRsimBands.s),xx,'pchip'),'r-','linew',0.5);
total_events = sum(S1>0&S2>200);
[S1 S2];
xlabel('S1 photoelectrons');
ylabel('log$_{10}$(S2/S1)');
set(gcf,'renderer','painters');
hold off;

S1_out = S1(kut);
logS2S1_out = yy(kut);
