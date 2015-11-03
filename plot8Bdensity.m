% density plot can be made after calling NRbandsim
%
%

figure(502);clf;
	axes('Position',[0.15 0.15 0.80 0.70]);
	% set the axes
	h=plot(-1,0,'ko');hold on;
    % surf plot
    xbinz = [0:0.5:10]; xbinzc = (xbinz(1:end-1)+xbinz(2:end))/2;

    hsurf=mySurf(S1,log10(S2./S1),xbinz,[0:0.10:3.0]);	
    %shading interp;
	
	hold on;

	% minS2phe phe
	h=plot([0.1:0.1:10],log10(200./[0.1:0.1:10]),'k-','linew',1);
    h=plot([0.1:0.1:10],log10(100./[0.1:0.1:10]),'k-','linew',1);
    h=plot([0.1:0.1:10],log10(50./[0.1:0.1:10]),'k-','linew',1);

	switch experiment
	case 'XENON10'
		%h=erbc( xebands.NR.be(1:end-1) , xebands.NR.mu , xebands.NR.mu_err , 'bb');set(h,'lineW',1.5);
		h=stairs( xebands.NR.be(1:end-1) , xebands.NR.mu , 'b-');set(h,'lineW',1.5);
	case 'LUX'
		load('./matFiles/LUXbands.mat');
		h=stairs(bands.NR.bc,bands.NR.mc,'k-');set(h,'lineW',1.5);
	end
	
	
	axis([0 xbinz(end) 0.5 2.5]); ax=axis;

	
	xlabel('S1 [phe]');
	ylabel(['log$_{10}$(S2/S1)']);

	
	set(gca,'Layer','top');
	set(gcf,'renderer','painters')
	set(gcf,'color','white')

    colorbar;