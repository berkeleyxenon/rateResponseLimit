% plot output of getSigmaWN.m
%
% 130820 pfs - constantly evolving, pls make your own copy if you want to edit


define_rainbow;
hleg=[];leg=[];

figure(1); clf;
	limitsPlot=axes('Position',[0.15    0.15    0.80    0.75]);% [left, bottom, width, height] %subplot(2,1,1);
	hold on;
	if 1 % 2010 Hooper et al.
		load('2010.Hooper.CoGeNT.mat');
		%hcol = [1 1 1]*0.7;
		%h=fill([hooper.Mchi hooper.Mchi(1)],[hooper.Schi hooper.Schi(1)],hcol,'FaceColor',hcol,'EdgeColor','k','FaceAlpha',0.5);
		h=plot([hooper.mchi hooper.mchi(1)],1e-36*[hooper.schi hooper.schi(1)],'g-','lineW',1);
	end
	if 0 % 2010 Hooper et al.
		load('2010.Hooper.CoGeNT.DAMA.allowed.mat');
		%hcol = [1 1 1]*0.7;
		%h=fill([hooper.Mchi hooper.Mchi(1)],[hooper.Schi hooper.Schi(1)],hcol,'FaceColor',hcol,'EdgeColor','k','FaceAlpha',0.5);
		h=plot([hooper.Mchi hooper.Mchi(1)],[hooper.Schi hooper.Schi(1)],'k-','lineW',1);
	end	
	if 0
		load XENON10_PRL.mat
		hold on;
		%xx=[10:1:1e3];yy=interp1(XENON10_PRL.m_chi,XENON10_PRL.sig_n,xx,'cubic');
		%plot(xx,yy,'r-');
		plot(XENON10_PRL.m_chi,XENON10_PRL.sig_n,'k-');
	end
	

	if 0 % CDMS II low-energy analysis
		%ahmed2010.Mchi = [4.532 4.712 5.038 5.312 5.586 5.912 6.358 6.967 7.464 8.021 8.759 9.513 10.311 11.083 11.983];
		%ahmed2010.Schi = [[9.906 6.591 3.628 2.437 1.733 1.232]*1e-40 [8.593 5.938 4.642 3.662 2.836 2.302 1.923 1.716 1.461]*1e-41];
		%h=plot([ahmed2010.Mchi],[ahmed2010.Schi],'k--');
		load 2010.CDMS.moore.limits.mat;
		h=plot([moore2010.Mchi],[moore2010.Schi],'m-');
		hleg(end+1)=h;
		leg{end+1}=['CDMS II, PRL 106 131302 (2011)'];
		h=text(13,1.5e-41,'CDMS','rotation',-12,'FontS',12,'Color','m','interpreter','latex');		
	end
	if 0 % CDMS low-threshold analysis
		load 2010.CDMS.bunker.limits.mat;
		h=plot([akerib2010.Mchi],[akerib2010.Schi],'k:','lineW',1);
		hleg(end+1)=h;
		leg{end+1}=['Ref. [11]'];
	end
	if 0 % XENON100 100 days
		load 2011.XENON100.100days.mat;
		h=plot(m_x,sigma_x,'r-');
		hleg(end+1)=h;
		leg{end+1}=['XENON100, PRL 107 131302 (2011)'];
	end
	if 1 % XENON10 S2-only
		%load 2011.XENON10.S2.mat;
		load 2013.XENON10.S2.mat
		h=loglog(xenon10.m_x,xenon10.sigma_x,'k-','lineW',1);%set(h,'color',cols.ora);%,'markerfacecolor',[1 1 1]);	
		hleg(end+1)=h;
		leg{end+1}=['XENON10, PRL 107 051301 (2011)'];
	end
	if 0 % CRESST 2002
		load('2002.Angloher.CRESST.mat');
		h=plot(angloher.mchi,angloher.schi,'g-');
		hleg(end+1)=h;
		leg{end+1}=['CRESST, Astropart. Phys. 18 43 (2002)'];
		%angloher.mchi = [1 2 3 4 5 10 15 20];
		%angloher.schi = [7.86e-38 1.05e-38 4.74e-39 2.99e-39 2.22e-39 1.26e-39 1.00e-39 9.32e-40];
		%save('2002.Angloher.CRESST.mat','angloher');
	end
	if 1
		xenon100.m_x =     [6     7      8     9   10  12 16 20  25  30  40  50  60  80  100 200];
		xenon100.sigma_x = [7.5e5 1.5e4  2.2e3 675 288 85 19 7.5 4.2 3.0 2.2 2.0 2.1 2.2 2.8 4.5]*1e-45;
		%ee=[6:0.1:100];
		%h=plot(ee,interp1(xenon100.m_x,xenon100.sigma_x,ee),'b-','linew',1.5);
		h=plot(xenon100.m_x,xenon100.sigma_x,'b-','linew',1);
		%set(h,'Color',[0.2 0.5 0.9]);
	end
	if 0
		xallow.m = [6.0 6.25 6.5 7.0 7.5 8.0 8.5 9.0 8.5 8.0 7.5 7.0 6.5 6.25 6.0];
		xallow.s = [10  13   12  7.5 4.3 2.2 1.0 0.5 0.4 0.6 1.0 2.0 4.5 7   10]*1e-43;
		plot(xallow.m,xallow.s,'k-');
	
	end
	if 0
		h=text(10,2.1e-40,'DAMA','rotation',-8,'FontS',fs,'interpreter','latex');		
		h=text(15.5,0.45e-40,'CoGeNT','rotation',0,'FontS',fs,'interpreter','latex');
			h=scribe.arrow;
			set(h,'X',[0.79 0.70]);
			set(h,'Y',[1 1]*0.585);
    end
    if 1
        lux.m_x =     [7 7.5 8 9 10 12 14 18 25 30 40 50 60 80 100 1e3];
        lux.sigma_x = [0.76 0.02 0.001 2.079e-4 7.405e-5 2.189e-5 9.410e-6 4.059e-6 2.368e-6 1.885e-6 1.574e-6 1.711e-6 1.764e-6 2.014e-6 2.363e-6 1.943e-5]*1e-40;
        xx = [7:0.5:150];
        yy = interp1(lux.m_x,lux.sigma_x,xx,'cubic',0);
        %h=plot(lux.m_x,lux.sigma_x,'ko');
        h=plot(xx,yy,'k--');
    end
    if 1
    	cdmse.m_x = [10 12 14 16 18 20 23 24 25 27 30 34 35 37 42 50 52 60 71 73 100 150];
    	cdmse.sigma_x = [841.5 231.9 88.28 51.43 34.75 26.90 21.72 20.46 18.64 15.55 12.75 10.65 10.09 8.688 6.446 4.764 4.484 3.924 3.643 3.363 3.223 3.784]*1e-44;
        xx = [7:0.5:150];
        yy = interp1(cdmse.m_x,cdmse.sigma_x,xx,'cubic',0);
        h=plot(xx,yy,'r-');
    
    end
	
	if 1 % CDMS Si 2013
		cdms2013.m_x=[5.726 5.650 5.842 6.188 6.662 7.495 8.573 10.455 11.770 13.517 15.850 19.045 20.332 20.757 20.409 19.038 17.280 10.392 8.933 7.333 6.530 5.726];
		cdms2013.sigma_x=[4.476e-40 3.658e-40 1.297e-40 5.394e-41 1.934e-41 6.997e-42 3.161e-42 1.534e-42 1.167e-42 1.001e-42 9.445e-43 1.082e-42 1.298e-42 1.522e-42 1.888e-42 2.772e-42 3.961e-42 2.173e-41 4.300e-41 1.274e-40 2.543e-40 4.476e-40];
		plot(cdms2013.m_x,cdms2013.sigma_x,'c-','lineW',2);
	end
	

	if exist('sn','var')
		% run several trial experiments
		%clear all;for iii=1:10;getSigmaWN;sn(:,iii)=S{ee}.sigma_n;end;plotLimits;%_2013
		clear mu sigma;
		for mm=1:size(sn,1)
			this_sn = sort(sn(mm,:)); %this_sn = this_sn(2:end-1);
			[N,X]=hist(this_sn,10);
			[mu(mm),sigma(mm)] = wmean(N,X);
		end
		% was previously using mean(sn,2), but wmean is better
		%h=plot(S{ee}.m_chi,mu,'m-');
		plot(S{ee}.m_chi,mean(sn,2),'m-');
	else
		% this trial
		h=plot(S{ee}.m_chi,S{ee}.sigma_n,'k-'); 
	end

	% projected S2 only, average of several trial experiments	
	if 0 % LZ, 131016 checking pp effect w 36 trials
		mm=[4.0 4.5000    5.0000    6.0000    7.0000    8.0000   10.0000   12.0000];
		% with pp, no resolution
		plot([4.2 4.5 5 6 7 8 10 12],1e-43*[6.3213 1.8061 0.4650 0.0940 0.0285 0.0132 0.0038 0.0019],'gp');
		% with pp
		plot(mm,1e-43*[1.19 0.5272 0.2427  0.0729  0.0268 0.0116 0.0038 0.0019],'m-');
		% no pp
		plot(mm,1e-43*[1.00 0.4655 0.2148  0.0621 0.0233 0.0108  0.0033 0.0016],'m--');
		% Matt
		plot([4.2 4.5 5 6 7 8 9 10],[1.0e-42 2.0e-43 3.9e-44 6.5e-45 2.41e-45 1.0e-45 5.35e-46 3.2e-46],'mo');
		h=text(4.5,3.0e-45,'LZ (projected S2-only)','rotation',-25,'FontS',12,'interpreter','latex');
	end
	
	if 0
	% 20130927 LUX projected S2 only, average of 10 trial experiments
		mx = [4.00 4.50 5.00 5.50 6.00 6.50 7.00 8.00 9.00 10.00 11.00 12.00 14.00 16.00 20.00];
		sx = 1e-41*[0.5066    0.1990    0.0890    0.0459    0.0301    0.0216    0.0162    0.0106    0.0077    0.0060    0.0050     0.0043    0.0034    0.0028    0.0023];
		h=plot(mx,sx,'b--');
	% 20130920 LZ projected S2 only, average of 10 trial experiments
		mx = [4.00 4.50 5.00 5.50 6.00 6.50 7.00 8.00 9.00 10.00 11.00 12.00];
		sx = [2.20e-44 1.00e-44 5.63e-45 3.64e-45 2.59e-45 1.97e-45 1.59e-45 1.00e-45 5.21e-46 3.06e-46 2.02e-46 1.44e-46];
		h=plot(mx,sx,'g--');
	end
		
	if 0
		xdata = [xenon10.m_x(13:21) [7.5 8 9 10 12 14] reverse(xenon10.m_x(13:27)) xenon10.m_x(13)];
		ydata = [xenon10.sigma_x(13:21) [0.02 0.001 2.079e-4 7.405e-5 2.189e-5 9.410e-6 ]*1e-40 reverse(xenon10.sigma_x(13:27))*1e-3 xenon10.sigma_x(13)];
		h=fill(xdata,ydata,'k');set(h,'facecolor',[1 1 1]*0.8,'edgecolor',[1 1 1]);
	end

	if 1
		h=text(35,2.1e-43,'CDMS +','rotation',-17,'FontS',12,'Color','r','interpreter','latex');		
		h=text(58,7.0e-44,'Edelweiss','rotation',-4,'FontS',12,'Color','r','interpreter','latex');		

		h=text(50,3.5e-45,'XENON100','rotation',+5,'FontS',12,'Color','b','interpreter','latex');		

		h=text(7,1.9e-40,sprintf('CoGeNT\n+ DAMA'),'rotation',-16,'FontS',12,'Color','g','interpreter','latex');		
		h=text(21,1.0e-42,'CDMS Si','rotation',-10,'FontS',12,'Color','c','interpreter','latex');		

		%h=text(4.7,5.0e-40,'XENON10','rotation',-54,'FontS',12,'interpreter','latex');		
		h=text(21,3.0e-42,'XENON10*','rotation',-4,'FontS',12,'interpreter','latex');		

		h=text(30,3.5e-46,'LUX','rotation',-5,'FontS',12,'interpreter','latex');		
		h=text(40,3.0e-46,'(projected)','rotation',+3,'FontS',12,'interpreter','latex');		

		%annotation(gcf,'arrow',[0.21 0.21],[0.27 0.42]);

	
	end
	% set plot
	yax=[1e-47 1e-39];
	xax = [2 110];
	set(gca,'xsc','log'); %
	set(gca,'ysc','log'); %
	axis([xax yax]); ax=axis;
	h=text(14,0.2*ax(3),'$m_{\chi}$ ~~~[GeV]');
	set(gca,'xtick',[2:1:10 20:10:100],'xtickl',[]);	
		for xtik = [5 10 20 50 100];
			h=text(xtik,ax(3)*0.5,dis('%1.0f',xtik),'HorizontalAlignment','center');
		end
	
	ylabel('$\sigma_n$ ~~~[cm$^2$]');
	set(gca,'ytick',10.^[log10(yax(1)):log10(yax(2))]);
	set(gca,'Layer','top');
	set(gcf,'Renderer','painters');

	box on;
	grid on;



if 0
save2pdf('sigma_v_m');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if 0 % for snowmass layover
	figure(2);clf;
	h=loglog(S{ee}.m_chi,mean(sn,2),'--'); set(h,'Color',[0 1 0.2]);
	axis([0.5 1e4 1e-49 1e-39]);
	makeTransparent;
	
	for ii=1:size(sn,1)
	dis('%2.1f %3.2e',S{ee}.m_chi(ii) , mean(sn(ii,:),2));
	end
end

