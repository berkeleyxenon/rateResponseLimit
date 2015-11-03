%function plotQypts;
%
% just does shit

pathh='.';

%% LUX 2015 preliminary
x=[ 0.7654    1.1556    1.4881    2.0331    2.7844    3.8275    5.4339     7.8016   11.1963   16.1432];
y=[ 8.8842    7.4423   10.4829    9.0010    8.2603    8.0731    7.7844     7.5754    6.5290    5.8848];
h=plot(x,y,'ro');set(h,'markerf','r');hold on;
plot([x ; x],[y-[2.0179    1.7826    1.8164    0.7728    0.5029    0.3807    0.2231    0.1552    0.1163    0.1294 ]; y+[ 2.3897    2.3226    1.7866    0.7933    0.6103    0.3576    0.3035     0.2171    0.1559    0.1436]],'r-');
%%% plot Yale Qy
c1 = [0.2 0.5 0.9];
load([pathh '/dataFiles/yale_Qy.mat']);
	for ii=1:length(yale.Enr)
		h=plot(yale.Enr(ii)+[-yale.Enr_err(ii) +yale.Enr_err(ii)] , yale.Qy(ii)*[1 1] , 'b-');set(h,'Color',c1);
		hold on;
		h=plot(yale.Enr(ii)*[1 1] , yale.Qy(ii)+[-yale.Qy_err(ii) +yale.Qy_err(ii)] , 'b-');set(h,'Color',c1);
	end
	h=plot(yale.Enr , yale.Qy , 'bp','markerSize',14); set(h,'markerEdgeColor',c1,'markerFaceColor',[1 1 1]);

%%% plot XENON10 Qy
switch -1 % use 0 !!
case 0
	load([pathh '/dataFiles/2009_XENON10_NR_Qy_multipleScatterMethod.mat']);
	h=erb(Q.Er,Q.Qy,[Q.Qy_err_m Q.Qy_err_p],'k.');
	hold on;
	h=plot(Q.Er,Q.Qy,'kv','markers',10,'markerFaceColor',[1 1 1]);
case 1			
	load([pathh './dataFiles/2009_XENON10_NR_Qy.mat']);
	[tmp_leff_1,tmp_Er_1]=getLeff(Q.keVr,'XENON10');
	%[tmp_leff_2,tmp_Er_2]=getLeff(Q.keVr,'XENON10_2010_db8');
	[tmp_leff_2,tmp_Er_2]=getLeff(Q.keVr,'XENON10'); % no change from published version

	h=erb(Q.keVr.*tmp_leff_1./tmp_leff_2,Q.Qy,Q.Qy_err,'k.');
	hold on;
	h=plot(Q.keVr.*tmp_leff_1./tmp_leff_2,Q.Qy,'kv');set(h,'markers',10,'markerFaceColor',[1 1 1]*0);
end

if 1 % PRL 97 081302 (2006) data
	nevis.x = [25 35 45 55 75 105 115];
	nevis.yl = [4.70 4.00 3.60 3.37 2.9 2.5 2.35];
	nevis.yh = [5.45 4.40 3.90 3.62 3.2 2.7 2.57];
	cwru.x = [35 45 55 75 105 115];
	cwru.yl = [3.40 3.20 2.91 2.55 2.25 2.04];
	cwru.yh = [3.70 3.40 3.20 2.80 2.4 2.17];
	ms = 7;
	%h=plot( nevis.x , nevis.yl , 'ko','markerFaceColor','k');set(h,'markers',ms);
	h=plot( nevis.x , nevis.yh , 'ks','markerFaceColor','k');set(h,'markers',ms);
	h=plot( cwru.x , cwru.yl , 'ko');	set(h,'markers',ms);
	%h=plot( cwru.x , cwru.yh , 'ks');	set(h,'markers',ms);
	%leg{end+1} = 'Phys. Rev. Lett. 081302 (2006)'; hleg(end+1) = h;
end

