%function plotLeffpts;
%
% just does shit
pathh='.';

%% LUX 2015 preliminary
x=[1.1501    2.0014    3.2133    4.5595    6.0467    7.6988    9.5670     11.6442   13.9701];
y=[3.9711    5.5415    5.1854    6.2147    5.8170    7.1102    7.4740     7.4799    8.9885];
h=plot(x,y*0.015,'ro');set(h,'markerf','r');hold on;
plot([x;x],0.015*[y-[0.7417    0.4930    0.4976    0.3450    0.3992    0.5429    0.4391  0.4235    0.6975];y+[0.9300    0.4949    0.4319    0.5102    0.3226    0.3554    0.4034   0.5549    0.5042]],'r-');hold on;

%%% Yale Leff 2010
c1 = [0.2 0.5 0.9];
	load([pathh '/dataFiles/yale_leff.mat']);
	%h=erb(yale.Enr,yale.leff,[yale.leff_err_m yale.leff_err_p],'bb',0);
	%dis('need a better erb function');
	for i=1:length(yale.Enr)
		h=plot([yale.Enr_m(i) yale.Enr_p(i)],yale.leff(i)*[1 1],'b-'); set(h,'Color',c1);
		hold on;
		h=plot(yale.Enr(i)*[1 1],yale.leff(i)+[-yale.leff_err_m(i) yale.leff_err_p(i)],'b-');set(h,'Color',c1);
	end
	h=plot(yale.Enr,yale.leff,'bp','markerSize',14);set(h,'markerEdgeColor',c1,'markerFaceColor',[1 1 1]);

%%% Columbia Leff 2009
c2 = [0.8 0.2 0.1];
	load([pathh '/dataFiles/columbia_leff.mat']);
	for i=1:length(columbia.Enr)
		h=plot([columbia.Enr_m(i) columbia.Enr_p(i)],columbia.leff(i)*[1 1],'r-'); set(h,'Color',c2);
		h=plot(columbia.Enr(i)*[1 1],columbia.leff(i)+[-columbia.leff_err_m(i) columbia.leff_err_p(i)],'r-'); set(h,'Color',c2);
	end
	h=plot(columbia.Enr,columbia.leff,'rd','markerSize',10); set(h,'markerEdgeColor',c2,'markerFaceColor',[1 1 1]);

%%% Columbia Leff 2011
	load([pathh '/dataFiles/plante2011.mat']);
	%plante2011.Enr =      [3.00 3.91 5.00 6.48 8.40 10.67 14.80 55.61];
	%plante2011.Enr_m =    [2.39 3.19 4.19 5.51 7.09 9.12  13.46 46.26];
	%plante2011.Enr_p =    [3.63 4.61 5.81 7.48 9.63 12.23 16.07 64.20];
	%plante2011.leff =     [.088 .095 .098 .121 .139 .143  .144  .268];
	%plante2011.leff_pm =  [.015 .016 .014 .010 .011 .011  .008  .014];
	
	for i=1:length(plante2011.Enr)
		h=plot([plante2011.Enr_m(i) plante2011.Enr_p(i)],plante2011.leff(i)*[1 1],'k-');
		h=plot(plante2011.Enr(i)*[1 1],plante2011.leff(i)+[-plante2011.leff_pm(i) plante2011.leff_pm(i)],'k-');
	end
	h=plot(plante2011.Enr,plante2011.leff,'ks','markerSize',10); set(h,'markerFaceColor',[1 1 1]);

%%% XENON10 Leff
if 0
	load([pathh '/dataFiles/xe10.mat']);
	xe10.err=sqrt(xe10.stat_err.^2);%+xe10.s1_err.^2+xe10.s_err.^2);
	h=erb(xe10.Enr,xe10.leff,xe10.err,'k.');
	hold on;
	h=plot(xe10.Enr,xe10.leff,'kv');set(h,'markerSize',10,'markerFaceColor',[1 1 1]);
end
