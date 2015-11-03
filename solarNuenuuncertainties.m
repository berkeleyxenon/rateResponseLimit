%
%
% 150812 pfs


flux.pp = [5.94 5.94 5.94 5.99 5.99 6.06 6.05];

flux.Be7 = [4.86 4.84 4.88 4.89 4.84 4.34 4.38];

flux.B8 = [5.79 5.74 5.87 5.83 5.69 4.51 4.59];


ff=fieldnames(flux);
for ii=1:length(ff)
stat = std(flux.(ff{ii}))/mean(flux.(ff{ii}));
syst = (max(flux.(ff{ii}))-min(flux.(ff{ii})))/mean(flux.(ff{ii}));

dis('%s 1sigma: %3.1f%%, syst: %3.1f%%',ff{ii},stat*100,syst*100);

    
end
