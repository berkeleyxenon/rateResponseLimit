function [Nout] = PoissonConvolution(raw,bins,c,bin_flag);
% [Nout] = PoissonConvolution(raw,bins,c,bin_flag);
%
% inputs:
% raw: input data
% bins: the energy bins you want for the histogram
% c: Poisson stats -> c=1
%		worse res. than Poisson -> c > 1
%		better res. than Poisson -> c < 1
% bin_flag=0 (default) 'raw' is already binned (raw = bincounts)
% bin_flag=1 raw data will be binned into 'bins' by histConvolution
%
% outputs:
% Nout: number of counts in each bin, post-convolution
%
% v1.0 pfs - should make 'c' an input arg
% v2.0 pfs 071030 - totally re-wrote it -- note, c is relative to Poisson R
% v3.1 rjg 071120 - include temporary store for resolution function to improve code speed
% v3.2 pfs 080129 - bin counts are too high below bin center (relative to Poisson) -- trying to fix this


% These are the globals that store lookup table
global hC_cscaling hC_bins hC_table

%raw=raw(:);
%bins=bins(:);

if nargin < 4 || isempty(bin_flag)
	bin_flag = 0; % default
end
if nargin<3
	c = 1; % Poisson R
end

if bin_flag == 0 % expect already histogram'd data
	Nin=raw;
else % expecting rawdata
	Nin = hist(raw,bins);
end

cscaling = 1.0 / c;
if c~=1
%	dis('Using non-std cscaling = %4.1f',cscaling);
end

Nout = zeros(size(Nin));
% determine ii_start
tmp=find(bins>=0);
ii_start = tmp(1);
dbins = bins(2)-bins(1);


iimode = 3; % How to perform resolution adjustment
			% 2 - conventional
			% 3 - lookup table (faster)

switch iimode  % (See above)
case 3 % Use look up table

	% Need to check if lookup table needs revising, due to change in calling parameters
	if (isempty(hC_cscaling) || cscaling ~= hC_cscaling) ...
		|| (isempty(hC_bins) || (bins(2)-bins(1))~= (hC_bins(2)-hC_bins(1)) || any(bins ~= hC_bins))
		%  Something has changed so calculate lookup table  hc_table
		% 	Note that it must be a square matrix
			dis('PoissonConvolution.m - recalculating convolution table');
			hC_table = zeros(length(Nin),length(Nin));
			hC_cscaling = cscaling;
			hC_bins = bins;
			bins_range = [bins-dbins bins(end)]; % gives best agreement w Poisson dist.
			% Entries < ii_start get left zero which is fine				
			for ii = ii_start:length(bins)	
				tmp2 = 0:ceil(cscaling.*(bins(end)+dbins/2));	% phe values to model
				tmp = interp1( ...  % Look up CDF at bin edges coverted to phe
					tmp2 ...
					,poisscdf(tmp2,cscaling.*bins(ii)) ...
					,cscaling.*bins_range ... % orig
					...,'cubic' ... leads to lumpy interpolation
					,'spline' ...
					,0 ... 	% Returns zero when out of range
					);

				interpd_poiss_pdf = diff(tmp);
				
				hC_table(:,ii) = interpd_poiss_pdf;

				if 0 % debug
					interpd_poiss_pdf(ii_start:10) 
					figure(1e6);clf;
					subplot(2,1,1);
						plt( bins,tmp(1:end-1),'bo-');
						axis([-1 6 0 1.1]);
						title(dis('ii = %d',ii));
					subplot(2,1,2);
						plt( bins,interpd_poiss_pdf,'ro-');
						set(gca,'ysc','log');
						axis([-1 6 1e-4 1.1]);
					title('YOU HAVE A DEBUG FLAG SET');
					drawnow;
					%keyboard;
				end
			end
	end		

	% ACTUAL CONVOLUTION DONE HERE:
	Nout = hC_table * Nin';
	Nout = Nout';

case {1,2}

	for ii = ii_start:length(bins)
		switch iimode
			case 1 % Original discrete look up - had normalization issues
				y_poiss = poisspdf( ...
				round(cscaling.*bins) , cscaling.*bins(ii) ...
				); % Poisson for 3 phe/keVee - need also to round to nearest integer which is a bit artificial
			%Òkeyboard

			case 2 % Do resolution calc explicitly
				tmp2 = 0:ceil(cscaling.*(bins(end)+dbins/2));	% phe values to model
				tmp = interp1( ...  % Look up CDF at bin edges coverted to phe
					tmp2 ...
					,poisscdf(tmp2,cscaling.*bins(ii)) ...
					,cscaling.*[bins-dbins/2 bins(end)+dbins/2] ...
					,'linear' ...
					,0 ... 	% Returns zero when out of range
					);
				y_poiss = diff(tmp);

			end
			Nout = Nout + Nin(ii) .* y_poiss;
	end
end	

