typ=0.31;
NRbandsim; plotResult;
bibe_S1 = S1(kut);
bibe_yy = yy(kut);

typ=0.3;
NRbandsim;
plotResult;
ybe_S1 = S1(kut);
ybe_yy = yy(kut);

typ=0;
source = 1;
NRbandsim;
plotResult;
b8_S1 = S1(kut);
b8_yy = yy(kut);

close all;

figure();
hold on;
plot(xx,interp1(NRsimBands.bc,(NRsimBands.m+3*NRsimBands.s),xx,'pchip'),'r-','linew',0.5);
plot(xx,interp1(NRsimBands.bc,NRsimBands.m,xx,'pchip'),'r-','linew',2);
plot(xx,interp1(NRsimBands.bc,(NRsimBands.m-3*NRsimBands.s),xx,'pchip'),'r-','linew',0.5);

%%
num_bins = 10;
%bibe_S1 = bibe_S1(bibe_S1 < 7);
%bibe_yy = bibe_yy(bibe_S1 < 7);
[minBiBe_S1 minBiBeS1Ind] = min(bibe_S1);
[maxBiBe_S1 maxBiBeS1Ind] = max(bibe_S1);
dBi = (maxBiBe_S1 - minBiBe_S1)/(num_bins + 1);
biBins = minBiBe_S1:dBi:maxBiBe_S1;
bibe_max = zeros(1, num_bins);
bibe_min = zeros(1, num_bins);
bibe_x = zeros(1, num_bins);
for i = 1:num_bins
    bibe_min(i) = min(bibe_yy(bibe_S1>biBins(i) & bibe_S1 < biBins(i + 1)));
    bibe_max(i) = max(bibe_yy(bibe_S1>biBins(i) & bibe_S1 < biBins(i + 1)));
    bibe_x(i) = mean([biBins(i) biBins(i + 1)]);
end
bibe_x = [minBiBe_S1 bibe_x maxBiBe_S1 fliplr(bibe_x) minBiBe_S1];
bibe_y = [bibe_yy(minBiBeS1Ind) bibe_min bibe_yy(maxBiBeS1Ind) fliplr(bibe_max) bibe_yy(minBiBeS1Ind)];

%ybe_S1 = ybe_S1(ybe_S1 < 8);
%ybe_yy = ybe_yy(ybe_S1 < 8);
[minYBe_S1 minYBeS1Ind] = min(ybe_S1);
[maxYBe_S1 maxYBeS1Ind] = max(ybe_S1);
dY = (maxYBe_S1 - minYBe_S1)/(num_bins + 1);
yBins = minYBe_S1:dY:maxYBe_S1;
ybe_max = zeros(1, num_bins);
ybe_min = zeros(1, num_bins);
ybe_x = zeros(1, num_bins);
for i = 1:num_bins
    ybe_min(i) = min(ybe_yy(ybe_S1> yBins(i) & ybe_S1 < yBins(i + 1)));
    ybe_max(i) = max(ybe_yy(ybe_S1> yBins(i) & ybe_S1 < yBins(i + 1)));
    ybe_x(i) = mean([yBins(i) yBins(i + 1)]);
end
ybe_x = [minYBe_S1 ybe_x maxYBe_S1 fliplr(ybe_x) minYBe_S1];
ybe_y = [ybe_yy(minYBeS1Ind) ybe_min ybe_yy(maxYBeS1Ind) fliplr(ybe_max) ybe_yy(minYBeS1Ind)];

%%
[minB8_S1 minB8S1Ind] = min(b8_S1);
[maxB8_S1 maxB8S1Ind] = max(b8_S1);
dB8 = (maxB8_S1 - minB8_S1)/(num_bins + 1);
b8Bins = minB8_S1:dB8:maxB8_S1;
b8_max = zeros(1, num_bins);
b8_min = zeros(1, num_bins);
b8_x = zeros(1, num_bins);
for i = 1:num_bins
    b8_min(i) = min(b8_yy(b8_S1> b8Bins(i) & b8_S1 < b8Bins(i + 1)));
    b8_max(i) = max(b8_yy(b8_S1> b8Bins(i) & b8_S1 < b8Bins(i + 1)));
    b8_x(i) = mean([b8Bins(i) b8Bins(i + 1)]);
end
b8_x = [minB8_S1 b8_x maxB8_S1 fliplr(b8_x) minB8_S1];
b8_y = [b8_yy(minB8S1Ind) b8_min b8_yy(maxB8S1Ind) fliplr(b8_max) b8_yy(minB8S1Ind)];

plot(bibe_x, bibe_y, 'b-', ybe_x, ybe_y, 'm-', b8_x, b8_y, 'k-', 'linew', 2);
xlabel('S1 (phd)', 'fontsize', 20);
ylabel('log10(S2/S1)', 'fontsize', 20);
legend('NR + 3\sigma', 'NR mean', 'NR - 3\sigma', '^{205}BiBe', '^{88}YBe', '^{8}B');
xlim([0 10]);
set(gca,'FontSize',20);
save_graphic('~/PhotoneutronSources.pdf', [10 6], 'pdf');
