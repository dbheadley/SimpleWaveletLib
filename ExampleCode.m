%% enter your matlab work directory here
matlabWork = 'C:\Users\dbh60\Documents\MATLAB\';

%% load test data
load([matlabWork 'SimpleWaveletLib\TestData.mat'])
dt = 0.001; % sample step size in seconds
freqs = 2.^(0:0.125:8); % frequencies to calculate in Hz


%% Caclulate wavelet spectrogram
blaFileName = WaveletSpec(blaSamp,dt,freqs, 'SAVE', ...
  [matlabWork 'SimpleWaveletLib\TestWaveletBLASpec.wls'], 'NORMALIZE');
vhFileName = WaveletSpec(vhSamp,dt,freqs, 'SAVE', ...
  [matlabWork 'SimpleWaveletLib\TestWaveletVHSpec.wls'], 'NORMALIZE');

%% load your spectrograms as wavelet objects
blaObj = WaveletObj(blaFileName);
vhObj = WaveletObj(vhFileName);

%% calculate coherence between spectrograms
cohFileName = WaveletCoherence(blaObj,vhObj,'SAVE', ...
  [matlabWork 'SimpleWaveletLib\TestWaveletBLA-VHCoh.wls']);

%% load your coherence object
cohObj = WaveletObj(cohFileName);

%% Get spectral power and coherence
blaSpec = blaObj.LoadAmplitude('SEGS', [10 20]);
vhSpec = vhObj.LoadAmplitude('SEGS', [10 20]);
cohSpec = cohObj.LoadCoherence('SEGS', [10 20]);

%% plot all the data
figure;
subplot(3,1,1);
h = pcolor(blaSpec.Times,blaSpec.Freqs,blaSpec.Amplitude);
set(h,'LineStyle','none');
set(gca,'yscale','log');
set(gca,'TickDir','out');
grid on;

subplot(3,1,2);
h = pcolor(vhSpec.Times,vhSpec.Freqs,vhSpec.Amplitude);
set(h,'LineStyle','none');
set(gca,'yscale','log');
set(gca,'TickDir','out');
grid on;

subplot(3,1,3);
h = pcolor(cohSpec.Times,cohSpec.Freqs,cohSpec.Coherence);
set(h,'LineStyle','none');
set(gca,'yscale','log');
set(gca,'TickDir','out');
grid on;