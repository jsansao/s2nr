function  [S2NR, T_axis, S2NR_mean] = s2nr_call(filename, NFFT , ...
                                                RTH, alpha, sigma)
% S2NR call procedure 
%   usage:
%    [S2NR, T_axis, S2NR_mean] = s2nr_call(filename, NFFT , RTH,
%    alpha)
% NFFT: fft window size 
% RTH: reliability threshold
% alpha: fft window overlap 
% sigma: image segmentation threshold
%
% 
% Example call:
%
%   [S2NR, T_axis, S2NR_mean] = s2nr_call('vowel_a.wav', 512 , ...
%                                                0.6, 0.9, 0.1)
% 
% Joao SANSAO, Maurilio VIEIRA, Feb 2007-Abril 2011
%


% Window overlap
NOVERLAP = ceil(alpha * NFFT); 
%Parametros gerais:
OFFSET = 10; %image offset /pre-processing

% default parameters 
param.eps = 1e-16;
param.blksze1 = 5; 
param.thresh = sigma; 
param.blksze2 = 16;
param.gradientsigma = 1;
param.blocksigma = 5;	
param.orientsmoothsigma = 5;
param.windsze = 5; 
param.minWaveLength = 1;
param.maxWaveLength = 15;
param.kx = 0.15;
param.ky = 0.15;    
param.rthresh = RTH;
param.OFFSET = OFFSET;
param.NFFT = NFFT;
param.NOVERLAP = NOVERLAP;

[Y,FS] = wavread(filename);

param.TotalTime = (length(Y) - 1)/FS ;
param.FS = FS;

Y = Y / max(abs(Y));
		
%Generate spectrographic image
		
[AmpliNorm,AmpliNormLog] = GenerateSpectrum( Y, NFFT, FS, NOVERLAP);

sizeAmpli = size(AmpliNorm);
		
%Pre-processing : zero-padding 	
SpectrumInput = [ zeros(OFFSET,sizeAmpli(2)); AmpliNorm(1:(end),:)];
SpectrumInputLog = [ zeros(OFFSET,sizeAmpli(2)); AmpliNormLog(1:(end),:)];

% Main call, s2nr measurement
[S2NR,SNR] = s2nr_measurement(SpectrumInput, SpectrumInputLog,param);

T_axis = (0:(length(S2NR)-1))/FS;

S2NR_mean = 10 * log10 (mean(10 .^ (S2NR( floor(0.25*end):floor(0.75*end) )/10) ))



function [AmpliNorm,AmpliNormLog] = GenerateSpectrum( Y, NFFT, FS, NOVERLAP)
	
B = spectrogram(Y,NFFT,NOVERLAP,NFFT,FS);
Ampli = abs(B);
AmpliNorm = (normalise(Ampli));

AmpliLog = 20*log10(abs(B+1e-8));
AmpliNormLog = (normalise(AmpliLog));


