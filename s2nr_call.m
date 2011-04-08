function  HNR_medio = s2nr_call(filename, NFFT , RTH, alpha)
% Programa de calibracao do Detector de SNR/HNR
% Joao SANSAO, Maurilio VIEIRA, Fev. 2007
%Janela da Transformada de Fourier
NFFT
% Overlap das janelas
NOVERLAP = ceil(alpha* NFFT); % 0.9 ok 
%Flag para reproducao da voz testada => 1 reproduz
reproduz = 0;
%Parametro de Reamostragem:
reamostra = 0;
P =1;
Q =1;


%Parametros gerais:
OFFSET = 10; %offset da imagem/pre-processamento


% incluir verificação se houve passagem de parâmetro
param.eps = 1e-16;
param.blksze1 = 5; 
param.thresh = 1e-1; 
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
param.P = P;
param.Q = Q;
param.show = 0.50; 

%loop principal
DadosSaida = [];


[Ys,FS] = wavread(filename);

%etapa de reamostragem
% if (reamostra)
% 	FS = P * FS / Q; 
% 	Y = RESAMPLE(Ys,P,Q);
% else
% 	Y = Ys;
% end

Y = Ys;
param.TotalTime = (length(Y) - 1)/FS ;
param.FS = FS;

%normalizacao: maior valor absoluto unitario
Y = Y / max(abs(Y));
		
%Geracao do Espectro a ser analisado
		
[AmpliNorm,AmpliNormLog] = GeraEspectro( Y, NFFT, FS, NOVERLAP);

sizeAmpli = size(AmpliNorm);
		
%Pre-processamento : zero-padding 	
EntEspectro = [ zeros(OFFSET,sizeAmpli(2)); AmpliNorm(1:(end),:)];
EntEspectroLog = [ zeros(OFFSET,sizeAmpli(2)); AmpliNormLog(1:(end),:)];
% Chamada Principal, calculo do SNR e HNR
[HNR,SNR] = s2nr_measurement(EntEspectro, EntEspectroLog,param);

valHNR = 10 * log10 (mean(10 .^ (HNR( floor(0.25*end):floor(0.75*end) )/10) ))
valSNR = 10 * log10 (mean(10 .^ (SNR( floor(0.25*end):floor(0.75*end) )/10) ));

HNR_medio = valHNR;

SNR_medio = valSNR;

plot(HNR)


function [AmpliNorm,AmpliNormLog] = GeraEspectro( Y, NFFT, FS, NOVERLAP)
	
B = spectrogram(Y,NFFT,NOVERLAP,NFFT,FS);
Ampli = abs(B);
AmpliNorm = (normalise(Ampli));

AmpliLog = 20*log10(abs(B+1e-8));
AmpliNormLog = (normalise(AmpliLog));


