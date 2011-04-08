% Calculo_snr 
%

%
%
% Adaptado de Peter Kovesi ("Function to demonstrate use of fingerprint code")
% http://www.csse.uwa.edu.au/~pk
%
% Usage:  [HNR,SNR] =  testfin(im);
%
% Argument:   im -  Imagem do Espectro pre-processado a ser analisado
%
% Returns:    HNR - Relacao Harmonico Ruído Temporal
%             SNR - Relacao Sinal Ruido Temporal
%             mask  - Ridge-like regions of the image
%             reliability - 'Reliability' of orientation data
%
%
%   Joao SANSAO / Maurilio Nunes Vieira
%   Fevereiro 2007

    
function [HNR,SNR] =  s2nr_measurement(im, imlog,param)
    

FonteLegenda = 14; 


[FreqS TempS] = size(im);
%A = 40 * log10( (1:FreqS));
%B = repmat(  (A)' , 1 ,TempS);
%imagesc(B);

%pause;


      
    if nargin == 1
	eps = 1e-12;
	blksze1 = 4; 
	thresh = 1e-2;
	blksze2 = 16;
	FS = 22050;
	gradientsigma = 1;
	blocksigma = 5;	
	orientsmoothsigma = 5;
	windsze = 5; 
	minWaveLength = 1;
	maxWaveLength = 15;
 	kx = 0.05;
    	ky = 0.05;    
	rthresh = 0.6;
	OFFSET = 10;
	TotalTime = 3;
  else
	eps = param.eps;
	blksze1 = param.blksze1; 
	thresh = param.thresh;
	blksze2 = param.blksze2;
	FS = param.FS;
	gradientsigma = param.gradientsigma;
	blocksigma = param.blocksigma;	
	orientsmoothsigma =param.orientsmoothsigma;
	windsze = param.windsze; 
	minWaveLength = param.minWaveLength;
	maxWaveLength = param.maxWaveLength;
 	kx = param.kx;
    	ky = param.ky;    
	rthresh = param.rthresh;
	OFFSET = param.OFFSET;
	TotalTime = param.TotalTime;
    end
    
    % Identify ridge-like regions and normalise image
    % [normim, mask] = ridgesegment(B .* im, blksze1, thresh);
   [normim, mask] = ridgesegment(imlog, blksze1, thresh);    
    % Determine ridge orientations 
    %[orientim, reliability] = ridgeorientation(im, gradientsigma,...
    %                                              blocksigma, ...
    %                                              orientsmoothsigma)
   [orientim, reliability] = ridgeorient(normim, gradientsigma, blocksigma, orientsmoothsigma);

    
    % Determine ridge frequency values across the image

    %[freqim, medianfreq] =  ridgefreq(im, mask, orientim, blksze, windsze, ...
    %                                 minWaveLength, maxWaveLength)


    [freq, medfreq] = ridgefreq(normim, mask, orientim, blksze2, windsze, minWaveLength, maxWaveLength);


imagesc(freq);    
%pause;
    
    % Actually I find the median frequency value used across the whole
    % fingerprint gives a more satisfactory result...
    freq = medfreq .* mask;



    % Now apply filters to enhance the ridge pattern
    
    %newim =  ridgefilter(im, orientim, freqim, kx, ky, showfilter) 
   
    newim = ridgefilter(normim, orientim, freq, kx, ky, 0);
    
    % Binarise, ridge/valley threshold is 0
    
binim = newim > 0;
%[newim2, binim] = ridgesegment(newim, blksze1, thresh)
%figure
%imagesc(newim2);colorbar;
%pause;

    % Display binary image for where the mask values are one and where
    % the orientation reliability is greater than rthresh
    
    MascaraR = reliability>rthresh;
    Mascara = binim.*mask .* MascaraR; 

    MascaraNew = Mascara;
    
    ImM = double(MascaraNew) .* double (im);
    
    ImMNot = double(not(MascaraNew)) .* double(im);
    ImMNot = medfilt2(ImMNot);

    SImH =   sum( ImM .^2 );
    SImM =   sum( im .^2 );
    SImMNot = sum( ImMNot .^2 );

    Temp = (SImM+eps) ./ (SImMNot+eps) ;
    Temp2 = (SImH+eps) ./ (SImMNot+eps) ;
    Temp_filt = medfilt1(Temp, 16);
    Temp2_filt = medfilt1(Temp2, 16);
    HNR = 10 * log10 (Temp2_filt);
    SNR = 10 * log10 (Temp_filt);

    [Y,X] = size(normim);
    FundoEscala = FS/2;
    EixoY = ((0 : (Y - OFFSET - 1)) / (Y - OFFSET)) * FundoEscala;
    EixoX = ( 0 : X - 1) / (X ) * TotalTime;

    OFFSET_1 = OFFSET + 1;

    %N_escala = 1024;    
    %escala_cor = 0:(1/N_escala):1;
    %a = 1-[ escala_cor; escala_cor; escala_cor]'; 


    %    figure(1)
    %    subplot(2,2,1),imagesc( EixoX, EixoY, im(OFFSET_1:end,:) ); 
    %    colormap(a); ylabel('Frequency(Hz)','FontSize', FonteLegenda); axis xy;
    %    axis([ 0 TotalTime 0 FundoEscala/2]); 
    %title('Normalized Image','FontSize', FonteLegenda );
    %    title('Normalized Image (A_{n})','FontSize', FonteLegenda );


    %    subplot(2,2,2),imagesc(EixoX, EixoY, imlog(OFFSET_1:end,:));
    %    colormap(a); ylabel('Frequency(Hz)','FontSize', FonteLegenda); axis xy;
    %    axis([ 0 TotalTime 0 FundoEscala/2]);
    %title('Log Normalized Image','FontSize', FonteLegenda );
    %    title('Log Normalized Image (A_{log})','FontSize', FonteLegenda );


    %    subplot(2,2,3),imagesc(EixoX, EixoY, normim(OFFSET_1:end,:));
    %    colormap(a); ylabel('Frequency(Hz)','FontSize', FonteLegenda); axis xy;
    %	axis([ 0 TotalTime 0 FundoEscala/2]);
    %title('Normalized Image','FontSize', FonteLegenda );
    %    title('Normalized Image after segmentation (A_{s})','FontSize', FonteLegenda );

    %    subplot(2,2,4),imagesc(EixoX, EixoY, mask(OFFSET_1:end,:));
    %    colormap(a); ylabel('Frequecy(Hz)','FontSize', FonteLegenda); axis xy;
    %	axis([ 0 TotalTime 0 FundoEscala/2]);
    %title('Normalized Image','FontSize', FonteLegenda );
    %    title('Mask (M_{s})','FontSize', FonteLegenda );
   
    %subplot(2,2,4),imagesc(EixoX, EixoY,reliability(OFFSET_1:end,:));
    %colormap(1-a); ylabel('Frequencia(Hz)','FontSize', FonteLegenda); axis xy;
    %axis([ 0 TotalTime 0 FundoEscala/2]);title('Confiabilidade (R)','FontSize', FonteLegenda );

    %   saveas(gca, 'fig1.eps', 'eps');


    %    figure(2)
    %    subplot(2,2,1),imagesc(EixoX, EixoY, mask(OFFSET_1:end,:));
    %    colormap(a); ylabel('Frequency(Hz)','FontSize', FonteLegenda);axis xy;
    %    axis([ 0 TotalTime 0 FundoEscala/2]);title('Segmentation mask (M_{s})','FontSize', FonteLegenda );

    %    subplot(2,2,2),imagesc(EixoX, EixoY, binim(OFFSET_1:end,:));
    %    colormap(a); ylabel('Frequency(Hz)','FontSize', FonteLegenda);axis xy;
    %    axis([ 0 TotalTime 0 FundoEscala/2]);title('Filter mask (M_{G})','FontSize', FonteLegenda );

    %    subplot(2,2,3),imagesc(EixoX, EixoY, MascaraR(OFFSET_1:end,:));
    %    colormap(a); ylabel('Frequency(Hz)','FontSize', FonteLegenda);axis xy;
    %    axis([ 0 TotalTime 0 FundoEscala/2]); title('Reliability mask (M_{R})','FontSize', FonteLegenda );

    %    subplot(2,2,4),imagesc(EixoX, EixoY, Mascara(OFFSET_1:end,:));
    %    colormap(a); ylabel('Frequency(Hz)','FontSize', FonteLegenda);axis xy;
    %    axis([ 0 TotalTime 0 FundoEscala/2]);title('Signal mask(M_{signal})','FontSize', FonteLegenda );
    %   saveas(gca, 'fig2.eps', 'eps');


    %    figure(3)
    %    subplot(2,2,1),imagesc(EixoX, EixoY, reliability(OFFSET_1:end,:));
    %    colormap(a); ylabel('Frequency(Hz)','FontSize', FonteLegenda);axis xy;
    %    axis([ 0 TotalTime 0 FundoEscala/2]);title('Reliability (R)','FontSize', FonteLegenda );

    %    subplot(2,2,2),imagesc(EixoX, EixoY, MascaraR(OFFSET_1:end,:));
    %    colormap(a); ylabel('Frequency(Hz)','FontSize', FonteLegenda);axis xy;
    %    axis([ 0 TotalTime 0 FundoEscala/2]);title('Reliability mask (M_{R})','FontSize', FonteLegenda );

    %    subplot(2,2,3),imagesc(EixoX, EixoY, newim(OFFSET_1:end,:));
    %    colormap(a); ylabel('Frequency(Hz)','FontSize', FonteLegenda);axis xy;
    %    axis([ 0 TotalTime 0 FundoEscala/2]); title('Filtered image','FontSize', FonteLegenda );

    %    subplot(2,2,4),imagesc(EixoX, EixoY, binim(OFFSET_1:end,:));
    %    colormap(a); ylabel('Frequency(Hz)','FontSize', FonteLegenda);axis xy;
    %    axis([ 0 TotalTime 0 FundoEscala/2]);title('Filter mask (M_{G})','FontSize', FonteLegenda );

    %   saveas(gca, 'fig3.eps', 'eps');

    %    plotridgeorient(orientim, 3, imlog,22);
    %    axis xy; xlabel('time (s)'); ylabel('Frequency'); 
    %    saveas(gca, 'orientacao.eps', 'eps');


