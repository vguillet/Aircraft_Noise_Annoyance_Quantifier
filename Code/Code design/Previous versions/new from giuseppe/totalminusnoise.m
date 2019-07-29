%size the noise powers matrix as the signal power matrix
psd1noise_sized = psd1noise(:,1:size(psd1,2));
%subtract the noise from the signal
polished = psd1-psd1noise_sized;
%filter the low powers
polished(polished<0.00000001)=0;
%transform in decibels
decibels_polished = 10*log10(polished/((2*10^(-5))^2));

%complete spectrogram before de-dopplerization
figure,
surf(t,f,decibels_polished,'Edgecolor','none')
view(2)
colormap jet
colorbar
