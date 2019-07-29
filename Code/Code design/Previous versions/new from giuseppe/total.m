%Run the aircraft flyover signal named "signal"
signal = signal.';

[s,f,t] = spectrogram(signal,800,[],[],40000);
T=t.';
F=f;
psd1 = ((t(2)-t(1))^2/t(end))*(abs(s).^2);
decibels = 10*log10(psd1/((2*10^(-5))^2));

%comment out this part if you don't want to print the spectrogram of the noise
figure,
surf(t,f,decibels,'Edgecolor','none')
view(2)
colormap jet
colorbar

