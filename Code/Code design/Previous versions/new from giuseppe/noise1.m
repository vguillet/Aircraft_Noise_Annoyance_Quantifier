%Run the background noise signal named "noise"
noise = noise.';

[s,f,t] = spectrogram(noise,800,[],[],40000);
T=t.';
F=f;
psd1noise = ((t(2)-t(1))^2/t(end))*(abs(s).^2);
decibelsnoise = 10*log10(psd1noise/((2*10^(-5))^2));

%comment out this part if you don't want to print the spectrogram of the noise 
figure,
surf(t,f,decibelsnoise,'Edgecolor','none')
view(2)
colormap jet
colorbar

