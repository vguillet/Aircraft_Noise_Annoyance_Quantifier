%------------------- Fast Fourier Transform of raw data -----------------%
[s,f,t] = spectrogram(signal,800,[],[],40000);
psd = ((t(2)-t(1))^2/t(end))*(abs(s).^2);
%decibels = 10*log10(psd/((2*10^(-5))^2));

%1/3 octave bands boundaries
b = [0,44.7,56.2,70.8,89.1,112,141,178,224,282,355,446,562,708,891,1122,1413,1778,2239,2818,3548,4467,5623,7079,8913,11220];

%------------------ Average sound pressure level ------------------------%
apsd = zeros(24,size(t,2)); % avg power in W
spl = zeros(24,size(t,2));  % avg spl in dB 
for i = 1:size(t,2)
    j = 1;
    for k = 1:(length(b)-1)
        count = 0;
        band = 0;
        while (f(j)>=b(k)) && (f(j)<b(k+1))
            band = band + psd(j,i);
            j = j + 1;
            count = count + 1;
        end
        avgband = band/count;
        apsd(k,i) = avgband;
        spl(k,i) = 10*log10(apsd(k,i)/((2*10^(-5))^2));
    end
end

spl(1,:) = []; %cuts off the frequencies < lowest band boundary


%--------------------- A-level weighted matrix --------------------------%
%central frequencies of each 1/3 octave band
centralf = [50.45, 63.5, 79.95, 100.55, 126.5, 159.5, 201, 253, 318.5, 400.5, 504, 635, 799.5, 1006.5, 1267.5, 1595.5, 2008.5, 2528.5, 3183, 4007.5, 5045, 6351, 7996, 10066.5]';

dLa = zeros(size(centralf));    % A-weighting correction
for m = 1:length(centralf)
    dLa(m) = -145.528 + 98.262*log10(centralf(m)) - 19.509*(log10(centralf(m)))^2 + 0.975*(log10(centralf(m)))^3;
end

La = zeros(size(spl));          % A-weighting corrected spl
oaspl = zeros(length(t),1);     % Overall A-weighted sound pressure level
for p = 1:size(t,2)
    sumLa = 0;
    for q = 1:size(centralf)
        La(q,p) = spl(q,p) + dLa(q,:);
        if ~isnan(spl(q,p))
            sumLa = sumLa + 10^(La(q,p)/10); % gives the summation inside the oaspl formula
        end
    end
    oaspl(p,:) = 10*log10(sumLa);
end

ospl = zeros(length(t),1);      % Overall sound pressure level uncorrected
for p = 1:size(t,2)
    sumos = 0;
    for q = 1:size(centralf)
        if ~isnan(spl(q,p))
            sumos = sumos + 10^(spl(q,p)/10);
        end
    end
    ospl(p,:) = 10*log10(sumos);
end


%plotting
time = t';
figure
plot(time,oaspl)
hold on
plot(time,ospl)
hold off

