%Welcome to the "Sounds Correction and Analysis Software " by CD8&co.

%User Manual:
%{
This version is optimised for the first 13 datasets with reference numbers up to and
including reference_nb 13 (for testing data not included in, a manual input mode is built into the software).

For general use, follow the steps provided bellow:
In the same order:
    1)Look up in the "flyovers" excel spreadsheet the reference number
      matching the flyover to be investigated
    2)Enter the reference number if available in the program options and
      run after having selected other relevant options

    3)If no reference numbers are available, set the reference number to 0
      and manualy input the relevant variables in the manual input section
      after having imported the data to be analysed and renamed it to "signal"

WARNING: The sharpness metric only works on MATLAB R2018a and up! If you
         dont have matlab 2018a, go to advanced options and set sharpness 
         to false
%}

%==========================================================================
%Program options:
%==========================================================================

reference_nb = 1;

%{
Can be found in the excel spreadsheet "flyovers"
(automaticaly selects all the relevant parameters matching the data tested
in the program)

For manual setup of the initial input parameters, set reference_nb = 0
%}

%-----------------------------------------------------
%General options
%-----------------------------------------------------
%Make true the option you are interested in:
%Only select one option per run (only one true allowed for plots).

Plot_oaspl_ospl = false;

Plot_r = false;                 %Plot r versus time

Plot_Correction = true;         %Plot a spectogram of the data corrected
Plot_Save_Correction = false;   %Plot and save to source folder a spectogram of the data corrected

Plot_Tonality = false;          %Plot a graph of the tonality
Plot_Save_Tonality = false;     %Plot and save to source folder a graph of the tonality

Print_aspl = false;             %Print aspl on the screen

Print_Loudness = false;         %Prints phonresult and phonav

%-----------------------------------------------------
%Advanced options
%-----------------------------------------------------
%Change to false to switch off parts of the code
A_weighted = true;
Tonality = true;
Metrics = true;
Sharpness = true;

%Play the datafile tested as audio
output_sound = false;




%--------------------------------------------------------------------------
%Manual input mode (ignore if data tested reference_nb not equal to 0)
%--------------------------------------------------------------------------
%Fill in manualy the following:

h0= 0;      %initial altitude matching aircraft data set [m]
v0= 0;      %velocity matching aircraft data set [m/s]

manual_tonal1 = [0,0]; %Tonal components
manual_tonal2 = [0,0]; %Tonal components




%Run program!






%DO NOT MODIFY MAIN PROGRAM STRUCTURE WITHOUT KNOWING EXACTLY WHAT YOU ARE
%DOING!

%FUNCTIONS CAN BE EDITED, BUT ANY MODIFICATION TO THE OUTPUTS/INPUTS NEED
%TO BE REPORTED!

%==========================================================================
% Main Program structure:
%==========================================================================

load coredata.mat
disp("Dataset reference_nb selected:")
disp(reference_nb)

%-----------------------------------------------------
% Inputs:
%-----------------------------------------------------

dt = 1/40000;          %sampling rate
r0 = 1;                %virtual distance from the aircraft for geometrical spreading correction
c = 343;
theta = 0;             %input angle of decent matching aircraft data set

if reference_nb == 0
    disp ("Manual input mode selected")
else
    %aircraft data set selection based on reference number
    signallst = {signal1,signal2,signal3,signal4,signal5,signal6,signal7,signal8,signal9,signal10,signal11,signal12,signal13};
    signal = signallst{1,reference_nb};
    
    %initial altitude matching aircraft data set [m]
    h = [60.44,58.79,67.45,59.74,58.01,53.20,57.92,63.40,60.96,52.74,63.28,57.53,63.08];
    h0 = h(reference_nb);
    %velocity matching aircraft data set [m/s]
    v = [68.46,66.82,86.73,71.97,58.90,62.45,57.97,75.12,86.40,61.58,72.34,82.92,73.63];
    v0 = v(reference_nb);   
end

if output_sound == true
    sound(signal,40000)
end

bn2 = @bndouble;
bnsignal = bn2(nsignal);

%-----------------------------------------------------
% A-weighted metric computation:
%-----------------------------------------------------

if A_weighted == true
    al = @aweighted;
    [oaspl,ospl] = al(signal, Plot_oaspl_ospl);
end

%-----------------------------------------------------
% Corrections:
%-----------------------------------------------------

%Applying first background noise correction (frequency) to signal
f1 = @BNfc;
signal = f1(signal);

%Applying fast fourier transform to aircraft signal
f2 = @FFTs;
[psd1,T,F] = f2(signal);

%Applying fast fourier transform to background noise signal
f3 = @FFTbns;
bnpsd1 = f3(bnsignal);

%Applying second background noise correction (psd/bnpsd1)
f4 = @BNdc;
decibels = f4(psd1,bnpsd1);

%Applying doppler effect correction
f5 = @DSC;
[X,index,xq,yq,vq,dT] = f5(decibels,T,F,psd1,c,h0,v0);

%Applying Geometrical spreading effect correction
f6 = @GSc;
[vq,r] = f6 (index,v0,theta,xq,h0,dT,vq,r0,Plot_r);

%Applying Atmospheric correction:
f7 = @ATMc;
yq = f7(yq,vq,r);

%Perform final correction treatment for code compatibility with the metric section
f8 = @corrfin;
[time,frequency,power] = f8(xq,yq,vq,Plot_Correction,Plot_Save_Correction);

%-----------------------------------------------------
% Metrics:
%-----------------------------------------------------

if Tonality == true
%Tonality
    m1 = @tonality;
    [K,Tonal] = m1(power,frequency,time,reference_nb,Plot_Tonality,Plot_Save_Tonality,manual_tonal1,manual_tonal2);
end

if Metrics == true
%epnl calculation   
    m1 = @epnlm;
    epnl_value = m1(power,frequency,time,bands);
    
%aspl calculation
    m2 = @asplcalc;
    aspl = m2(power,frequency,index,Print_aspl);

%Loudness metric:
    m3 = @loudness;
    [phonresult,phonav] = m3(aspl,Print_Loudness);
end

if Sharpness == true
%Sharpness metric:
    m4 = @sharpness;
    S = m4(phonresult); 
end

disp("Signal analysis finished")








%==========================================================================
% Core code functions:
%==========================================================================

function bnsignal = bndouble(nsignal)

bnsignal = [nsignal,nsignal];

end

%-----------------------------------------------------
% Corrections:
%-----------------------------------------------------

%Background noise (frequency) correction:
function signal = BNfc(signal)

fc = 300;       %fc represents the 'cutoff frequency' and can be changed.
fs = 40000;
fcnorm = fc/(fs/2);
[b,a] = butter(4,fcnorm,'high');
% Obtaining corrected data

signal  = filtfilt(b,a,signal);

end

%Fast Fourier Transform (signal):
function [psd1,T,F] = FFTs(signal)

[s,f,t] = spectrogram(signal,800,[],[],40000);
T=t.';
F=f;
psd1 = ((t(2)-t(1))^2/t(end))*(abs(s).^2);

end

%Fast Fourier Transform (bnsignal):
function bnpsd1 = FFTbns(bnsignal)

[s,~,t] = spectrogram(bnsignal,800,[],[],40000);
bnpsd1 = ((t(2)-t(1))^2/t(end))*(abs(s).^2);

end

%Background noise (decibels) correction:

function decibels = BNdc(psd1,bnpsd1)

%size the noise powers matrix as the signal power matrix
bnpsd1_sized = bnpsd1(:,1:size(psd1,2));
%subtract the noise from the signal
polished = psd1-bnpsd1_sized;
%filter the low powers
polished(polished<0.00000001)=0;
%transform in decibels
decibels = 10*log10(polished/((2*10^(-5))^2));

end

%Doppler shift correction:
function [X,index,xq,yq,vq,dT] = DSC(decibels,T,F,psd1,c,h0,v0)

%-------------- CORRECTED FREQUENCIES ------------%
V=zeros(size(decibels,2),1);
dT=T(2)-T(1);
dF=F(2)-F(1);
Fmat=zeros(size(decibels,1),size(decibels,2));
ALPHA=zeros(size(decibels,2),1);
[~,index] = max(sum(psd1,1)); %Find maximum power at all times.
idxMAX= round(index-h0/c/0.01);%Approximation of passing time index.
%display (idxMAX)
%-------------- DE-DOPLERIZE --------------%

for j=1:1:size(decibels,1)
    t=0;
    i=1;
while t<=T(idxMAX)
    t = t+dT;
    DeltaT=T(idxMAX)-t;
    alpha=abs(atan(h0/(v0*DeltaT)));
    ALPHA(i)=alpha;
    dr = v0*dT*cos(alpha);
    V(i,1)=-dr/dT;
    Fmat(j,i)=1/(1-V(i)/343)*F(j);
    i=i+1;
end
while t<T(end)
    t = t+dT;
    DeltaT=t-T(idxMAX);
    alpha =abs(atan(h0/(v0*DeltaT)));
    ALPHA(i)=alpha;
    dr = v0*dT*sin(alpha);
    V(i,1)=dr/dT;
    Fmat(j,i)=1/(1-V(i)/343)*F(j);
    i=i+1;
end
end
a=size(F,1)*size(T,1);%number of triplets
X=ones(a,3);
k=1;
for i=1:1:size(T,1)
    for j=1:1:size(F,1)
        X(j+(i-1)*size(F,1),1)=T(k);
    end
k=k+1;
end
for i=1:1:size(T,1)
    for j=1:1:size(F,1)
        X(j+(i-1)*size(F,1),2)=Fmat(j,i);
    end
end
for i=1:1:size(T,1)
    for j=1:1:size(F,1)
        X(j+(i-1)*size(F,1),3)=decibels(j,i);
    end
end

%----------------------------------------------%

[xq,yq] = meshgrid(0:dT:T(end), 0:dF:20000); %Creates X and Y axes for the grid.
vq = griddata(X(:,1),X(:,2),X(:,3),xq,yq); %Relates Z to X and Y axes.

index=idxMAX;
end

%Geometrical spreading correction:
function [r,timelist,x] = rdet(index,v0,theta,xq,h0,dT,Plot_r)

timelist = zeros(1,size(xq,2));

for i=1:size(xq,2)
    tnew = xq(1,i);
    timelist(1,i) = tnew;
end

y=zeros(1,length(timelist));
x=zeros(1,length(timelist));
r=zeros(1,length(timelist));

%_____________________________v calculation:___________________________

vy=v0*sin(theta*(pi/180));         %speed positive down
vx=v0*cos(theta*(pi/180));

%_____________________________x calculation:_________________________

%find x=0 position
timeatx0 = index*dT;

%finding most negative x distance from microphone
xneg = -(timeatx0*vx);

%finding x distance from microphone for every datapoint
for i=1:length(timelist)

    xi=xneg+(((i-1)*dT)*vx);
    x(1,i) = abs(xi); %Appending all the x positions in list x
end

%_____________________________y calculation:___________________________

ymax=(vy*timeatx0)+h0;          %finding max height at start of descent

for i=1:length(timelist)

    yi=ymax-(((i-1)*dT)*vy);
    y(1,i)= yi;      %Appending all the y positions in list y
end

%____________________r magnitude calculation:______________________

%Computing r length using pythagoras
for i=1:length(timelist)

   ri=sqrt(x(i)*x(i)+y(i)*y(i));
   r(1,i) = ri;
end

%_________________________________________________
%PLOTTING
if Plot_r == true
    plot(timelist,r);
end

end
function vq = GS(vq,r,r0)

vqcorr = zeros(size(vq));

for i = 1:size(vq,2)
    for j = 1:size(vq,1)
        vqcorr(j,i) = vq(j,i) + 20*log10(r(i)/r0);
    end
end

vqcorr(vqcorr<60) = -inf;
vq = vqcorr;
end

function [vq,r] = GSc (index,v0,theta,xq,h0,dT,vq,r0,Plot_r)

%Calculating r(distance aircraft-mic) for every point
f6 = @rdet;
r = f6(index,v0,theta,xq,h0,dT,Plot_r);

%Applying geometrical spreading correction
f7 = @GS;
vq = f7(vq, r, r0);

end %The geometrical spreading correciton being done in two part, a function was created to combine the two steps (r determination and GS calculation)as one for ease of use

%Atmospheric correction:
function yq= ATMc(yq,vq,r)

T = 293;
h = 0.70;
f_ro = 24 + (4.04*(10^4)*h*(0.02+h)/(0.391+h));
f_rn = 9+280*h;

alpha = zeros(size(yq,1),1);
alpha_check = zeros(size(yq,1),1);
for k = 1:size(yq,1)
    term1 = (1.84*(10^-11));
    term2 = 0.01275*(exp(-2239.1/T))/(f_ro+((yq(k,1))^2)/f_ro);
    term3 = 0.1068*exp(-3352/T)/(f_rn+((yq(k,1))^2)/f_rn);
    alpha_check(k) = 869*((yq(k,1))^2)*(term1+term2+term3);
end

% replace alpha by alpha_check
alpha = alpha_check;
alpha = alpha/1000;

for f = 1:size(alpha)
    for g = 1:size(r,2)
        dist = r(g);
        atmeffect = alpha(f)*dist;
        vq(f,g) = vq(f,g) + atmeffect;
    end
end

end

%Final corrections treatment:
function [time,frequency,power] = corrfin(xq,yq,vq,Plot_Correction,Plot_Save_Correction)
%Renaming variables
time = xq;
frequency = yq;
power = vq;

%_________________________________________________
%PLOTTING

if Plot_Correction == true
    mesh(time,frequency,power);  %Creates surface.
    view(2);                     %Makes the view 2D.
    colormap jet;                %Selects color.
    colorbar;
    xlabel('Time [s]')
    ylabel('Frequency [Hz]')
    zlabel('Power [dB]')
    title('Corrected data spectrogram')
end

if Plot_Save_Correction == true
    fcorr = figure();
    mesh(time,frequency,power);  %Creates surface.
    view(2);                     %Makes the view 2D.
    colormap jet;                %Selects color.
    colorbar;
    xlabel('Time [s]')
    ylabel('Frequency [Hz]')
    zlabel('Power [dB]')
    title('Corrected data spectrogram')
    saveas(fcorr, 'Spectogram_corrected.png');
end

end


%-----------------------------------------------------
% Metrics:
%-----------------------------------------------------

%(Pre-correction metrics)
%A-weigted level metric:
function [oaspl,ospl] = aweighted(signal,Plot_oaspl_ospl)

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
if Plot_oaspl_ospl == true
time = t';
figure
plot(time,oaspl)
hold on
plot(time,ospl)
hold off
end

end


%(Post-correction metrics)
%EPNL metric:
function epnl_value = epnlm(power, frequency, time, bands)

avspltot = zeros(24,size(time,2));
for h = 1:size(time,2)
    
    pow = power(:,h);
    freq = frequency(:,h); 
    freqsorted = sort(freq);

    % make a matrix with frequency in the first column and power in the second
    % column, sorted so that frequency goes from smallest value to biggest,
    % keeping the appropriate power in the second column
    freqpowsorted = [freqsorted,zeros(length(freqsorted),1)];
    for n = 1:length(freqsorted)
        x = freqsorted(n,1);
        i = find(freq == x);
        freqpowsorted(n,2) = pow(i,1);
    end

    % delete rows where power = -inf
    % the matrix shortens after every row removal so the "r" is there to
    % account for it
    r = 0;
    for m = 1:length(freqpowsorted)
        if freqpowsorted((m-r),2) < 0 
            freqpowsorted((m-r),:)=[];
            r = r + 1;
        end
    end

    lowfreq = [44.7,56.2,70.8,89.1,112,141,178,224,282,355,446,562,708,891,1122,1413,1778,2239,2818,3548,4467,5623,7079,8913,11220];
    lowfreq = transpose(lowfreq);
    
    avspl = zeros(24,1);
    plist = zeros(length(freqpowsorted),(length(lowfreq)-1)); 
    % plist has 24 columns for every band
    % for every band, pressures will be calsulated for points which frequencies
    % lye in the band, where every p is a sum of the previous and the current p
    lenlist = zeros(24,1);
    % for every band, find the frequencies and calculate the average spl 
    for k = 1:(length(lowfreq)-1)

        p = 0;
        for t = 1:length(freqpowsorted)
            if freqpowsorted(t,1)>lowfreq(k,1) && freqpowsorted(t,1)<lowfreq((k+1),1)
                powerr = 10^(freqpowsorted(t,2)/10);
                p = p + sqrt(powerr);
                plist(t,k) = p;
            end
        end

        len = nnz(plist(:,k)); % number of non zero elements, this is the number of points that lye in this band     
        lenlist(k,1) = len;
        pav = max(plist(:,k))/len; % final pressure sum divided by number of points
        
        avspl(k,1) = 20*log10(pav); % average sound pressure level

    end
    
   
    avsplcol = avspl(:,1);
    avspltot(:,h) = avsplcol;
    
end

avspltot(isnan(avspltot)) = 0;
for row = 1:size(avspltot,2)
    for cell = 1:size(avspltot,1)
        if avspltot(cell,row) < 0
            avspltot(cell,row)=0;
        end
    end
end   

pnlt = zeros(size(time,2),1);
for b = 1:size(time,2)
    %________________________aspl to pnl_____________________________ 

    % aspl to noys 
    band = table2array(bands(:,2));
    L1 = table2array(bands(:,4));
    L2 = table2array(bands(:,6));
    L3 = table2array(bands(:,8));
    Lc = table2array(bands(:,10));
    L4 = table2array(bands(:,12));
    M1 = table2array(bands(:,5));
    M2 = table2array(bands(:,7));
    M3 = table2array(bands(:,9));
    M4 = table2array(bands(:,11));
  

    N_list = zeros((numel(band)-1),1);
    for n = 1:(numel(band)-1)    
        if L1(n,1)<avspltot(n,b)<L2(n,1)
            N = 0.1*(10^(M1(n,1)*(avspltot(n,b)-L1(n,1))));
            N_list(n,1) = N; 
        elseif L2(n,1)<avspltot(n,b)<L3(n,1) 
            N = 10^(M2(n,1)*(avspltot(n,b)-L3(n,1)));
            N_list(n,1) = N; 
        elseif L3(n,1)<avspltot(n,b)<Lc(n,1) 
            N = 10^(M3(n,1)*(avspltot(n,b)-L3(n,1)));
            N_list(1,n) = N; 
        elseif Lc(n,1)<avspltot(n,b)<150 
            N = 10^(M4(n,1)*(avspltot(n,b)-L4(n,1)));
            N_list(n,1) = N; 
        end
    end

    % summing the noys   
    Nmax = max(N_list);
    N_list(isnan(N_list)) = 0; 
    Nsum = sum(N_list);
    Ntot = 0.85*Nmax + 0.15*Nsum;

    % total number of noys to pnl 
    pnltot = 40 + 10/(log10(2))*(log10(Ntot));

    %______________________ tone corrected pnl _________________________

    % step 1 - computing arithmetic difference in sound pressure level
   
    dspl = zeros(numel(band),1);
    apr = zeros(numel(band),1);
    dpr = zeros(numel(band),1);
    
    for n = 4:(numel(band)-1)
        apr(n,1) = 10^(avspltot(n,b)/20);   % apr - average pressure
    end
    
    r = 0;
    nonzerolist = apr;
    for w = 1:length(nonzerolist)
        if nonzerolist(w-r) == 1 || nonzerolist(w-r) == 0
            nonzerolist(w-r) = [];
            r = r + 1;
        end
    end

    dpr = zeros(length(nonzerolist),1);
    for cell = 2:length(nonzerolist)
        dpr(cell,1) = abs(nonzerolist(cell) - nonzerolist(cell-1));
    end
    
    i = 0;
    for r = 1:length(apr)
        if apr(r,1) ~= 0 && apr(r,1) ~= 1
            i = i + 1;
            dspl(r,1) = 20*log10(dpr(i,1)); % difference in sound pressure levels
            if i == 1
                dspl(r,1) = 0;
            end
        end
    end
    
    
    % step 2 - values of dspl bigger than 5
    % step 3 - choosing aspls
    column = avspltot(:,b);
    aspl = [column,zeros(length(column),1)]; % adding a column of zeros to aspl list

    for n = 5:numel(band)
        if abs(dspl(n,1)-dspl((n-1),1))> 5 
            if dspl(n,1)>0 && dspl(n,1)>dspl((n-1),1)
                aspl(n,2) = 1;
            elseif dspl(n,1)<=0 && dspl((n-1),1)>0
                aspl((n-1),2) = 1;
            end
        end
    end

    % step 4 - new adjusted spl

    aspl = [aspl,zeros(length(aspl),1)];

    for n = 3:(numel(band)-2)
        if aspl(n,2) == 1 
            aspl(n,3) = 0.5*(aspl((n-1),1)) + aspl((n+1),1);
        else 
            aspl(n,3) = aspl(n,1);
        end
    end

    for n = (numel(band)-1)
        if aspl(24,2) == 1
            aspl(24,3) = aspl(23,1) + dspl(23,1);
        else 
            aspl(24,3) = aspl(24,1);
        end
    end

    % step 5 - add dspl(3) and dspl(25)
    
    dspl1 = zeros(numel(band),1);
    pres = zeros(numel(band),1);
    dpres = zeros(numel(band),1);
    for n = 4:(numel(band)-1)
        pres(n,1) = 10^((aspl(n,3))/20);  % pres - pressure
        dpres(n,1) = abs(pres(n,1)-pres((n-1),1));   % arithmetic difference
        dspl1(n,1) = 20*log10(dpres(n,1));  % difference in sound pressure levels
    end
    
    for n = 3
        dspl1(n,1) = dspl1((n+1),1);
    end

    for n = length(band)
        dspl1(n,1) = dspl1((n-1),1);
    end
    
    for cell = 1:length(dspl1)
        if dspl1(cell,1) < 0
            dspl1(cell,1) = 0;
        end
    end
    
    % step 6 - arithmetic average of 3 adjacent slopes

    dspl3 = zeros(numel(band),1);
    pnew = zeros(numel(band),1);
    p3 = zeros(numel(band),1);
    for n = 3:23
        pnew(n,1) = 10^(dspl1(n,1)/20);                     % pnew - new pressure
        if pnew(n,1) == 1 || pnew((n+1),1) == 1 || pnew((n+2),1) == 1
            p3(n,1) = 0;
        else
            p3(n,1) = (pnew(n,1) + pnew((n+1),1) + pnew((n+2),1))/3;    % arithmetic average
            dspl3(n,1) = 20*log10(p3(n,1));                     % logarithmic average
        end
    end

    % step 7 - final 1/3 octave band spls

    aspl = [aspl,zeros(length(aspl),1)];

    for n = 3
       aspl(n,4) = aspl(n,1);
    end

    for n = 4:23
        aspl(n,4) = aspl((n-1),4) + dspl3((n-1),1);
    end

    % step 8 - differences between aspl(1,n) and aspl(4,n)

    F = zeros(numel(band),1);
    for n = 3:24
       F(n,1) = aspl(n,1) - aspl(n,4);
    end

    % step 9 - tone correction factors

    C = zeros(numel(band),1);
    for n = 3:11
       if 1.5 <= F(n,1) < 3 
          C(n,1) = F(n,1)/3 - 0.5;
       elseif 3 <= F(n,1) < 20 
          C(n,1) = F(n,1)/6;
       elseif 20 <= F(n,1) 
          C(n,1) = 10/3;
       end
    end

    for n = 12:21
       if 1.5 <= F(n,1) < 3 
          C(n,1) = F(n,1)/3*2 - 1;
       elseif 3 <= F(n,1) < 20
          C(n,1) = F(n,1)/3;
       elseif F(n,1) >= 20 
          C(n,1) = 20/3;
       end
    end

    for n = 22:24
       if 1.5 <= F(n,1) < 3
          C(n,1) = F(n,1)/3 - 0.5;
       elseif 3 <= F(n,1) < 20
          C(n,1) = F(n,1)/6;
       elseif F(n,1) >= 20
          C(n,1) = 10/3;
       end
    end

    % step 10 - tone corrected pnl
    
    pnlt(b,1) = pnltot + max(C(:,1));
    
end

antilog_pnlt = zeros(length(pnlt),1);
for g = 1:length(pnlt)
    antilog_pnlt = 10^(pnlt(g,1)-2);
end
antilog_pnlt_sum = sum(antilog_pnlt);
t1 = 0;
t2 = time(1,size(time,2));
epnl_value = log10(antilog_pnlt_sum/length(pnlt))-13;

end

%Tonality metric:
function [K,Tonal] = tonality(power,frequency,time,reference_nb,Plot_Tonality,Plot_Save_Tonality,manual_tonal1,manual_tonal2)

%2017-10-17_11-00-03
freqtone = [1060, 1175, 2875, 3360, 3910, 4060, 5080, 5250, 6325, 7425, 9645, 10150, 10320, 10940, 15600];
freqtonb = [940, 1075, 2775, 3245, 3790, 3930, 4950, 5125, 6225, 7275, 9540, 10050, 10220, 10820, 15500];

% Tonal components
if reference_nb == 0
    freqtone = manual_tonal1;
    freqtonb = manual_tonal2;
elseif reference_nb == 1
    % 2017-08-14_13-13-48
    freqtone = [449 546 860 1133 1406 1758 2109 2773 3125 3398 3555 3789 3945 4141 4297 4531 4766 4883 5000 5117 5312 5664 5977 6367 6445 7305];
    freqtonb = [352 449 703 976 1289 1602 1992 2578 3008 3281 3438 3672 3828 4023 4129 4453 4648 4766 4883 5039 5195 5547 5859 6250 6328 7188];

elseif reference_nb == 2
    %2017-08-14_13-15-16
    freqtone = [429 507 781 930 1094 1250 1386 1484 1641 1914 2636 2734 2891 3125 3320 3555 3633 3828 4004 4102 4297 4492 4727 5195 5508 6328 6641 6836 7305 7519 7617 7773 7988 8086 8359 8672 8984 9648];
    freqtonb = [312 391 625 820 977 1133 1289 1386 1523 1797 2539 2636 2773 2969 3203 3398 3477 3672 3906 4004 4141 4336 4570 5078 5352 6250 6523 6719 7148 7422 7519 7656 7891 7988 8203 8555 8828 9531];

elseif reference_nb == 3
    %2017-08-14_13-16-48
    freqtone = [468 859 1133 1445 2109 2383 2815 3164 3320 4062 4492 4688 4844 5000 5195 5508 5664 6055 6367 6719 7031 7387];
    freqtonb = [351 742 977 1328 1992 2266 2659 3047 3242 3906 4375 4570 4727 4883 5078 5391 5547 5898 6250 6641 6875 7266];

elseif reference_nb == 4
    %2017-08-14_13-19-05
    freqtone = [312 430 703 937 1172 1445 1680 1914 2070 2266 2656 2891 3125 3359 3594 3750 3906 4141 4375 4570 4727 4961 5195 5390 5781 6094 6289 6406 6601 6914 7148 7422 7695 7851 8164 8320 8476 9140 9336];
    freqtonb = [234 351 547 742 977 1211 1562 1758 1953 2148 2500 2734 2969 3242 3438 3633 3789 3984 4219 4414 4609 4844 5078 5273 5625 5938 6172 6289 6445 6797 7031 7305 7578 7734 8047 8164 8359 8984 9180];

elseif reference_nb == 5
    %2017-08-14_13-22-04
    freqtone = [546 6015 976 1211 1484 1836 3047 3242 3476 3672 3945 4140 4375 4531 4765 4922 5117 5468 5586 5781 6133 6250 6445 6679 7070 7305 7695 8008 8164 8437 9258];
    freqtonb = [390 5859 820 1094 1328 1719 2930 3125 3320 3555 3828 4023 4258 4375 4648 4805 4961 5312 5469 5664 6016 6133 6289 6562 6953 7188 7578 7891 8008 8281 9141];

elseif reference_nb == 6
    %2017-08-14_13-23-36
    freqtone = [9648 9023 8711 8437 8281 8125 7656 7031 6836 6601 6328 6015 5742 5547 5195 4883 4687 4336 4719 3984 3320 3086 1914 1679 1445 1054 742 507];
    freqtonb = [9492 8906 8594 8320 8164 8008 7500 6875 6719 6445 6172 5898 5625 5391 5078 4727 4570 4219 4602 3828 3203 2969 1758 1523 1289 937 586 351];

elseif reference_nb == 7
    %2017-08-14_13-25-04
    freqtone = [507 742 900 1055 1406 1953 2461 2734 3008 3203 3430 3750 3945 4336 4570 4648 4844 4961 5430 5820 6133 6758 7070 7344 7773 8393 9219 10120 10550 11310 11950 12770 13710 14060 15630];
    freqtonb = [390 625 800 900 1290 1797 2344 2617 2891 3047 3281 3516 3789 4219 4414 4570 4648 4844 5273 5623 5938 6602 6914 7227 7617 8281 9102 10000 10430 11210 11880 12660 13630 13980 15510];
    
elseif reference_nb == 8
    %2017-08-14_13-28-58
    freqtone = [650, 800, 950, 1550, 1900, 2650, 3200, 3675, 4775, 5000, 5625, 5900, 6300, 6950, 7700, 8400, 9100, 9850, 10350, 10580, 10770, 12080, 12580];
    freqtonb = [550, 700, 850, 1400, 1800, 2550, 3090, 3560, 4675, 4880, 5525, 5750, 6200, 6830, 7580, 8300, 8970, 9750, 10250, 10470, 10670, 11920, 12460];

elseif reference_nb == 9
    %2017-10-17_10-46-21
    freqtone = [420, 703, 937, 1133, 1211, 1406, 1523, 1680, 1836, 2617, 2800, 3086, 3350, 3550, 3775, 3900, 4125, 4300, 4650, 4850, 5120, 5400, 5800, 6025, 6200, 7150, 8750, 9350, 10250, 11100, 11580];
    freqtonb = [350, 547, 820, 976, 1172, 1289, 1445, 1602, 1758, 2539, 2700, 3008, 3200, 3440, 3675, 3775, 4025, 4175, 4525, 4725, 5000, 5250, 5700, 5900, 6100, 7050, 8625, 9250, 10100, 11000, 11450];

elseif reference_nb == 10
    %2017-10-17_10-50-19
    freqtone = [450, 600, 750, 900, 1400, 1800, 2500, 3000, 3250, 3350, 3800, 4200, 4450, 4850, 5100, 5450, 5725, 5950, 6300, 6450, 6850, 7050, 7350, 7700, 8100, 8350, 8800, 9525, 10150, 10300, 10550, 11275, 12125, 12750, 13150, 14150, 14280, 14700];
    freqtonb = [300, 500, 650, 800, 1300, 1700, 2400, 2900, 3100, 3250, 3675, 4100, 4300, 4700, 5000, 5300, 5625, 5850, 6200, 6350, 6700, 6950, 7200, 7550, 8000, 8250, 8650, 9425, 10050, 10200, 10450, 11225, 11975, 12600, 13000, 14050, 14180, 14600];

elseif reference_nb == 11
    %2017-10-17_10-54-49
    freqtone = [570, 775, 1040, 1300, 2425, 2600, 3125, 4250, 4500, 4900, 5780, 5930, 7100, 7400, 8950, 10600, 10800];
    freqtonb = [425, 625, 940, 1170, 2300, 2500, 3000, 4100, 4350, 4750, 5650, 5780, 6980, 7250, 8850, 10500, 10650];

elseif reference_nb == 12
    %2017-10-17_10-58-10
    freqtone = [475, 1270, 1425, 1725, 1875, 2325, 2575, 2910, 3300, 3450, 3975, 4335, 4520, 5575, 5825, 6825, 7125, 8100, 8300, 9025, 9175, 9650, 9990, 10550, 10800, 11130, 11900];
    freqtonb = [350, 1125, 1325, 1625, 1725, 2225, 2450, 2780, 3200, 3350, 3850, 4200, 4420, 5425, 5725, 6675, 7025, 7950, 8250, 8925, 9075, 9550, 9870, 10450, 10650, 11000, 11750];

elseif reference_nb == 13
    %2017-10-17_11-00-03
    freqtone = [475, 1270, 1425, 1725, 1875, 2325, 2575, 2910, 3300, 3450, 3975, 4335, 4520, 5575, 5825, 6825, 7125, 8100, 8300, 9025, 9175, 9650, 9990, 10550, 10800, 11130, 11900];
    freqtonb = [350, 1125, 1325, 1625, 1725, 2225, 2450, 2780, 3200, 3350, 3850, 4200, 4420, 5425, 5725, 6675, 7025, 7950, 8250, 8925, 9075, 9550, 9870, 10450, 10650, 11000, 11750];

end

numfreqtone = round(freqtone/39.0625)+1;
numfreqtonb = round(freqtonb/39.0625)+1;
freqtone = frequency(numfreqtone,1);
freqtonb = frequency(numfreqtonb,1);
freqton = (freqtone+freqtonb)/2;
numfreqton = round(freqton/39.0625)+1;

power(~isfinite(power))=0;
power(isnan(power))=0;

% Filtering tonal components out of power / not used%

freqtonttab = zeros(length(freqton),length(freqton));
powton = [];
powtont = [];
powtonttot = [];
numfreqtont = zeros(length(freqton),length(freqton));

for i = 1:length(freqton)
powton = [powton; power(numfreqton(i),:)];
freqtont = [freqtonb(i):39.0625:freqtone(i)];
powtont = [];
    for j = 1:length(freqtont)
    freqtonttab(i,j) = freqtont(j);
    numfreqtont(i,j)= freqtont(j) /39.0625+1;
    powtont = [powtont; power(numfreqtont(i,j),:)];
    end
powtonttot = [powtonttot; sum(powtont)];
end

powtonttottot = sum(powtonttot,'omitnan');
powtot = sum(power,'omitnan');
powres = powtot-powtonttottot;


%Calculation EHS

Ehs_tab = 3.64*(freqton/1000).^(-0.8)-6.5*exp(-0.6*(freqton/1000-3.3).^2)+10^(-3)*(freqton/1000).^4;

for i = 1:length(Ehs_tab)
    if Ehs_tab(i) < 0
        Ehs_tab(i) = 0;
    end
end


% Calculation z

zall = 13*atan(0.76*frequency(:,1)/1000)+3.5*atan(frequency(:,1)/7500).^2;
z= 13*atan(0.76*freqton/1000)+3.5*atan(freqton/7500).^2;
zb = 13*atan(0.76*freqtonb/1000)+3.5*atan(freqtonb/7500).^2;
ze = 13*atan(0.76*freqtone/1000)+3.5*atan(freqtone/7500).^2;

%Calculation dz

dz = ze-zb;
dztab = [];
for i =1:length(power(1,:))
    dztab = [dztab dz];
end
dztab = transpose(dztab);

% Calculation EGr

zbb = zb-0.5;
zee = ze+0.5;
freqz = zeros(length(freqton),4);

for j = 1:length(zbb)
for i = 1:length(zall)
    if zbb(j) > zall(i)
        freqz(j,1) = frequency(i,1);
        freqz(j,3) = i;
    end
end
end

for j = 1:length(zee)
for i = 1:length(zall)
    if zee(j) > zall(i)
        freqz(j,2) = frequency(i,1);
        freqz(j,4) = i;
    end
end
end
powerz = zeros(1,length(power(1,:)));

EGr = [];
for i =1:length(zee)
for j =freqz(i,3):freqz(i,4)
    powerz = powerz + power(j,:);
    if j == freqz(i,4)
     powerz = powerz/length(freqz(i,3):freqz(i,4));
    end
end
EGr = [EGr,transpose(powerz)];
end
EGr = transpose(EGr);


% Calculation A_Ek and dLi

L_Ek = zeros(length(freqton),length(freqton));
A_Ek = zeros(length(freqton),length(freqton));
dLirl = zeros(length(power(1,:)),length(freqton));
freqtonrl = zeros(length(power(1,:)),length(freqton));
Atot = zeros(length(freqton),length(power(1,:)));
L_Ektab = [];
A_Ektab = [];

for a = 1:length(powton(1,:))
 for i = 1:length(freqton)
     Atot(i,a) = 0;
   for k = 1:length(freqton)
       if freqton(i) < freqton(k)
           s = 27;
           L_Ek(i,k) = powton(k,a) - s*(z(k)-z(i));
           A_Ek(i,k) = 10^(L_Ek(i,k)/20);
           Atot(i,a) = Atot(i,a) + A_Ek(i,k);
       elseif freqton(i) > freqton(k)
           s = -24-(230/freqton(k)+0.2*(powton(k,a)));
           L_Ek(i,k) = powton(k,a) - s*(z(k)-z(i));
           A_Ek(i,k) = 10^(L_Ek(i,k)/20);
           Atot(i,a) = Atot(i,a) + A_Ek(i,k);
       end
   end
   dLi = powton(i,a)-10*log(Atot(i,a)^2+Ehs_tab(i)+EGr(i,a));
   if dLi < 0
       dLi = 0;
   end
   dLirl(a,i) = dLi;
   freqtonrl(a,i) = freqton(i);
 end
   L_Ektab = [L_Ektab;L_Ek];
   A_Ektab = [A_Ektab;A_Ek];
end

%Calculate al w's%

hold off;

w3 = (1-exp(-(dLirl/15))).^0.29;
w2 = ((1+0.2*(freqtonrl*1/700+700./freqtonrl).^2).^(0.5)).^(-.29);
w1 = 0.13./(0.13+dztab);
w = (w1.^(1/0.29).*w2.^(1/0.29).*w3.^(1/0.29)).^2;
wT = sum(w,2).^(0.5);

%wGr = transpose(1-sum(power2,1,'omitnan')./sum(powtot,1,'omitnan'));
%wGr = transpose((1-powres./powtot));
wGr = 1;

% Calculate tonality
K = 1.09*wT.^0.29.*wGr.^0.79;

Tonal = median(K);


%_________________________________________________
%PLOTTING

if Plot_Tonality == true
    plot(time(1,:),K);
    ylim([0 0.5]);
    xlabel('Time [s]')
    ylabel('Tonality')
    title('Tonality vs Time')
end

if Plot_Save_Tonality == true
    ftonality = figure();
    plot(time(1,:),K);
    ylim([0 0.5]);
    xlabel('Time [s]')
    ylabel('Tonality')
    title('Tonality vs Time')
    saveas(ftonality,'Tonality_plot.png');
end
end

%aspl calculation
function aspl = asplcalc(power,frequency,index,Print_aspl)
%________average sound pressure level at t=index (for loudness)________

pow = power(:,index);
freq = frequency(:,index);
freqsorted = sort(freq);

% make a matrix with frequency in the first column and power in the second
% column, sorted so that frequency goes from smallest value to biggest,
% keeping the appropriate power in the second column
freqpowsorted = [freqsorted,zeros(length(freqsorted),1)];
for n = 1:length(freqsorted)
    x = freqsorted(n,1);
    i = find(freq == x);
    freqpowsorted(n,2) = pow(i,1);
end

% delete rows where power = -inf
% the matrix shortens after every row removal so the "r" is there to
% account for it
r = 0;
for m = 1:513
    if freqpowsorted((m-r),2) < 0
        freqpowsorted((m-r),:)=[];
        r = r + 1;
    end
end

lowfreq = [44.7,56.2,70.8,89.1,112,141,178,224,282,355,446,562,708,891,1122,1413,1778,2239,2818,3548,4467,5623,7079,8913,11220];
lowfreq = transpose(lowfreq);
avspl = zeros(24,1);
plist = zeros(length(freqpowsorted),(length(lowfreq)-1));
% plist has 24 columns for every band
% for every band, pressures will be calsulated for points which frequencies
% lye in the band, where every p is a sum of the previous and the current p
lenlist = zeros(24,1);
% for every band, find the frequencies and calculate the average spl
for k = 1:(length(lowfreq)-1)

    p = 0;
    for t = 1:length(freqpowsorted)
        if freqpowsorted(t,1)>lowfreq(k,1) && freqpowsorted(t,1)<lowfreq((k+1),1)
            p = p + 0.00002*sqrt((freqpowsorted(t,2))/(10^(-12)));
            plist(t,k) = p;
        end
    end

    len = nnz(plist(:,k)); % number of non zero elements, this is the number of points that lye in this band
    lenlist(k,1) = len;
    pav = max(plist(:,k))/len; % final pressure sum divided by number of points
    avspl(k,1) = 20*log10(pav/0.00002); % average sound pressure level

end
aspl = avspl;

if Print_aspl == true
    disp (aspl)
end

end

%Loudness metric:
function [phonresult,phonav] = loudness(aspl,Print_Loudness)

% input data points
f = [20 25 31.5 40 50 63 80 100 125 160 200 250 315 400 500 630 800 ...
     1000 1250 1600 2000 2500 3150 4000 5000 6300 8000 10000 12500];

af = [0.532 0.506 0.480 0.455 0.432 0.409 0.387 0.367 0.349 0.330 0.315 ...
      0.301 0.288 0.276 0.267 0.259 0.253 0.250 0.246 0.244 0.243 0.243 ...
      0.243 0.242 0.242 0.245 0.254 0.271 0.301];

Lu = [-31.6 -27.2 -23.0 -19.1 -15.9 -13.0 -10.3 -8.1 -6.2 -4.5 -3.1 ...
       -2.0  -1.1  -0.4   0.0   0.3   0.5   0.0 -2.7 -4.1 -1.0  1.7 ...
        2.5   1.2  -2.1  -7.1 -11.2 -10.7  -3.1];

Tf = [ 78.5  68.7  59.5  51.1  44.0  37.5  31.5  26.5  22.1  17.9  14.4 ...
       11.4   8.6   6.2   4.4   3.0   2.2   2.4   3.5   1.7  -1.3  -4.2 ...
       -6.0  -5.4  -1.5   6.0  12.6  13.9  12.3];

%Error Trapping
phonvalues = 0:200;
phontable = zeros(200,29);
for i = 1:2000
    k = i/10;
    phon = k;
    if((phon < 0) || (phon > 200))
        disp('Phon value out of bounds!')
        spl = 0;
        freq = 0;
    else
        %Setup user-defined values for equation
        Ln = phon;

        %Deriving sound pressure level from loudness level (iso226 sect 4.1)
        Af=4.47E-3 * (10.^(0.025*Ln) - 1.15) + (0.4*10.^(((Tf+Lu)/10)-9 )).^af;
        Lp=((10./af).*log10(Af)) - Lu + 94;

        %Return user data
        spl = Lp;
        freq = f;
        for j = 1:29
            phontable(i,j) = spl(1,j);
        end
    end

end

% compare Aspl to Phonlines
phonresult = zeros(24,1);
 for i = 1:24
     for j = 1:2000
         if aspl(i) <= phontable(j,(i+4))
             phonresult(i,1) = j/10;
             break
         end

     end
 end

 %append frequencies to phonresult
 for i = 1:24
     phonresult(i,2) = f(i+4);
 end

rr = 0;
phontot = 0;
for i = 1:24
    if phonresult(i,1) ~= 0
        rr = rr + 1;
        phontot = phontot + phonresult(i,1);
    end
end
phonav = phontot/rr;

if Print_Loudness == true
    disp(phonresult)
    disp(phonav)
end

end

%Sharpness metric:
function S = sharpness(phonresult)

r=0;
for m = 1:length(phonresult)
    if phonresult((m-r),1) == 0
       phonresult((m-r),:)=[];
       r = r + 1;
    end
end

miguel = phonresult;

%defining g(z) function
syms g(z)
g(z)= piecewise(z<=16, 1, z>16, 0.066*exp(0.171*z));

%defining the z that goes into the g(z) function
z_values_tab=zeros(length(miguel),1);
g_values_tab=zeros(length(miguel),1);
N_prime_tab=zeros(length(miguel),1);
%N_tab=zeros(length(miguel),1);
top_tab=zeros(length(miguel),1);
c=0.11;

for i=1:length(miguel)
    z_value=13*atan(0.76*miguel(i,2)/1000)+3.5*atan(miguel(i,2)/7500).^2;
    z_values_tab(i,1)=z_value;

    %calculating a g(z), and N_prime function for every z_value
    g_value=vpa(g(z_value));
    g_values_tab(i,1)=g_value;

    %define Loudness function HERE:
    N_prime= miguel(i,1);
    N_prime_tab(i,1)=N_prime;

    %defining N (bottom of sharpness function): integrate N'(z) between 0 and 24
    N = trapz(z_value,N_prime_tab);

    %defining the top of sharpness function
    top= g_values_tab(i)*N_prime_tab(i)*z_values_tab(i);
    top_tab(i,1)=top;
    a = trapz(z_value,top_tab);

end

%sharpness function
S = ((a/N)*c);
%plot(z_values_tab,top_tab)
%plot(z_values_tab,N_prime_tab)
end
