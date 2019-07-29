%Welcome to the sounds correction software by CD8&co.
%To use, follow the steps provided bellow:
%!!!in the same order!!!:
    %1)Import background noise data, rename it to "bnsignal"
    %2)Import pressure data to be tested (should be named signal automaticaly).
    %3)Change h0,v0,theta etc... to match tested dataset
    %4)Select in the Tonal metric definition the matching Tonal components
    %5)Select what plot should be shown/ploted
    %4)Run and be patient, the program takes about 20~30s to generate diagram

%==========================================================================
%Inputs:
%==========================================================================

%signal = signal.';
%bnsignal = bnsignal.';

c = 343;
h0 = 63.4;             %input initial altitude matching aircraft data set
v0 = 75.12;            %input velocity matching aircraft data set 
theta = 3;             %input angle of decent matching aircraft data set          
dt = 1/40000;
r0 = 1;

%==========================================================================
%Program options:
%==========================================================================

%Make true the plots you are interested in:
Plot_Correction = false;
Save_Correction = false;

Plot_Tonality = false;
Save_Tonality = false;

%==========================================================================
% Main code:
%==========================================================================

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
%f6 = @GSc;
%vq = f6 (index,v0,theta,xq,h0,dT,vq,r0);

%Perform final correction treatment for code compatibility
f7 = @corrfin;
[time,frequency,power] = f7(xq,yq,vq,Plot_Correction,Save_Correction);

%-----------------------------------------------------
% Metrics:
%-----------------------------------------------------

%Tonality
m1 = @tonality;
K = m1(power,frequency,time,Plot_Tonality,Save_Tonality);

%Loudness metric:
%m2 = @loudness;
%phonresult = m2(aspl);

%Sharpness metric:
%m3 = @sharpness;
%S = m3(phonresult);


%==========================================================================
% Core functions:
%==========================================================================

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
function [r,timelist,x] = rdet(index,v0,theta,xq,h0,dT)

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
%plot(timelist,r)

end
function vq = GS(vq,r,r0)

vq = zeros(size(vq));

for i = 1:size(vq,2)
    for j = 1:size(vq,1)
        vq(j,i) = vq(j,i) + 20*log10(r(i)/r0);
    end
end

vq(vq<60) = -inf;
end

function vq = GSc (index,v0,theta,xq,h0,dT,vq,r0)

%Calculating r(distance aircraft-mic) for every point
f6 = @rdet;
r = f6(index,v0,theta,xq,h0,dT);

%Applying geometrical spreading correction
f7 = @GS;
vq = f7(vq, r, r0);

end %The geometrical spreading correciton being done in two part, a function was created to combine the two steps (r determination and GS calculation)as one for ease of use

%Final corrections treatment:
function [time,frequency,power] = corrfin(xq,yq,vq,Plot_Correction,Save_Correction)
%Renaming variables
time = xq;
frequency = yq;
power = vq;

%_________________________________________________
%PLOTTING
f = figure();

if Plot_Correction == true
    mesh(time,frequency,power);  %Creates surface.
    view(2);                     %Makes the view 2D.
    colormap jet;                %Selects color.
    colorbar;
end

if Save_Correction == true
    saveas(f, 'Spectogram_corrected.png');
end

end

%-----------------------------------------------------
% Metrics:
%-----------------------------------------------------

%Tonality metric:
function K = tonality(power,frequency,time,Plot_Tonality,Save_Tonality)

% Tonal components

% % 2017-08-16_11-16-07
% freqtone = [470 600 700 900 1100 1550 1750 1900 2300 2500 2850 3000 3350 3700 4250 4400 5200 5600 6450 6650 7900 8200 9250 10550]; 
% freqtonb = [300 500 600 750 1000 1400 1650 1800 2200 2350 2750 2850 3200 3550 4100 4250 5050 5450 6350 6550 7800 8100 9100 10400];
%
% % 2017-08-14_13-13-48
% freqtone = [400 500 600 800 1000 1580 2100 2640 3520 3960 5560 6367 6600 7109]; 
% freqtonb = [240 400 500 700 800 1500 2000 2600 3400 3800 5480 6300 6520 7070];
% 
% %2017-08-14_13-15-16
% freqtone = [1211 1406 1750 2188 2550 2800 3400 5078 5625 6650]; 
% freqtonb = [1055 1250 1641 2070 2450 2695 3300 4961 5500 6550];
% 
% %2017-08-14_13-16-48
% freqtone = [781 900 1250 1700 2150 2500 2700 3400 4350 4600]; 
% freqtonb = [664 800 1133 1550 2050 2380 2550 3300 4250 4500];
% 
% 2017-08-14_13-19-05
% freqtone = [530 800 1275 1500 2109 2500 2800 2950 3281 3516 3950 4300 4450 4700 5508 5900 6680 7200 7500 7900 8300 9200];
% freqtonb = [430 650 1150 1375 1992 2350 2700 2850 3164 3398 3850 4150 4300 4600 5391 5750 6560 7100 7400 7800 8200 9100];
%
%freqtone = (450:300:15000);
%freqtonb = (300:300:14900);
% %2017-08-14_13-22-04
% freqtone = [1211 1523 1758 2266 3000 4453 5600 6700 6900 8008]; 
% freqtonb = [1094 1328 1641 2148 2900 4336 5500 6600 6797 7852];
% 
% %2017-08-14_13-23-36
freqtone = [742 1719 3300 4050 4550 4922 5500 5989 7109 8711];
freqtonb = [586 1602 3200 3900 4450 4800 5400 5703 6992 8600];

% %2017-08-14_13-25-04
% freqtone = [1211 1523 1758 2266 3000 4453 5600 6700 6900 8008]; 
% freqtonb = [1094 1328 1641 2148 2900 4336 5500 6600 6797 7852];


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

%_________________________________________________
%PLOTTING
f = plot(time(1,:),K);

if Plot_Tonality == true
    ylim([0 0.5]);
    xlim([8 15]);
end   
 
if Save_Tonality == true
    saveas(f,'Tonality_plot.png')
end

end

%Loudness metric:
function phonresult = loundness(aspl)

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
phonvalues = 0:90;
phontable = zeros(90,29);
for i = 1:900
    k = i/10;
    phon = k;
    if((phon < 0) | (phon > 90))
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
phonresult = zeros(29,1);
 for i = 5:29
     for j = 1:900     
         if aspl(i) <= phontable(j,i)
             phonresult(i,1) = j/10;
             break
         end
        
     end
 end
 
end

%Sharpness metric:
function S = sharpness(phonresult)

%defining g(z) function
syms g(z)
g(z)= piecewise(z<=16, 1, z>16, 0.066*exp(0.171*z));

%defining the z that goes into the g(z) function
z_values_tab=zeros(length(phonresult),1);
g_values_tab=zeros(length(phonresult),1);
N_prime_tab=zeros(length(phonresult),1);
N_tab=zeros(length(phonresult),1);
top_tab=zeros(length(phonresult),1);
c=0.11;

for i=1:length(phonresult)
    z_value=13*atan(0.76*phonresult(i,1)/1000)+3.5*atan(phonresult(i,1)/7500).^2;
    z_values_tab(i,1)=z_value;
    
    %calculating a g(z), and N_prime function for every z_value
    g_value=vpa(g(z_value));
    g_values_tab(i,1)=g_value;
    
    %define Loudness function HERE:
    N_prime= phonresult(i,2);
    N_prime_tab(i,1)=N_prime;
    
    %defining N (bottom of sharpness function): integrate N'(z) between 0 and 24 
    N = trapz(z_value,N_prime_tab);
    
    %defining the top of sharpness function 
    top= g_values_tab(i)*N_prime_tab(i)*z_values_tab(i);
    top_tab(i,1)=top;
    a = trapz(z_value,top_tab);
    
    %sharpness function
    S= vpa((a/N)*c);
end


%plot(z_values_tab,top_tab)
%plot(z_values_tab,N_prime_tab)

end



