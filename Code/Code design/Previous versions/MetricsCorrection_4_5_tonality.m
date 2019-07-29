%Welcome to the sounds correction software by CD8&co.
%To use, follow the steps provided bellow:
%!!!in the same order!!!:
    %1)Import background noise data, rename it to "bnsignal"
    %2)Import pressure data to be tested (should be named signal automaticaly).
    %3)Run and be patient, the program takes about 20~30s to generate diagram

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
%Main code:
%==========================================================================

%-----------------------------------------------------
%Corrections:
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

%Calculating r(distance aircraft-mic) for every point
f6 = @rdet;
[r,timelist,x] = f6(index,v0,theta,xq,h0,dT);

%Applying geometrical spreading correction
f7 = @GSC;
vqcorr = f7(vq, r, r0);

%Renaming variables
time = xq;
frequency = yq;
power = vqcorr;

%PLOTTING
mesh(time(:,1000:1200),frequency(:,1000:1200),power(:,1000:1200))  %Creates surface.
view(2)                     %Makes the view 2D.
colormap jet                %Selects color.
colorbar





%-----------------------------------------------------
%Metrics:
%-----------------------------------------------------

% Tonality
tonal = @tonality;
K = tonal(power,frequency,time);

% Plot
figure = plot(time(1,:),K);
ylim([0 0.5]);
xlim([8 15]);
% 
% saveas(figure,'2017-08-14_13-16-48.png')

%==========================================================================
%Core functions:
%==========================================================================

%-----------------------------------------------------
%Corrections:
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

%plot(timelist,r)

end
function vqcorr = GSC(vq,r,r0)

vqcorr = zeros(size(vq));

for i = 1:size(vq,2)
    for j = 1:size(vq,1)
        vqcorr(j,i) = vq(j,i) + 20*log10(r(i)/r0);
    end
end

end

%-----------------------------------------------------
%Metrics:
%-----------------------------------------------------

function K = tonality(power,frequency,time)

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

end
