%_____________________________Inputs:______________________________________
signal = signal.';
%bnsignal = --

c = 343;
h0 = 63.4;             %input initial altitude matching aircraft data set
v0 = 75.12;            %input velocity matching aircraft data set 
theta = 3;             %input angle of decent matching aircraft data set          
dt = 1/40000;
r0 = 1;


%Main code:----------------------------------------------------------------
%--------------------------------------------------------------------------
%Corrections:

%Applying first background noise correction (frequency) to signal
f1 = @BNfc;
signal = f1(signal);

% test = @Tonc;
% signal = test(signal);
%Applying fast fourier transform to aircraft signal
f2 = @FFTs;
[decibels,T,F,psd1] = f2(signal);

%Applying fast fourier transform to background noise signal
f3 = @FFTbns;
[bndecibels,bnpsd1] = f3(bnsignal);

%Applying second background noise correction (decibels) to decibels using
%bndecibels
% f4 = @BNdc;
% decibels = f4(bndecibel);

%Applying doppler effect correction
f5 = @DSC;
[X,index,xq,yq,vq,dT] = f5(decibels,T,F,psd1,c,h0,v0);

%Calculating r(distance aircraft-mic) for every point
f6 = @rdet;
r = f6(index,v0,theta,dT);

%Applying geometrical spreading correction
f7 = @GSC;
vqcorr = f7(vq, r, r0);

%Renaming variables
time = xq;
frequency = yq;
power = vqcorr;

%PLOTTING
mesh(time,frequency,power)  %Creates surface.
view(2)                     %Makes the view 2D.
colormap jet                %Selects color.

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%Metrics:





%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%Core functions:


%Corrections:
%________________Background noise (frequency) correction:__________________
function signal = BNfc(signal)
%This bit of code corrects for background noise with functions
%like 'butter' and 'filtfilt' in matlab. 
%fc represents the 'cutoff frequency' and can be changed. 

fc = 300;
fs = 40000;
fcnorm = fc/(fs/2);
[b,a] = butter(4,fcnorm,'high');
% Obtaining corrected data 

signal  = filtfilt(b,a,signal); 

end

function signal = Tonc(signal)
freqtone = [650 900 1250 1950 2800 4250 5400 6000 6500 9100 11700]; 
freqtonb = [300 750 1000 1700 2700 4100 5300 5900 6300 8950 11500];
freqton = (freqtone+freqtonb)/2;
datan = zeros(1,1);
fN = 20000;
N_order = 4;
for i = 1:length(freqton)
    [b,a] = butter(N_order,[freqtonb(i) freqtone(i)]/fN);
    data = filtfilt(b,a,signal);
    data(isnan(data))=0;
    datan = datan + data;
end
 signal = signal - datan;

end
%________________Fast Fourier Transform (signal):__________________________
function [decibels,T,F,psd1] = FFTs(signal)
[s,f,t] = spectrogram(signal,800,[],[],40000);
T=t.';
F=f;
psd1 = ((t(2)-t(1))^2/t(end))*(abs(s).^2);
decibels = 10*log10(psd1/((2*10^(-5))^2));
end

%________________Fast Fourier Transform (bnsignal):________________________
function [bndecibels,bnpsd1] = FFTbns(bnsignal)
[s,~,t] = spectrogram(bnsignal,800,[],[],40000);
bnpsd1 = ((t(2)-t(1))^2/t(end))*(abs(s).^2);
bndecibels = 10*log10(bnpsd1/((2*10^(-5))^2));
end

%________________Background noise (decibels) correction:___________________

function decibels = BNdc(decibels,bndecibels)

sum = 0;
for idx = 1:numel(bndecibels)
    element = bndecibels(idx);
    sum = sum + element;
    average = sum/idx ;
end

%Now taking this average and setting the ones that are below
%this to 0. 

for idx = 1:numel(decibels)     %These should be the array of the aircraft decibels 
    if decibels(idx) < average
        decibels(idx) = 0; 
    end
end
end

%________________Doppler shift correction:_________________________________
function [X,index,xq,yq,vq,dT] = DSC(decibels,T,F,psd1,c,h0,v0)

%-------------- CORRECTED FREQUENCIES ------------%
V=zeros(size(decibels,2),1);
dT=T(2)-T(1);
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

%---------------- Max decibel locating for r determination-----------%
[~,index] = max(sum(decibels,1));
index=index-h0/343*100;

%----------------------------------------------%

[xq,yq] = meshgrid(0:.05:T(end), 0:5:20000); %Creates X and Y axes for the grid.
vq = griddata(X(:,1),X(:,2),X(:,3),xq,yq); %Relates Z to X and Y axes.
end

%________________Geometrical spreading correction:_________________________
function r = rdet(index,v0,theta,dT)

timelist = 0.01:dT:24.99;   

y=zeros(1,length(timelist));
x=zeros(1,length(timelist));
r=zeros(1,length(timelist));

%_____________________________v calculation:___________________________

vy=v0*sin(theta*(pi/180));         %speed positive down
vx=v0*cos(theta*(pi/180));

%_____________________________x calculation:_________________________

%find x=0 position
timeatx0 = index*0.01;

%finding most negative x distance from microphone
xneg = -(timeatx0*vx);

%finding x distance from microphone for every datapoint
for i=1:length(timelist)
    
    xi=xneg+((i*dT)*vx);
    x(1,i) = abs(xi); %Appending all the x positions in list x
end

%_____________________________y calculation:___________________________

ymax=vy*timeatx0;          %finding max height at start of descent

for i=1:length(timelist)
   
    yi=ymax-(vy*(i*dT)*sin(theta));
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


%Metrics:


