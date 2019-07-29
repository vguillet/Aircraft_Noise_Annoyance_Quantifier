%_____________________________Inputs:______________________________________
signal = signal.';

c = 343;
h0 = 63.4;             %input initial altitude matching aircraft data set
v0 = 75.12;            %input velocity matching aircraft data set 
theta = 3;             %input angle of decent matching aircraft data set          
dt = 1/40000;
r0 = 1;

%__________________________________________________________________________
%Main code:

f1 = @BNC;
signal = f1(signal);

f2 = @UFC;
[decibels,T,F,psd1] = f2(signal);

f3 = @DSC;
[X,index,xq,yq,vq,dT] = f3(decibels,T,F,psd1,c,h0,v0);

f4 = @rdet;
r = f4(index,v0,theta,dT);

%f3 = @GSC;
%r = f3(signal, r, r0);

%---------------- PLOTTING --------------------%
mesh(xq,yq,vq) %Creates surface.
view(2) %Makes the view 2D.
colormap jet %Selects color.

%____________________Background noise correction:__________________________
function signal = BNC(signal)
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

%________________Uncorrected Frequencies conversion________________________
function [decibels,T,F,psd1] = UFC(signal)
[s,f,t] = spectrogram(signal,800,[],[],40000);
T=t.';
F=f;
psd1 = ((t(2)-t(1))^2/t(end))*(abs(s).^2);
decibels = 10*log10(psd1/((2*10^(-5))^2));
end

%______________________Doppler shift correction:___________________________
function [X,index,xq,yq,vq,dT] = DSC(decibels,T,F,psd1,c,h0,v0)

t = 0;
%-------------- CORRECTED FREQUENCIES ------------%
V=zeros(size(decibels,2),1);
dT=T(2)-T(1);
FNew=F;
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
[b,index] = max(sum(decibels,1));
index=index-h0/343*100;

%----------------------------------------------%

[xq,yq] = meshgrid(0:.05:T(end), 0:5:20000); %Creates X and Y axes for the grid.
vq = griddata(X(:,1),X(:,2),X(:,3),xq,yq); %Relates Z to X and Y axes.
end

%____________________Geometrical spreading correction:_____________________
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
    x(1,i) = [abs(xi)]; %Appending all the x positions in list x
end

%_____________________________y calculation:___________________________

ymax=vy*timeatx0;          %finding max height at start of descent

for i=1:length(timelist)
   
    yi=ymax-(vy*(i*dT)*sin(theta));
    y(1,i)= [yi];      %Appending all the y positions in list y
end

%____________________r magnitude calculation:______________________

%Computing r length using pythagoras
for i=1:length(timelist)
    
   ri=sqrt(x(i)*x(i)+y(i)*y(i));
   r(1,i) = [ri];
end

plot(timelist,r)

end
function signal = GSC(signal, r, r0)


end





