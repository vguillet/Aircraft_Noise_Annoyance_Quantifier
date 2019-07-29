%_____________________________Inputs:______________________________________
signal = signal.';

c = 343;
h0 = 63.40;             %input initial altitude matching aircraft data set
v0 = 75.12;             %input velocity matching aircraft data set 
theta = 3;                %input angle of decent matching aircraft data set          
dt = 1/40000;
r0 = 1;
%__________________________________________________________________________
%Main code:

f1 = @BNC;
signal = f1(signal);

%f2 = @rdet;
%r = f2(signal,h0,v0,theta,dt);

%f3 = @GSC;
%r = f3(signal, r, r0);

f4 = @DSC;
X = f4(signal,c,h0,v0);

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

%____________________Geometrical spreading correction:_____________________
function r = rdet(signal,h0,v0,theta,dt)

testlist = signal;

y=zeros(1,length(testlist));
x=zeros(1,length(testlist));
r=zeros(1,length(testlist));
t=zeros(1,length(testlist));

%_____________________________v calculation:_______________________________

vy=v0*sin(theta*(pi/180));         %speed positive down
vx=v0*cos(theta*(pi/180));

%_____________________________t calculation:_______________________________

t = length(testlist)*dt;

for i=1:length(testlist) %Creating t step list
    
    ti=dt*i;
    t(1,i)= ti;      %Appending all the t steps in list t
end

%_____________________________x calculation:_______________________________

%find max position
[~,I] = max(testlist);

%finding max x distance from microphone
xmax = (I*dt)*vx;

%finding x distance from microphone for every datapoint
for i=1:length(testlist)
    
    xi=xmax-((i*dt)*vx);
    x(1,i) = abs(xi); %Appending all the x positions in list x
end

%_____________________________y calculation:_______________________________

ymax=vy*(I*dt)+h0;          %finding max height at start of descent

for i=1:length(testlist)
   
    yi=ymax-(vy*(i*dt)*sin(theta));
    y(1,i)= yi;      %Appending all the y positions in list y
end

%____________________r magnitude calculation:______________________

%Computing r length using pythagoras
for i=1:length(testlist)
    
   ri=sqrt(x(i)*x(i)+y(i)*y(i));
   r(1,i) = ri;
end

%plot(t,r)

end
function signal = GSC(signal, r, r0)

for i=1:length(signal)
    signali = signal(i) + 20*log(r(i)/r0);
    signal(1,i)= signali;
end
end
%______________________Doppler shift correction:___________________________
function X = DSC(signal,c,h0,v0)

t = 0;

%----------- Uncorrected Frequencies -----------%
[s,f,t] = spectrogram(signal,800,[],[],40000);
T=t.';
F=f;
psd1 = ((t(2)-t(1))^2/t(end))*(abs(s).^2);
decibels = 10*log10(psd1/((2*10^(-5))^2));

% figure,
% surf(t,f,decibels,'Edgecolor','none')
% % view(2)
% colormap jet
% colorbar

%-------------- CORRECTED FREQUENCIES ------------%
V=zeros(size(decibels,2),1);
dt=T(2)-T(1);
FNew=F;
Fmat=zeros(size(decibels,1),size(decibels,2));
ALPHA=zeros(size(decibels,2),1);
[b,index] = max(sum(psd1,1)); %Find maximum power at all times.
idxMAX= round(index-h0/c/0.01);%Approximation of passing time index.
%display (idxMAX)

%-------------- DE-DOPLERIZE --------------%

for j=1:1:size(decibels,1)
    t=0;
    i=1;
while t<=T(idxMAX)
    t = t+dt;
    DeltaT=T(idxMAX)-t;
    alpha=abs(atan(h0/(v0*DeltaT)));
    ALPHA(i)=alpha;
    dr = v0*dt*cos(alpha);
    V(i,1)=-dr/dt;
    Fmat(j,i)=1/(1-V(i)/343)*F(j);
    i=i+1;
end 
while t<T(end)
    t = t+dt;
    DeltaT=t-T(idxMAX);
    alpha =abs(atan(h0/(v0*DeltaT)));
    ALPHA(i)=alpha;
    dr = v0*dt*sin(alpha);
    V(i,1)=dr/dt;
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

%---------------- PLOTTING --------------------%

[xq,yq] = meshgrid(0:.05:T(end), 0:5:20000); %Creates X and Y axes for the grid.
vq = griddata(X(:,1),X(:,2),X(:,3),xq,yq); %Relates Z to X and Y axes.

mesh(xq,yq,vq) %Creates surface.
view(2) %Makes the view 2D.
colormap jet %Selects color.
end




