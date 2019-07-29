signal = signal.';

%----------- PARAMETERS ------------%

c = 343; %Change depending on temperature
h0 = 63.48; %Change for each signal
v0 = 75.12; %change for each signal
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

[xq,yq] = meshgrid(0:.1:T(end), 0:30:20000); %Creates X and Y axes for the grid.
vq = griddata(X(:,1),X(:,2),X(:,3),xq,yq); %Relates Z to X and Y axes.

mesh(xq,yq,vq) %Creates surface.
view(2) %Makes the view 2D.
colormap jet %Selects color.

