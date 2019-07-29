%rename "decibels_polished" to "decibels" in order to just use "decibels" in the
%following code.

decibels = decibels_polished;

%de-dopplerize
h0 = 63.4;
v0 = 75.12;
t = 0;
V=zeros(size(decibels,2),1);
dt=T(2)-T(1);
Fmat=zeros(size(decibels,1),size(decibels,2));
ALPHA=zeros(size(decibels,2),1);
[b,index] = max(sum(psd1,1));
index=index-h0/343*100;
disp(index);
idxMAX=round(index);%find on data


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

decibelscorr = decibels;

for i=1:1:size(decibelscorr,2)
    for j=1:1:size(decibelscorr,1)
        
    end    
end    

[xq,yq] = meshgrid(0:0.01:T(end), 0:5:20000);
vq = griddata(X(:,1),X(:,2),X(:,3),xq,yq);

mesh(xq,yq,vq)
% hold on
% plot3(x,y,v,'o')
% xlim([0 25])
% ylim([0 20000])
view(2)
colormap jet