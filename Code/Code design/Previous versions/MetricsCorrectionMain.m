%_____________________________Inputs:______________________________________
signal = signal.';

%Parameters to doublecheck/redefine:

c = 343;
h0 = 63.48;             %input initial altitude matching aircraft data set
v0 = 75.12;             %input velocity matching aircraft data set 
theta=3                 %input angle of decent matching aircraft data set          

%_____________________________Doppler shift correction:____________________
t = 0;

[s,f,t] = spectrogram(signal,800,[],[],40000);
T=t.';
F=f;
psd1 = ((t(2)-t(1))^2/t(end))*(abs(s).^2);
decibels = 10*log10(psd1/((2*10^(-5))^2));

% figure,
% surf(decibels,'Edgecolor','none')
% % view(2)
% colormap jet
% colorbar

V=zeros(size(decibels,2),1);
dt=T(2)-T(1);
FNew=F;
Fmat=zeros(size(decibels,1),size(decibels,2));
ALPHA=zeros(size(decibels,2),1);
[b,index] = max(sum(decibels,1));
idxMAX= round(index-h0/c/0.01);%find on data
display (idxMAX)
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

% tri = delaunay(X(:,1),X(:,2));
% %plot(X(:,1),X(:,2),'.')
% [r,c] = size(tri);
% disp(r)
% h = trisurf(tri,X(:,1),X(:,2),X(:,3));
% %axis vis3d
% l = light('Position',[-50 -15 29]);
% %set(gca,'CameraPosition',[208 -50 0])
% shading interp
% colorbar EastOutside
% view(2)

%[xq,yq] = meshgrid(0:.1:25, 0:5:20000);
%vq = griddata(X(:,1),X(:,2),X(:,3),xq,yq);

%mesh(xq,yq,vq)
%hold on
%plot3(X(:,1),X(:,2),X(:,3),'o')
%xlim([0 25])
%ylim([0 20000])
%view(2)
%colormap jet

[b,index] = max(sum(decibels,1));
index=index-h0/343*100;

%_____________________________r determination:_________________________

timelist = 0.01:dt:24.99;   

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
    
    xi=xneg+((i*dt)*vx);
    x(1,i) = [abs(xi)]; %Appending all the x positions in list x
end

%_____________________________y calculation:___________________________

ymax=vy*timeatx0;          %finding max height at start of descent

for i=1:length(timelist)
   
    yi=ymax-(vy*(i*dt)*sin(theta));
    y(1,i)= [yi];      %Appending all the y positions in list y
end

%____________________r magnitude calculation:______________________

%Computing r length using pythagoras
for i=1:length(timelist)
    
   ri=sqrt(x(i)*x(i)+y(i)*y(i));
   r(1,i) = [ri];
end

plot(timelist,r)