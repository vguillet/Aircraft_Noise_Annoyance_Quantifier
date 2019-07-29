%_____________________________r determination:_________________________

%Parameters to doublecheck/redefine:

v0 = 68.46;             %input velocity matching aircraft data set
h0 = 60.44;            %input initial altitude matching aircraft data set
theta = 5;             %input angle of decent matching aircraft data set
dt=1/40000;

testlist=signal;    %input pressure dataset


%______________________________________________________________________
%______________________________________________________________________

y=zeros(1,length(testlist));
x=zeros(1,length(testlist));
r=zeros(1,length(testlist));
t=zeros(1,length(testlist));

%_____________________________v calculation:___________________________

vy=v0*sin(theta*(pi/180));         %speed positive down
vx=v0*cos(theta*(pi/180));

%_____________________________t calculation:___________________________

t=length(testlist)*dt;

for i=1:length(testlist) %Creating t step list
    
    ti=dt*i;
    t(1,i)= ti;      %Appending all the t steps in list t
end

%_____________________________x calculation:_________________________

%find max position
[M,I] = max(testlist);

%finding max x distance from microphone
xmax = (I*dt)*vx;

%finding x distance from microphone for every datapoint
for i=1:length(testlist)
    
    xi=xmax-((i*dt)*vx);
    x(1,i) = abs(xi); %Appending all the x positions in list x
end

%_____________________________y calculation:___________________________

ymax=vy*(I*dt)+h0;          %finding max height at start of descent

for i=1:length(testlist)
   
    yi=ymax+(vy*(i*dt)*sin(theta));
    y(1,i)= [yi];      %Appending all the y positions in list y
end

%____________________r magnitude calculation:______________________

%Computing r length using pythagoras
for i=1:length(testlist)
    
   ri=sqrt(x(i)*x(i)+y(i)*y(i));
   r(1,i) = ri;
end

plot(t,y)