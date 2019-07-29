freqtone = [650 1600 4250 5200 5600 6700 3000 8050 9450]; 
freqtonb = [550 1200 4100 5100 5500 6500 2500 7950 9300];
freqton = (freqtone+freqtonb)/2;
freqtonttab = zeros(end,end);
numfreqton = freqton /5+1;
freq = [1:5:20000];
powton = [];
data = zeros(end,end);

fN = 20000;
N_order = 4;
% for i = 1:length(freqton)
%     [b,a] = butter(N_order,[freqtonb(6) freqtone(6)]/fN);
%     data_f = filtfilt(b,a,power);
%     data_f(isnan(data_f))=0;
%     data = data + data_f;
% end
% data = data*1.1;
% %mesh(time,frequency,power);
% %hold on;
% power2 = power+data;
% power3 = power-data;
% mesh(time,frequency,power3);
%mesh(time(100:3000,51:321),frequency(100:3000,51:321),power3(100:3000,51:321));
title('power2');
view(2)

powtont = [];
powtonttot = [];
numfreqtont = zeros(end,end);
for i = 1:length(freqton)
powton = [powton; power(numfreqton(i),:)];
freqtont = [freqtonb(i):5:freqtone(i)];
powtont = [];
    for j = 1:length(freqtont)
    freqtonttab(i,j) = freqtont(j);
    numfreqtont(i,j)= freqtont(j) /5+1;
    powtont = [powtont; power(numfreqtont(i,j),:)];
    end
powtonttot = [powtonttot; sum(powtont)];
end
powtonttottot = sum(powtonttot,'omitnan');
powtot = sum(power,'omitnan');
powres = powtot-powtonttottot;



totlen = length(freqton)*length(powton(1,:));
L_Ek = zeros(end,end);
A_Ek = zeros(end,end);
stab = zeros(end,end);
dLirl = zeros(end,end);
freqtonrl = zeros(end,end);
Atot = zeros(end,end);
stottab = [];
L_Ektab = [];
A_Ektab = [];
Ehs_tab = 3.64*(freqton/1000).^(-0.8)-6.5*exp(-0.6*(freqton/1000-3.3).^2)+10^(-3)*(freqton/1000).^4;
for i = 1:length(Ehs_tab)
    if Ehs_tab(i) < 0
        Ehs_tab(i) = 0;
    end
end
z= 13*atan(0.76*freqton/1000)+3.5*atan(freqton/7500).^2;
zb = 13*atan(0.76*freqtonb/1000)+3.5*atan(freqtonb/7500).^2;
ze = 13*atan(0.76*freqtone/1000)+3.5*atan(freqtone/7500).^2;
dz = ze-zb;
Ehs2_tab = 3.64*(freq/1000).^(-0.8)-6.5*exp(-0.6*(freq/1000-3.3).^2)+10^(-3)*(freq/1000).^4;
%plot(Ehs2_tab,freq);

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
   dLi = powton(i,a)-10*log(Atot(i,a)^2+Ehs_tab(i)+10);
   if dLi > 0
       dLirl(a,i) = dLi;
       freqtonrl(a,i) = freqton(i);
   end
 end
   L_Ektab = [L_Ektab;L_Ek];
   A_Ektab = [A_Ektab;A_Ek];
end

w3 = (1-exp(-(dLirl/15))).^0.29;
w2 = ((1+0.2*(freqtonrl*1/700+700./freqtonrl).^2).^(0.5)).^(-.29);
w1 = 0.13./(0.13+dz);
w = (w1.^(1/0.29).*w2.^(1/0.29).*w3.^(1/0.29)).^2;
wT = sum(w,2).^(0.5);
wGr = transpose(1-powres./powtot);
K = 1.09*wT.^0.29;
plot(time(1,:),K);
ylim([0.1 0.32]);