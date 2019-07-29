t = time(:,250);
p = power(:,250);
f = frequency(:,250);
ptab = zeros(end,end);
%power(:, all(isnan(power), 1) ) = [];
ftab = zeros(end,end);
ttab = zeros(end,end);
qtab = zeros(end,end);
q = mean(p)
y = median(p)
for i = 1:length(power(:,1))
    p = power(i,:);
    q = nanmean(p);
    qtab(i,1) = q;
end
%plot(q,time(1,:))
for i = 1:length(time(1,:))
    p = power(:,i);
    f = frequency(:,i);
    x = 1.2*mean(p);

    [pton,y] = findpeaks(p,'MinPeakDistance',225,'MinPeakHeight',100,'MinPeakProminence',25,'MinPeakWidth',0);
%      p1 = min(pton);
% tallPeakIndexes = pton > 140; % Whatever value you want.
% pton(tallPeakIndexes) = 0;    
 % Or you can use someThreshold instead of minValue.
    fton = frequency(y,i);
    tton = time(y,i);
    for j = 1:length(pton)
        ptab(j,i) = pton(j);
        ftab(j,i) = fton(j);
        ttab(j,i) = tton(j);
    end
end
for i = 1:length(ftab(1,:))
    if unique(A
end
plot = stem3(ttab,ftab,ptab,'LineStyle','none');
%plot(f,p);
%plot(ptab,ftab)
%mesh(time,frequency,power);
%view(2)    

%findpeaks(frequency,power));
% Fil_tab=zeros(1,end);
% a = 1;
% j = 1;
% z = 1;
% b = 1;
% c =1;
% d =1;
% tonalc_tab=zeros(end,end);
% dp_tab=zeros(1,end);
% Fre2_tab=zeros(1,end);
% 
% for i=1:length(X)
%     if X(i,3) > 40
%         Fil_tab(3,a) = X(i,3);
%         Fil_tab(2,a) = X(i,2);
%         Fil_tab(1,a) = X(i,1);
%         a = a + 1;
%     end
% end
% z_values_tab=zeros(1,length(Fre_tab));
% CBW_tab=zeros(1,length(Fre_tab));

% for i=1:length(Fre_tab)
%     z_value=13*atan(0.76*Fre_tab(i)/1000)+3.5*atan(Fre_tab(i)/7500).^2;
%     z_values_tab(1,i)=z_value;
%     CBW = 25+75*(1+1.4*(Fre_tab(i)/1000)^2)^0.69;
%     CBW_tab(1,i)=CBW;    
% end
% 
% for i=2:length(Fre_tab)
%  dF = Fre_tab(i)-Fre_tab(j);
%     if dF >= CBW_tab(j)
%      dp = Pow_tab(i)-Pow_tab(j);
%      dp_tab(1,z)=dp;
%      Fre2_tab(1,z)=Fre_tab(i);
%      if Fre_tab(z)-Fre_tab(b) >= CBW_tab(j)
%          if (dp_tab(1,b) > 7) && (dp_tab(1,z) < -7)
%                 tonalc_tab(1,c) = Pow_tab(i);
%                 tonalc_tab(2,c) = Fre_tab(i);
%                 tonalc_tab(3,d) = Fre_tab(j);
%                  tonalc_tab(4,d) = Pow_tab(i);
%                  tonalc_tab(5,d) = dp_tab(1,b);
%                  tonalc_tab(6,d) = dp_tab(1,z);
%                  d = d +1;
%                  tonalc_tab(5,d) = dp_tab(1,b);
%                  tonalc_tab(6,d) = dp_tab(1,z);
%                  tonalc_tab(3,d) = Fre_tab(i);
%                  tonalc_tab(4,d) = Pow_tab(i);
%                  d = d + 1;
% 
%                  c = c + 1;
%                 
%          end
%          b = b + 1;
%      end
%      z = z + 1;
%      j = j + 1;
%     end
% end
% 
% hold off;
% %plot(Fre2_tab,dp_tab)
% %hold on;
% semilogx(Fre_tab,Pow_tab);
% for i = 2:2:length(tonalc_tab)
%     %hold on;
%     j =i-1;
%     %plot(tonalc_tab(3,j:i),tonalc_tab(4,j:i))
%     %hold on;
%     %plot(tonalc_tab(3,j:i),tonalc_tab(5,j:i))
%     hold on;
%     %plot(tonalc_tab(3,j:i),tonalc_tab(6,j:i))
% end
% 
% % %Main Code
% % Ehs_tab = zeros(1,end);
% % for 1:length(Fre_tab)
% %     Ehs_tab(1,i) = 3.64*(Fre_tab(1,i)/1000)^(-0.8)
% % end

% syms s(f_i)
% s(f_i) = piecewise(tonalc_tab(2,i) <= tonalc_tab(2,k), 27, tonalc_tab(2,i) > tonalc_tab(2,k), -24-(230/tonalc_tab(2,k))+0.2*tonalc_tab(1,k));
% 
% L_Ek(f_i) = tonalc_tab(1,k) - s(z_k - z_i)
% 
% A_Ek(f_i) = 10*exp(L_Ek(f_i)/20)

