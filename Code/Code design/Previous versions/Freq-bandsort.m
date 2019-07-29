%sup Wero! :)


lowfreq = [44.7,56.2,70.8,89.1,112,141,178,224,282,355,446,562,708,891,1122,1413,1778,2239,2818,3548,4467,5623,7079,8913,11220];
powert = power(:,index);
frequencyt = frequency(:,index); 
p0=20*10^(-6);

totindmax = 0;

for i = 1:length(lowfreq)-1
    x = 0;
    for j = 1:length(frequencyt)
        if (lowfreq(i) <= frequencyt(j,1) && frequencyt(j,1) < lowfreq(i+1))
            x=x+1;
        end
    end
    
    if x > totindmax
          totindmax = x;
    end 
end

indlst=zeros(length(lowfreq)-1,totindmax);


for k = 1:length(lowfreq)-1
    y = 0;
    for l = 1:length(frequencyt)
        if (lowfreq(k) <= frequencyt(l,1) && frequencyt(l,1) < lowfreq(k+1)) 
            y=y+1;
            indlst(k,y)= powert(l,1);
        end
        
    end
end

finalav=zeros(1,length(lowfreq));
q=0;
for m = 1:length(lowfreq)-1
    total=0;
    sum=0;
    for n = 1:totindmax

        if indlst(m,n)~=0
           total=total+1;
           sum=sum+indlst(m,n);
        end
    end
    av = sum/total
    spwav = 20*log(av/p0);
    q=q+1;
    finalav(1,q)= spwav;

end
