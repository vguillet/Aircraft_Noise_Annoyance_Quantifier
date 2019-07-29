function epnl_value = epnlm(power, frequency, time)
%____________________Effective Perceived Noise Level_____________________

%   1.  After applying corrections, you'll have tables "time", "power" and
%       "frequency" in the workspace.
%   2.  Upload file "bands.xlsx".
%   3.  Run and wait ~30s.
%   4.  The output of this code is the value called "epnl". It's just one
%       number.

%________average sound pressure level________

avspltot = zeros(24,size(time,2));
for h = 1:size(time,2)
    
    pow = power(:,h);
    freq = frequency(:,h); 
    freqsorted = sort(freq);

    % make a matrix with frequency in the first column and power in the second
    % column, sorted so that frequency goes from smallest value to biggest,
    % keeping the appropriate power in the second column
    freqpowsorted = [freqsorted,zeros(length(freqsorted),1)];
    for n = 1:length(freqsorted)
        x = freqsorted(n,1);
        i = find(freq == x);
        freqpowsorted(n,2) = pow(i,1);
    end

    % delete rows where power = -inf
    % the matrix shortens after every row removal so the "r" is there to
    % account for it
    r = 0;
    for m = 1:length(freqpowsorted)
        if freqpowsorted((m-r),2) < 0 
            freqpowsorted((m-r),:)=[];
            r = r + 1;
        end
    end

    lowfreq = [44.7,56.2,70.8,89.1,112,141,178,224,282,355,446,562,708,891,1122,1413,1778,2239,2818,3548,4467,5623,7079,8913,11220];
    lowfreq = transpose(lowfreq);
    
    avspl = zeros(24,1);
    plist = zeros(length(freqpowsorted),(length(lowfreq)-1)); 
    % plist has 24 columns for every band
    % for every band, pressures will be calsulated for points which frequencies
    % lye in the band, where every p is a sum of the previous and the current p
    lenlist = zeros(24,1);
    % for every band, find the frequencies and calculate the average spl 
    for k = 1:(length(lowfreq)-1)

        p = 0;
        for t = 1:length(freqpowsorted)
            if freqpowsorted(t,1)>lowfreq(k,1) && freqpowsorted(t,1)<lowfreq((k+1),1)
                powerr = 10^(freqpowsorted(t,2)/10);
                p = p + sqrt(powerr);
                plist(t,k) = p;
            end
        end

        len = nnz(plist(:,k)); % number of non zero elements, this is the number of points that lye in this band     
        lenlist(k,1) = len;
        pav = max(plist(:,k))/len; % final pressure sum divided by number of points
        
        avspl(k,1) = 20*log10(pav); % average sound pressure level

    end
    
   
    avsplcol = avspl(:,1);
    avspltot(:,h) = avsplcol;
    
end

avspltot(isnan(avspltot)) = 0;
for row = 1:size(avspltot,2)
    for cell = 1:size(avspltot,1)
        if avspltot(cell,row) < 0
            avspltot(cell,row)=0;
        end
    end
end   

pnlt = zeros(size(time,2),1);
for b = 1:size(time,2)
    %________________________aspl to pnl_____________________________ 

    % aspl to noys 
    band = table2array(bands(:,2));
    L1 = table2array(bands(:,4));
    L2 = table2array(bands(:,6));
    L3 = table2array(bands(:,8));
    Lc = table2array(bands(:,10));
    L4 = table2array(bands(:,12));
    M1 = table2array(bands(:,5));
    M2 = table2array(bands(:,7));
    M3 = table2array(bands(:,9));
    M4 = table2array(bands(:,11));
  

    N_list = zeros((numel(band)-1),1);
    for n = 1:(numel(band)-1)    
        if L1(n,1)<avspltot(n,b)<L2(n,1)
            N = 0.1*(10^(M1(n,1)*(avspltot(n,b)-L1(n,1))));
            N_list(n,1) = N; 
        elseif L2(n,1)<avspltot(n,b)<L3(n,1) 
            N = 10^(M2(n,1)*(avspltot(n,b)-L3(n,1)));
            N_list(n,1) = N; 
        elseif L3(n,1)<avspltot(n,b)<Lc(n,1) 
            N = 10^(M3(n,1)*(avspltot(n,b)-L3(n,1)));
            N_list(1,n) = N; 
        elseif Lc(n,1)<avspltot(n,b)<150 
            N = 10^(M4(n,1)*(avspltot(n,b)-L4(n,1)));
            N_list(n,1) = N; 
        end
    end

    % summing the noys   
    Nmax = max(N_list);
    N_list(isnan(N_list)) = 0; 
    Nsum = sum(N_list);
    Ntot = 0.85*Nmax + 0.15*Nsum;

    % total number of noys to pnl 
    pnltot = 40 + 10/(log10(2))*(log10(Ntot));

    %______________________ tone corrected pnl _________________________

    % step 1 - computing arithmetic difference in sound pressure level
   
    dspl = zeros(numel(band),1);
    apr = zeros(numel(band),1);
    dpr = zeros(numel(band),1);
    
    for n = 4:(numel(band)-1)
        apr(n,1) = 10^(avspltot(n,b)/20);   % apr - average pressure
    end
    
    r = 0;
    nonzerolist = apr;
    for w = 1:length(nonzerolist)
        if nonzerolist(w-r) == 1 || nonzerolist(w-r) == 0
            nonzerolist(w-r) = [];
            r = r + 1;
        end
    end

    dpr = zeros(length(nonzerolist),1);
    for cell = 2:length(nonzerolist)
        dpr(cell,1) = abs(nonzerolist(cell) - nonzerolist(cell-1));
    end
    
    i = 0;
    for r = 1:length(apr)
        if apr(r,1) ~= 0 && apr(r,1) ~= 1
            i = i + 1;
            dspl(r,1) = 20*log10(dpr(i,1)); % difference in sound pressure levels
            if i == 1
                dspl(r,1) = 0;
            end
        end
    end
    
    
    % step 2 - values of dspl bigger than 5
    % step 3 - choosing aspls
    column = avspltot(:,b);
    aspl = [column,zeros(length(column),1)]; % adding a column of zeros to aspl list

    for n = 5:numel(band)
        if abs(dspl(n,1)-dspl((n-1),1))> 5 
            if dspl(n,1)>0 && dspl(n,1)>dspl((n-1),1)
                aspl(n,2) = 1;
            elseif dspl(n,1)<=0 && dspl((n-1),1)>0
                aspl((n-1),2) = 1;
            end
        end
    end

    % step 4 - new adjusted spl

    aspl = [aspl,zeros(length(aspl),1)];

    for n = 3:(numel(band)-2)
        if aspl(n,2) == 1 
            aspl(n,3) = 0.5*(aspl((n-1),1)) + aspl((n+1),1);
        else 
            aspl(n,3) = aspl(n,1);
        end
    end

    for n = (numel(band)-1)
        if aspl(24,2) == 1
            aspl(24,3) = aspl(23,1) + dspl(23,1);
        else 
            aspl(24,3) = aspl(24,1);
        end
    end

    % step 5 - add dspl(3) and dspl(25)
    
    dspl1 = zeros(numel(band),1);
    pres = zeros(numel(band),1);
    dpres = zeros(numel(band),1);
    for n = 4:(numel(band)-1)
        pres(n,1) = 10^((aspl(n,3))/20);  % pres - pressure
        dpres(n,1) = abs(pres(n,1)-pres((n-1),1));   % arithmetic difference
        dspl1(n,1) = 20*log10(dpres(n,1));  % difference in sound pressure levels
    end
    
    for n = 3
        dspl1(n,1) = dspl1((n+1),1);
    end

    for n = length(band)
        dspl1(n,1) = dspl1((n-1),1);
    end
    
    for cell = 1:length(dspl1)
        if dspl1(cell,1) < 0
            dspl1(cell,1) = 0;
        end
    end
    
    % step 6 - arithmetic average of 3 adjacent slopes

    dspl3 = zeros(numel(band),1);
    pnew = zeros(numel(band),1);
    p3 = zeros(numel(band),1);
    for n = 3:23
        pnew(n,1) = 10^(dspl1(n,1)/20);                     % pnew - new pressure
        if pnew(n,1) == 1 || pnew((n+1),1) == 1 || pnew((n+2),1) == 1
            p3(n,1) = 0;
        else
            p3(n,1) = (pnew(n,1) + pnew((n+1),1) + pnew((n+2),1))/3;    % arithmetic average
            dspl3(n,1) = 20*log10(p3(n,1));                     % logarithmic average
        end
    end

    % step 7 - final 1/3 octave band spls

    aspl = [aspl,zeros(length(aspl),1)];

    for n = 3
       aspl(n,4) = aspl(n,1);
    end

    for n = 4:23
        aspl(n,4) = aspl((n-1),4) + dspl3((n-1),1);
    end

    % step 8 - differences between aspl(1,n) and aspl(4,n)

    F = zeros(numel(band),1);
    for n = 3:24
       F(n,1) = aspl(n,1) - aspl(n,4);
    end

    % step 9 - tone correction factors

    C = zeros(numel(band),1);
    for n = 3:11
       if 1.5 <= F(n,1) < 3 
          C(n,1) = F(n,1)/3 - 0.5;
       elseif 3 <= F(n,1) < 20 
          C(n,1) = F(n,1)/6;
       elseif 20 <= F(n,1) 
          C(n,1) = 10/3;
       end
    end

    for n = 12:21
       if 1.5 <= F(n,1) < 3 
          C(n,1) = F(n,1)/3*2 - 1;
       elseif 3 <= F(n,1) < 20
          C(n,1) = F(n,1)/3;
       elseif F(n,1) >= 20 
          C(n,1) = 20/3;
       end
    end

    for n = 22:24
       if 1.5 <= F(n,1) < 3
          C(n,1) = F(n,1)/3 - 0.5;
       elseif 3 <= F(n,1) < 20
          C(n,1) = F(n,1)/6;
       elseif F(n,1) >= 20
          C(n,1) = 10/3;
       end
    end

    % step 10 - tone corrected pnl
    
    pnlt(b,1) = pnltot + max(C(:,1));
    
end

antilog_pnlt = zeros(length(pnlt),1);
for g = 1:length(pnlt)
    antilog_pnlt = 10^(pnlt(g,1)-2);
end
antilog_pnlt_sum = sum(antilog_pnlt);
t1 = 0;
t2 = time(1,size(time,2));
epnl_value = log10(antilog_pnlt_sum/length(pnlt))-13;
disp(epnl_value)

end




















