%defining g(z) function
syms g(z)
g(z)= piecewise(z<=16, 1, z>16, 0.066*exp(0.171*z));

%defining the z that goes into the g(z) function
z_values_tab=zeros(length(miguel),1);
g_values_tab=zeros(length(miguel),1);
N_prime_tab=zeros(length(miguel),1);
N_tab=zeros(length(miguel),1);
top_tab=zeros(length(miguel),1);
c=0.11;

for i=1:length(miguel)
    z_value=13*atan(0.76*miguel(i,1)/1000)+3.5*atan(miguel(i,1)/7500).^2;
    z_values_tab(i,1)=z_value;
    
    %calculating a g(z), and N_prime function for every z_value
    g_value=vpa(g(z_value));
    g_values_tab(i,1)=g_value;
    
    %define Loudness function HERE:
    N_prime= miguel(i,2);
    N_prime_tab(i,1)=N_prime;
    
    %defining N (bottom of sharpness function): integrate N'(z) between 0 and 24 
    N = trapz(z_value,N_prime_tab);
    
    %defining the top of sharpness function 
    top= g_values_tab(i)*N_prime_tab(i)*z_values_tab(i);
    top_tab(i,1)=top;
    a = trapz(z_value,top_tab);
    
    %sharpness function
    S= vpa((a/N)*c);
end

%plot(z_values_tab,top_tab)
%plot(z_values_tab,N_prime_tab)