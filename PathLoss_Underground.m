clc
clear
d_min=0.05;
d_max = 0.55;
d_range = 0.05;
d_num = (d_max-d_min)/d_range;
PLa = zeros(1,d_num);
PLb = zeros(1,d_num);
eps_prime=13.25;
eps_2_prime = 2.18;
eps0 = 8.854*10^(-12);
f=700;
sigma=0.1;
frac_term = (eps_2_prime + (sigma/2*pi*f*eps0))/eps_prime;
deno_term = sqrt((eps_prime/2)*(1+sqrt(1+frac_term^2)));
alpha = (8.68*60*pi*((2*pi*f*eps0*eps_2_prime)+sigma))/deno_term;
R = ((1-sqrt(eps_prime))/(1+sqrt(eps_prime)))^2;
Rc = 10*log10((2*R)/(1+R));
mu_s=1; %ralative Permeability for soil
mu = mu_s*4*pi*10^(-7);
term1 = sqrt(1+(eps_2_prime/eps_prime)^2);
term2 = (mu*eps_prime)/2;
alpha2 = 2*pi*f*sqrt(term2*(term1-1));
beta = 2*pi*f*sqrt(term2*(term1+1));
r=0;
for d = d_min:d_range:d_max
    r=r+1;
    PLa(r) = alpha*d + Rc; 
    PLb(r) = 6.45+20*log10(d)+20*log10(beta)+8.69*alpha2*d;
end

dist = d_min:d_range:d_max;
figure()
plot(dist,PLa,dist,PLb,'LineWidth',3)
xlim([d_min d_max])
xlabel('Distance [m]','FontSize',20);
ylabel('Path Loss [dBm]','FontSize',20);
legend('CRIM-Fresnel', 'Modified-Friis')
set(gca,'FontSize',14) 

