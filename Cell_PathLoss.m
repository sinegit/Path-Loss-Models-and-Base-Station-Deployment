%%
%Path Loss and Outage Probability for Cat-NB2
clc
Pt_NBIoT = -16;% Transmitted power in dBW
P_min_NBIoT = -135;%Minimum power in dBW
fc_NBIoT = 0.7; %As per the Xylem project document
d_2D_min = 10;
d_2D_max = 5000;
PL_RMa_LOS = zeros(1, d_2D_max);
PL_RMa_NLOS = zeros(1, d_2D_max);
PL_UMi_LOS = zeros(1, d_2D_max);
PL_UMi_NLOS = zeros(1, d_2D_max);
Outage_Prob_RMa_LOS = zeros(1, d_2D_max);
Outage_Prob_RMa_NLOS = zeros(1, d_2D_max);
Outage_Prob_UMi_LOS = zeros(1, d_2D_max);
Outage_Prob_UMi_NLOS = zeros(1, d_2D_max);
sigma_RMa_LOS_NBIoT = 4; 
sigma_RMa_NLOS_NBIoT = 8; 
sigma_UMi_LOS_NBIoT = 4; 
sigma_UMi_NLOS_NBIoT = 7.82; 
h = 5;%In m
h_BS_RMa = 35;
h_BS_UMi = 25;
h_UT = 1.5;
h_E=1;
W=20;
c = 3*10^8;%Speed of light
d_BP = (2*pi*h_BS_RMa*h_UT*fc_NBIoT*(10^9))/c;
d_prime_BP = (4*(h_BS_UMi-h_E)*(h_UT-h_E)*fc_NBIoT*(10^9))/c;
for d_2D = d_2D_min:d_2D_max
    d_3D = sqrt((d_2D^2)+((h_BS_RMa-h_UT)^2));
    if ((d_2D >= d_2D_min) && (d_2D <= d_BP))
        term1 = 20*log10(40*pi*d_3D*fc_NBIoT/3);
        term2 = min (0.03*h^(1.72),10)*log10(d_3D);
        term3 = min (0.044*h^(1.72),14.77);
        term4 = 0.002*log10(h)*d_3D;
        PL_1 = term1 + term2 -term3 +term4;
        PL_RMa_LOS(d_2D) = PL_1;
    elseif  ((d_2D >= d_BP) && (d_2D <= d_2D_max))
        term1 = 20*log10(40*pi*d_BP*fc_NBIoT/3);
        term2 = min (0.03*h^(1.72),10)*log10(d_BP);
        term3 = min (0.044*h^(1.72),14.77);
        term4 = 0.002*log10(h)*d_BP;
        PL_1 = term1 + term2 -term3 +term4;
        PL_2 = PL_1 + 40*log10(d_3D/d_BP);
        PL_RMa_LOS(d_2D) = PL_2;
    else
        disp('ERROR: Inappropriate range')
    end
    term5 = (24.37-3.7*(h/h_BS_RMa)^2)*log10(h_BS_RMa);
    term6 = (43.42-3.1*log10(h_BS_RMa))*(log10(d_3D)-3);
    term7 = 20*log10(fc_NBIoT) - (3.2*(log10(11.75*h_UT))^2 - 4.97);
    PL_prime_RMa_NLOS = 161.04-7.1*log10(W)+7.5*log10(h)-term5+term6+term7;
    PL_RMa_NLOS(d_2D) = max(PL_RMa_LOS(d_2D),PL_prime_RMa_NLOS);
    Outage_Prob_RMa_LOS(d_2D) = 1 - qfunc((P_min_NBIoT - (Pt_NBIoT - PL_RMa_LOS(d_2D)))/sigma_RMa_LOS_NBIoT);
    Outage_Prob_RMa_NLOS(d_2D) = 1 - qfunc((P_min_NBIoT - (Pt_NBIoT - PL_RMa_NLOS(d_2D)))/sigma_RMa_NLOS_NBIoT);
    
    d_3D_UMi = sqrt((d_2D^2)+((h_BS_UMi-h_UT)^2));
    PL_1_UMi = 32.4 + 21*log10(d_3D_UMi) + 20*log10(fc_NBIoT);
    PL_2_UMi = 32.4 + 40*log10(d_3D_UMi) + 20*log10(fc_NBIoT)-9.5*log10(((d_prime_BP)^2) + ((h_BS_UMi-h_UT)^2));
    if ((d_2D >= d_2D_min) && (d_2D <= d_prime_BP))
        PL_UMi_LOS(d_2D) = PL_1_UMi;
    elseif  ((d_2D >= d_prime_BP) && (d_2D <= d_2D_max))
        PL_UMi_LOS(d_2D) = PL_2_UMi;  
    else
        disp('ERROR: Inappropriate range')
    end
    PL_prime_UMi_NLOS = 35.3*log10(d_3D_UMi) + 22.4+21.3*log10(fc_NBIoT) - 0.3*(h_UT-1.5);
    PL_UMi_NLOS(d_2D) = max(PL_UMi_LOS(d_2D),PL_prime_UMi_NLOS);
    Outage_Prob_UMi_LOS(d_2D) = 1 - qfunc((P_min_NBIoT - (Pt_NBIoT - PL_UMi_LOS(d_2D)))/sigma_UMi_LOS_NBIoT);
    Outage_Prob_UMi_NLOS(d_2D) = 1 - qfunc((P_min_NBIoT - (Pt_NBIoT - PL_UMi_NLOS(d_2D)))/sigma_UMi_NLOS_NBIoT);
end
d_2D_axis = d_2D_min:d_2D_max;
figure()
plot(d_2D_axis,PL_RMa_LOS(d_2D_axis),d_2D_axis,PL_RMa_NLOS(d_2D_axis),d_2D_axis,PL_UMi_LOS(d_2D_axis),d_2D_axis,PL_UMi_NLOS(d_2D_axis),'LineWidth',3)
xlim([d_2D_min d_2D_max])
xlabel('Distance [m]','FontSize',20);
ylabel('Path Loss [dB]','FontSize',20);
legend('Rural: LOS', 'Rural: NLOS','Urban: LOS','Urban: NLOS')
set(gca,'FontSize',14) 
set(gca, 'XScale', 'log', 'YScale', 'log')
grid on
figure()
plot(d_2D_axis,Outage_Prob_RMa_LOS(d_2D_axis),d_2D_axis,Outage_Prob_RMa_NLOS(d_2D_axis),d_2D_axis,Outage_Prob_UMi_LOS(d_2D_axis),d_2D_axis,Outage_Prob_UMi_NLOS(d_2D_axis),'LineWidth',3)
xlim([d_2D_min d_2D_max])
xlabel('Distance [m]','FontSize',20);
title('Outage Probability for Cat NB2','FontSize',20);
legend('Rural: LOS', 'Rural: NLOS','Urban: LOS','Urban: NLOS')
set(gca,'FontSize',14) 
set(gca, 'XScale', 'log', 'YScale', 'log')
grid on
%%
%Outage Probability for Cat-M2
clc
Pt_LTE = -10;% Transmitted power in dBW
P_min_LTE = -125; % Transmitted power in dBW
fc_NBIoT = 0.7; %As per the Xylem project document
d_2D_min = 10;
d_2D_max = 5000;
PL_RMa_LOS = zeros(1, d_2D_max);
PL_RMa_NLOS = zeros(1, d_2D_max);
PL_UMi_LOS = zeros(1, d_2D_max);
PL_UMi_NLOS = zeros(1, d_2D_max);
Outage_Prob_RMa_LOS = zeros(1, d_2D_max);
Outage_Prob_RMa_NLOS = zeros(1, d_2D_max);
Outage_Prob_UMi_LOS = zeros(1, d_2D_max);
Outage_Prob_UMi_NLOS = zeros(1, d_2D_max);
sigma_RMa_LOS_NBIoT = 4; 
sigma_RMa_NLOS_NBIoT = 8; 
sigma_UMi_LOS_NBIoT = 4; 
sigma_UMi_NLOS_NBIoT = 7.82; 
h = 5;%In m
h_BS_RMa = 35;
h_BS_UMi = 25;
h_UT = 1.5;
h_E=1;
W=20;
c = 3*10^8;%Speed of light
d_BP = (2*pi*h_BS_RMa*h_UT*fc_NBIoT*(10^9))/c;
d_prime_BP = (4*(h_BS_UMi-h_E)*(h_UT-h_E)*fc_NBIoT*(10^9))/c;
for d_2D = d_2D_min:d_2D_max
    d_3D = sqrt((d_2D^2)+((h_BS_RMa-h_UT)^2));
    if ((d_2D >= d_2D_min) && (d_2D <= d_BP))
        term1 = 20*log10(40*pi*d_3D*fc_NBIoT/3);
        term2 = min (0.03*h^(1.72),10)*log10(d_3D);
        term3 = min (0.044*h^(1.72),14.77);
        term4 = 0.002*log10(h)*d_3D;
        PL_1 = term1 + term2 -term3 +term4;
        PL_RMa_LOS(d_2D) = PL_1;
    elseif  ((d_2D >= d_BP) && (d_2D <= d_2D_max))
        term1 = 20*log10(40*pi*d_BP*fc_NBIoT/3);
        term2 = min (0.03*h^(1.72),10)*log10(d_BP);
        term3 = min (0.044*h^(1.72),14.77);
        term4 = 0.002*log10(h)*d_BP;
        PL_1 = term1 + term2 -term3 +term4;
        PL_2 = PL_1 + 40*log10(d_3D/d_BP);
        PL_RMa_LOS(d_2D) = PL_2;
    else
        disp('ERROR: Inappropriate range')
    end
    term5 = (24.37-3.7*(h/h_BS_RMa)^2)*log10(h_BS_RMa);
    term6 = (43.42-3.1*log10(h_BS_RMa))*(log10(d_3D)-3);
    term7 = 20*log10(fc_NBIoT) - (3.2*(log10(11.75*h_UT))^2 - 4.97);
    PL_prime_RMa_NLOS = 161.04-7.1*log10(W)+7.5*log10(h)-term5+term6+term7;
    PL_RMa_NLOS(d_2D) = max(PL_RMa_LOS(d_2D),PL_prime_RMa_NLOS);
    Outage_Prob_RMa_LOS(d_2D) = 1 - qfunc((P_min_LTE - (Pt_LTE - PL_RMa_LOS(d_2D)))/sigma_RMa_LOS_NBIoT);
    Outage_Prob_RMa_NLOS(d_2D) = 1 - qfunc((P_min_LTE - (Pt_LTE - PL_RMa_NLOS(d_2D)))/sigma_RMa_NLOS_NBIoT);
    
    d_3D_UMi = sqrt((d_2D^2)+((h_BS_UMi-h_UT)^2));
    PL_1_UMi = 32.4 + 21*log10(d_3D_UMi) + 20*log10(fc_NBIoT);
    PL_2_UMi = 32.4 + 40*log10(d_3D_UMi) + 20*log10(fc_NBIoT)-9.5*log10(((d_prime_BP)^2) + ((h_BS_UMi-h_UT)^2));
    if ((d_2D >= d_2D_min) && (d_2D <= d_prime_BP))
        PL_UMi_LOS(d_2D) = PL_1_UMi;
    elseif  ((d_2D >= d_prime_BP) && (d_2D <= d_2D_max))
        PL_UMi_LOS(d_2D) = PL_2_UMi;  
    else
        disp('ERROR: Inappropriate range')
    end
    PL_prime_UMi_NLOS = 35.3*log10(d_3D_UMi) + 22.4+21.3*log10(fc_NBIoT) - 0.3*(h_UT-1.5);
    PL_UMi_NLOS(d_2D) = max(PL_UMi_LOS(d_2D),PL_prime_UMi_NLOS);
    Outage_Prob_UMi_LOS(d_2D) = 1 - qfunc((P_min_LTE - (Pt_LTE - PL_UMi_LOS(d_2D)))/sigma_UMi_LOS_NBIoT);
    Outage_Prob_UMi_NLOS(d_2D) = 1 - qfunc((P_min_LTE - (Pt_LTE - PL_UMi_NLOS(d_2D)))/sigma_UMi_NLOS_NBIoT);
end
d_2D_axis = d_2D_min:d_2D_max;
figure()
plot(d_2D_axis,Outage_Prob_RMa_LOS(d_2D_axis),d_2D_axis,Outage_Prob_RMa_NLOS(d_2D_axis),d_2D_axis,Outage_Prob_UMi_LOS(d_2D_axis),d_2D_axis,Outage_Prob_UMi_NLOS(d_2D_axis),'LineWidth',3)
xlim([d_2D_min d_2D_max])
xlabel('Distance [m]','FontSize',20);
title('Outage Probability for Cat M2','FontSize',20);
legend('Rural: LOS', 'Rural: NLOS','Urban: LOS','Urban: NLOS')
set(gca,'FontSize',14) 
set(gca, 'XScale', 'log', 'YScale', 'log')
grid on
