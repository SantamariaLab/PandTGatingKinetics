%% Ion Channel Thermodynamics (OLD) %%
% Jake Miller %
% 3/8/2023 %

%----------Constants-----------------------------------------------------%

kb = 1.380649e-23;             % Boltzmann constant (m^2 kg / s^2 K)
h = 6.62607015e-34;            % Planck's constant (m^2 kg / s)
R= 8.3145;                     % Ideal gas constant (J / mol K)

%----------Parameters----------------------------------------------------%

DCp = [-2860,-4144,-3607];     % Change in heat capacity (J / mol K)
DSo = [63.1,59.5,87.4];        % Change in reference entropy (J / mol K)
DHo = [90110,88980,97180];     % Change in reference enthalpy (J / mol)
DVo = 19e-6;                   % Change in reference volume (m^3 / mol)
To = 298.15;                   % Reference temperature (K)
Po = 101325;                   % Reference pressure (Pa)
Da = [0,1e-8,-1e-8];           % Change in expansivity (m^3 / mol K)
Dk = [0,1e-13,-1e-13];         % Change in compressibility (m^3 / mol Pa)

%----------Isobaric------------------------------------------------------%
%GENERATES OLD FIGURE

% This explores entropy at different experimental temperatures

T_exp = transpose([0,20,40]+273.15);      % Experimental temperature
DSo_0 = -(R*T_exp(1,1).*log(kb*T_exp(1,1)/h)-DCp*(T_exp(1,1)-To-T_exp(1,1).*log(T_exp(1,1)/To))-DHo)./T_exp(1,1);
DSo_20 = -(R*T_exp(2,1).*log(kb*T_exp(2,1)/h)-DCp*(T_exp(2,1)-To-T_exp(2,1).*log(T_exp(2,1)/To))-DHo)./T_exp(2,1);
DSo_40 = -(R*T_exp(3,1).*log(kb*T_exp(3,1)/h)-DCp*(T_exp(3,1)-To-T_exp(3,1).*log(T_exp(3,1)/To))-DHo)./T_exp(3,1);
%DSo = DSo_20;

% This should show reducibility to the results of MMRT

T = transpose((0:100)+273.15);  % Temperature from 0 to 80 degrees Celsius
DG_T_Na = DG(Po, T, Po, To, DCp(1,1), DSo(1,1), DHo(1,1), DVo, Da(1,1), Dk(1,1));
DG_T_K = DG(Po, T, Po, To, DCp(1,2), DSo(1,2), DHo(1,2), DVo, Da(1,1), Dk(1,1));
DG_T_Ca = DG(Po, T, Po, To, DCp(1,3), DSo(1,3), DHo(1,3), DVo, Da(1,1), Dk(1,1));
k_T_Na = k(T, DG_T_Na, kb, h, R);
k_T_K = k(T, DG_T_K, kb, h, R);
k_T_Ca = k(T, DG_T_Ca, kb, h, R);

figure('Position',[100 100 400 300])
hold on
plot(T-273.15, k_T_Na,'r','LineWidth',1)
plot(T-273.15, k_T_K,'g','LineWidth',1)
plot(T-273.15, k_T_Ca,'b','LineWidth',1)
xlabel(['Temperature (' char(176) 'C)'])
ylabel('Rate Coefficient')
ylim([0 15])
title('Macromolecular Rate Theory')
legend('Na','K','Ca')
legend boxoff
box off
hold off

% Calculate average optimal T
[~, i_Na] = max(k_T_Na);
[~, i_K] = max(k_T_K);
[~, i_Ca] = max(k_T_Ca);
av_T_opt = (T(i_Na)+T(i_K)+T(i_Ca))/3-273.15

%----------Isothermal----------------------------------------------------%
%GENERATES OLD FIGURE

% This section should show reducibility to pressure results

P = (0:100)*8*101325;          % Pressure from 1 to 800 atmospheres
DG_P_Na = DG(P, To, Po, To, DCp(1,1), DSo(1,1), DHo(1,1), DVo, Da(1,1), Dk(1,1));
DG_P_K = DG(P, To, Po, To, DCp(1,2), DSo(1,2), DHo(1,2), DVo, Da(1,1), Dk(1,1));
DG_P_Ca = DG(P, To, Po, To, DCp(1,3), DSo(1,3), DHo(1,3), DVo, Da(1,1), Dk(1,1));
k_P_Na = k(To, DG_P_Na, kb, h, R);
k_P_K = k(To, DG_P_K, kb, h, R);
k_P_Ca = k(To, DG_P_Ca, kb, h, R);

% Double the volume parameter
DG_P_Na_2 = DG(P, To, Po, To, DCp(1,1), DSo(1,1), DHo(1,1), DVo*3, Da(1,1), Dk(1,1));
DG_P_K_2 = DG(P, To, Po, To, DCp(1,2), DSo(1,2), DHo(1,2), DVo*3, Da(1,1), Dk(1,1));
DG_P_Ca_2 = DG(P, To, Po, To, DCp(1,3), DSo(1,3), DHo(1,3), DVo*3, Da(1,1), Dk(1,1));
k_P_Na_2 = k(To, DG_P_Na_2, kb, h, R);
k_P_K_2 = k(To, DG_P_K_2, kb, h, R);
k_P_Ca_2 = k(To, DG_P_Ca_2, kb, h, R);

figure('Position',[100 100 800 300])
subplot(1,2,1)
hold on
plot(P/101325, k_P_Na, 'r','LineWidth',1)
plot(P/101325, k_P_K, 'g','LineWidth',1)
plot(P/101325, k_P_Ca, 'b','LineWidth',1)
xlabel('Pressure (atm)')
ylabel('Rate Coefficient')
ylim([0 2.25])
title('High Pressure Effects')
legend('Na','K','Ca')
legend boxoff
box off
hold off

subplot(1,2,2)
hold on
plot(P/101325, k_P_Na_2, 'r','LineWidth',1)
plot(P/101325, k_P_K_2, 'g','LineWidth',1)
plot(P/101325, k_P_Ca_2, 'b','LineWidth',1)
xlabel('Pressure (atm)')
ylabel('Rate Coefficient')
ylim([0 2.25])
title('High Pressure Effects - Triple Activation Volume')
legend('Na','K','Ca')
legend boxoff
box off
hold off

% If we consider adiabatic heating which is not exactly isothermal
%GENERATES OLD FIGURE

T_ah = P/101325/200+To;        % Adiabatic heating of 1 degree per 200 atm
DG_ah = DG(P, T_ah, Po, To, DCp(1,1), DSo(1,1), DHo(1,1), DVo, Da(1,1), Dk(1,1));
k_ah = k(T_ah, DG_ah, kb, h, R);

DG_ah_2 = DG(P, T_ah, Po, To, DCp(1,1), DSo(1,1), DHo(1,1), DVo*3, Da(1,1), Dk(1,1));
k_ah_2 = k(T_ah, DG_ah_2, kb, h, R);

figure('Position',[100 100 800 300])
subplot(1,2,1)
hold on
plot(P/101325, k_P_Na, 'r','LineWidth',1)
plot(P/101325, k_ah,'--','LineWidth',1)
xlabel('Pressure (atm)')
ylabel('Rate Coefficient')
ylim([0 2.25])
title('Sodium Pressure Effects')
legend('Isothermal','Adiabatic Heating')
legend boxoff
box off
hold off

subplot(1,2,2)
hold on
plot(P/101325, k_P_Na_2, 'r','LineWidth',1)
plot(P/101325, k_ah_2,'--','LineWidth',1)
xlabel('Pressure (atm)')
ylabel('Rate Coefficient')
ylim([0 2.25])
title('Sodium Pressure Effects - Triple Activation Volume')
legend('Isothermal','Adiabatic Heating')
legend boxoff
box off
hold off

%----------Constant Activation Volume------------------------------------%
%GENERATES OLD FIGURE

DG_PT_Na = DG(P, T, Po, To, DCp(1,1), DSo(1,1), DHo(1,1), DVo, Da(1,1), Dk(1,1));
DG_PT_K = DG(P, T, Po, To, DCp(1,2), DSo(1,2), DHo(1,2), DVo, Da(1,1), Dk(1,1));
DG_PT_Ca = DG(P, T, Po, To, DCp(1,3), DSo(1,3), DHo(1,3), DVo, Da(1,1), Dk(1,1));
k_PT_Na = k(T, DG_PT_Na, kb, h, R);
k_PT_K = k(T, DG_PT_K, kb, h, R);
k_PT_Ca = k(T, DG_PT_Ca, kb, h, R);

figure('Position',[100 100 1200 375])
subplot(1,3,1)
surf(P/101325, T-273.15, k_PT_Na)
shading interp
xlabel('Pressure (atm)')
xlim([0 800])
ylabel(['Temperature (' char(176) 'C)'])
zlabel('Rate Coefficient')
zlim([0 15])
caxis([0 15])
title('Na')
subplot(1,3,2)
surf(P/101325, T-273.15, k_PT_K)
shading interp
xlabel('Pressure (atm)')
xlim([0 800])
ylabel(['Temperature (' char(176) 'C)'])
zlabel('Rate Coefficient')
zlim([0 15])
caxis([0 15])
title('K')
subplot(1,3,3)
surf(P/101325, T-273.15, k_PT_Ca)
shading interp
xlabel('Pressure (atm)')
xlim([0 800])
ylabel(['Temperature (' char(176) 'C)'])
zlabel('Rate Coefficient')
zlim([0 15])
caxis([0 15])
title('Ca')

%----------Pressure-Temperature Dependence-------------------------------%
%GENERATES OLD FIGURE

DG_PT_0 = DG(P, T, Po, To, DCp(1,1), DSo(1,1), DHo(1,1), DVo, Da(1,1), Dk(1,1));
DG_PT_kp = DG(P, T, Po, To, DCp(1,1), DSo(1,1), DHo(1,1), DVo, Da(1,1), Dk(1,2));
DG_PT_kn = DG(P, T, Po, To, DCp(1,1), DSo(1,1), DHo(1,1), DVo, Da(1,1), Dk(1,3));
DG_PT_ap = DG(P, T, Po, To, DCp(1,1), DSo(1,1), DHo(1,1), DVo, Da(1,2), Dk(1,1));
DG_PT_an = DG(P, T, Po, To, DCp(1,1), DSo(1,1), DHo(1,1), DVo, Da(1,3), Dk(1,1));
DG_PT_pos = DG(P, T, Po, To, DCp(1,1), DSo(1,1), DHo(1,1), DVo, Da(1,2), Dk(1,2));
DG_PT_neg = DG(P, T, Po, To, DCp(1,1), DSo(1,1), DHo(1,1), DVo, Da(1,3), Dk(1,3));
DG_PT_pn = DG(P, T, Po, To, DCp(1,1), DSo(1,1), DHo(1,1), DVo, Da(1,2), Dk(1,3));
DG_PT_np = DG(P, T, Po, To, DCp(1,1), DSo(1,1), DHo(1,1), DVo, Da(1,3), Dk(1,2));

k_PT_0 = k(T, DG_PT_0, kb, h, R);
k_PT_kp = k(T, DG_PT_kp, kb, h, R);
k_PT_kn = k(T, DG_PT_kn, kb, h, R);
k_PT_ap = k(T, DG_PT_ap, kb, h, R);
k_PT_an = k(T, DG_PT_an, kb, h, R);
k_PT_pos = k(T, DG_PT_pos, kb, h, R);
k_PT_neg = k(T, DG_PT_neg, kb, h, R);
k_PT_pn = k(T, DG_PT_pn, kb, h, R);
k_PT_np = k(T, DG_PT_np, kb, h, R);

figure('Position',[100 100 800 750]);
hold on

subplot(3,3,1)
surf(P/101325, T-273.15, k_PT_0)
shading interp
grid off
xlabel('Pressure (atm)')
ylabel(['Temperature (' char(176) 'C)'])
zlabel('Rate Coefficient')
view(0,90)
xlim([0 800])
ylim([0 100])
zlim([0 15])
caxis([0 15])
title('\alpha 0, \kappa 0')

subplot(3,3,2)
surf(P/101325, T-273.15, k_PT_kp)
shading interp
grid off
xlabel('Pressure (atm)')
ylabel(['Temperature (' char(176) 'C)'])
zlabel('Rate Coefficient')
view(0,90)
xlim([0 800])
ylim([0 100])
zlim([0 15])
caxis([0 15])
title('\alpha 0, \kappa +')

subplot(3,3,3)
surf(P/101325, T-273.15, k_PT_kn)
shading interp
grid off
xlabel('Pressure (atm)')
ylabel(['Temperature (' char(176) 'C)'])
zlabel('Rate Coefficient')
view(0,90)
xlim([0 800])
ylim([0 100])
zlim([0 15])
caxis([0 15])
title('\alpha 0, \kappa -')

subplot(3,3,4)
surf(P/101325, T-273.15, k_PT_ap)
shading interp
grid off
xlabel('Pressure (atm)')
ylabel(['Temperature (' char(176) 'C)'])
zlabel('Rate Coefficient')
view(0,90)
xlim([0 800])
ylim([0 100])
zlim([0 15])
caxis([0 15])
title('\alpha +, \kappa 0')

subplot(3,3,7)
surf(P/101325, T-273.15, k_PT_an)
shading interp
grid off
xlabel('Pressure (atm)')
ylabel(['Temperature (' char(176) 'C)'])
zlabel('Rate Coefficient')
view(0,90)
xlim([0 800])
ylim([0 100])
zlim([0 15])
caxis([0 15])
title('\alpha -, \kappa 0')

subplot(3,3,5)
surf(P/101325, T-273.15, k_PT_pos)
shading interp
grid off
xlabel('Pressure (atm)')
ylabel(['Temperature (' char(176) 'C)'])
zlabel('Rate Coefficient')
view(0,90)
xlim([0 800])
ylim([0 100])
zlim([0 15])
caxis([0 15])
title('\alpha +, \kappa +')

subplot(3,3,9)
surf(P/101325, T-273.15, k_PT_neg)
shading interp
grid off
xlabel('Pressure (atm)')
ylabel(['Temperature (' char(176) 'C)'])
zlabel('Rate Coefficient')
view(0,90)
xlim([0 800])
ylim([0 100])
zlim([0 15])
caxis([0 15])
title('\alpha -, \kappa -')

subplot(3,3,6)
surf(P/101325, T-273.15, k_PT_pn)
shading interp
grid off
xlabel('Pressure (atm)')
ylabel(['Temperature (' char(176) 'C)'])
zlabel('Rate Coefficient')
view(0,90)
xlim([0 800])
ylim([0 100])
zlim([0 15])
caxis([0 15])
cb_5 = colorbar('Position', [0.93 0.075 0.02 0.85]);
cb_5.Label.String = 'Rate Coefficient';
title('\alpha +, \kappa -')

subplot(3,3,8)
surf(P/101325, T-273.15, k_PT_np)
shading interp
grid off
xlabel('Pressure (atm)')
ylabel(['Temperature (' char(176) 'C)'])
zlabel('Rate Coefficient')
view(0,90)
xlim([0 800])
ylim([0 100])
zlim([0 15])
caxis([0 15])
title('\alpha -, \kappa +')

hold off

% If we circumscribe to positive activation volumes
%GENERATES OLD FIGURE

DV_0 = DV(P, T, Po, To, DVo, Da(1,1), Dk(1,1));
DV_kp = DV(P, T, Po, To, DVo, Da(1,1), Dk(1,2));
DV_kn = DV(P, T, Po, To, DVo, Da(1,1), Dk(1,3));
DV_ap = DV(P, T, Po, To, DVo, Da(1,2), Dk(1,1));
DV_an = DV(P, T, Po, To, DVo, Da(1,3), Dk(1,1));
DV_pos = DV(P, T, Po, To, DVo, Da(1,2), Dk(1,2));
DV_neg = DV(P, T, Po, To, DVo, Da(1,3), Dk(1,3));
DV_pn = DV(P, T, Po, To, DVo, Da(1,2), Dk(1,3));
DV_np = DV(P, T, Po, To, DVo, Da(1,3), Dk(1,2));

DV_0_bool = (DV_0>=0);
DV_kp_bool = (DV_kp>=0);
DV_kn_bool = (DV_kn>=0);
DV_ap_bool = (DV_ap>=0);
DV_an_bool = (DV_an>=0);
DV_pos_bool = (DV_pos>=0);
DV_neg_bool = (DV_neg>=0);
DV_pn_bool = (DV_pn>=0);
DV_np_bool = (DV_np>=0);

k_pos_0 = k_PT_0.*DV_0_bool;
k_pos_kp = k_PT_kp.*DV_kp_bool;
k_pos_kn = k_PT_kn.*DV_kn_bool;
k_pos_ap = k_PT_ap.*DV_ap_bool;
k_pos_an = k_PT_an.*DV_an_bool;
k_pos_pos = k_PT_pos.*DV_pos_bool;
k_pos_neg = k_PT_neg.*DV_neg_bool;
k_pos_pn = k_PT_pn.*DV_pn_bool;
k_pos_np = k_PT_np.*DV_np_bool;

figure('Position',[100 100 800 750])
hold on

subplot(3,3,1)
surf(P/101325, T-273.15, k_pos_0)
shading interp
grid off
xlabel('Pressure (atm)')
ylabel(['Temperature (' char(176) 'C)'])
zlabel('Rate Coefficient')
view(0,90)
xlim([0 800])
ylim([0 100])
zlim([0 15])
caxis([0 15])
title('\alpha 0, \kappa 0')

subplot(3,3,2)
surf(P/101325, T-273.15, k_pos_kp)
shading interp
grid off
xlabel('Pressure (atm)')
ylabel(['Temperature (' char(176) 'C)'])
zlabel('Rate Coefficient')
view(0,90)
xlim([0 800])
ylim([0 100])
zlim([0 15])
caxis([0 15])
title('\alpha 0, \kappa +')

subplot(3,3,3)
surf(P/101325, T-273.15, k_pos_kn)
shading interp
grid off
xlabel('Pressure (atm)')
ylabel(['Temperature (' char(176) 'C)'])
zlabel('Rate Coefficient')
view(0,90)
xlim([0 800])
ylim([0 100])
zlim([0 15])
caxis([0 15])
title('\alpha 0, \kappa -')

subplot(3,3,4)
surf(P/101325, T-273.15, k_pos_ap)
shading interp
grid off
xlabel('Pressure (atm)')
ylabel(['Temperature (' char(176) 'C)'])
zlabel('Rate Coefficient')
view(0,90)
xlim([0 800])
ylim([0 100])
zlim([0 15])
caxis([0 15])
title('\alpha +, \kappa 0')

subplot(3,3,7)
surf(P/101325, T-273.15, k_pos_an)
shading interp
grid off
xlabel('Pressure (atm)')
ylabel(['Temperature (' char(176) 'C)'])
zlabel('Rate Coefficient')
view(0,90)
xlim([0 800])
ylim([0 100])
zlim([0 15])
caxis([0 15])
title('\alpha -, \kappa 0')

subplot(3,3,5)
surf(P/101325, T-273.15, k_pos_pos)
shading interp
grid off
xlabel('Pressure (atm)')
ylabel(['Temperature (' char(176) 'C)'])
zlabel('Rate Coefficient')
view(0,90)
xlim([0 800])
ylim([0 100])
zlim([0 15])
caxis([0 15])
title('\alpha +, \kappa +')

subplot(3,3,9)
surf(P/101325, T-273.15, k_pos_neg)
shading interp
grid off
xlabel('Pressure (atm)')
ylabel(['Temperature (' char(176) 'C)'])
zlabel('Rate Coefficient')
view(0,90)
xlim([0 800])
ylim([0 100])
zlim([0 15])
caxis([0 15])
title('\alpha -, \kappa -')

subplot(3,3,6)
surf(P/101325, T-273.15, k_pos_pn)
shading interp
grid off
xlabel('Pressure (atm)')
ylabel(['Temperature (' char(176) 'C)'])
zlabel('Rate Coefficient')
view(0,90)
xlim([0 800])
ylim([0 100])
zlim([0 15])
caxis([0 15])
cb_5 = colorbar('Position', [0.93 0.075 0.02 0.85]);
cb_5.Label.String = 'Rate Coefficient';
title('\alpha +, \kappa -')

subplot(3,3,8)
surf(P/101325, T-273.15, k_pos_np)
shading interp
grid off
xlabel('Pressure (atm)')
ylabel(['Temperature (' char(176) 'C)'])
zlabel('Rate Coefficient')
view(0,90)
xlim([0 800])
ylim([0 100])
zlim([0 15])
caxis([0 15])
title('\alpha -, \kappa +')

hold off

%----------Rate Coefficient----------------------------------------------%

function k = k(T, G, kb, h, R)
% Calculate the rate coefficient
    k = (kb*T/h).*exp(-G./(R*T));
end

%----------Energy Barrier------------------------------------------------%

function DG = DG(P, T, Po, To, DCp, DSo, DHo, DVo, Da, Dk)
% Calculate the change in Gibbs free energy across activation barrier
    DG = DCp*(T-To)-DCp*T.*log(T/To)-T*DSo+Da*(T-To).*(P-Po)+DVo*(P-Po)-Dk/2*(P-Po).^2+DHo;
end

%----------Activation Volume---------------------------------------------%

function DV = DV(P, T, Po, To, DVo, Da, Dk)
% Calculate the activation volume based on its P and T dependence
    DV = DVo + Da*(T-To)-Dk*(P-Po);
end

%------------------------------------------------------------------------%