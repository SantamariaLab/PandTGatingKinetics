%% Ion Channel Thermodynamics (OLD) %%
% Jake Miller %
% 3/27/2023 %

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

%----------Pressure Plots------------------------------------------------%
%GENERATES OLD FIGURE

P = (0:100)*8*101325;          % Pressure from 1 to 800 atmospheres

% This explores entropy at different experimental temperatures

T_exp = transpose([10,20,52.33]+273.15);      % Experimental temperature
DSo_lo = -(R*T_exp(1,1).*log(kb*T_exp(1,1)/h)-DCp*(T_exp(1,1)-To-T_exp(1,1).*log(T_exp(1,1)/To))-DHo)./T_exp(1,1);
DSo_md = -(R*T_exp(2,1).*log(kb*T_exp(2,1)/h)-DCp*(T_exp(2,1)-To-T_exp(2,1).*log(T_exp(2,1)/To))-DHo)./T_exp(2,1);
DSo_hi = -(R*T_exp(3,1).*log(kb*T_exp(3,1)/h)-DCp*(T_exp(3,1)-To-T_exp(3,1).*log(T_exp(3,1)/To))-DHo)./T_exp(3,1);
DSo = DSo_md;

% Activation energy at three temperatures for three ions
DG_lo_Na = DG(P, T_exp(1,1), Po, To, DCp(1,1), DSo(1,1), DHo(1,1), DVo, Da(1,1), Dk(1,1));
DG_lo_K = DG(P, T_exp(1,1), Po, To, DCp(1,2), DSo(1,2), DHo(1,2), DVo, Da(1,1), Dk(1,1));
DG_lo_Ca = DG(P, T_exp(1,1), Po, To, DCp(1,3), DSo(1,3), DHo(1,3), DVo, Da(1,1), Dk(1,1));
DG_md_Na = DG(P, T_exp(2,1), Po, To, DCp(1,1), DSo(1,1), DHo(1,1), DVo, Da(1,1), Dk(1,1));
DG_md_K = DG(P, T_exp(2,1), Po, To, DCp(1,2), DSo(1,2), DHo(1,2), DVo, Da(1,1), Dk(1,1));
DG_md_Ca = DG(P, T_exp(2,1), Po, To, DCp(1,3), DSo(1,3), DHo(1,3), DVo, Da(1,1), Dk(1,1));
DG_hi_Na = DG(P, T_exp(3,1), Po, To, DCp(1,1), DSo(1,1), DHo(1,1), DVo, Da(1,1), Dk(1,1));
DG_hi_K = DG(P, T_exp(3,1), Po, To, DCp(1,2), DSo(1,2), DHo(1,2), DVo, Da(1,1), Dk(1,1));
DG_hi_Ca = DG(P, T_exp(3,1), Po, To, DCp(1,3), DSo(1,3), DHo(1,3), DVo, Da(1,1), Dk(1,1));

% Rate coefficient at three temperatures for three ions
k_lo_Na = k(T_exp(1,1), DG_lo_Na, kb, h, R);
k_lo_K = k(T_exp(1,1), DG_lo_K, kb, h, R);
k_lo_Ca = k(T_exp(1,1), DG_lo_Ca, kb, h, R);
k_md_Na = k(T_exp(2,1), DG_md_Na, kb, h, R);
k_md_K = k(T_exp(2,1), DG_md_K, kb, h, R);
k_md_Ca = k(T_exp(2,1), DG_md_Ca, kb, h, R);
k_hi_Na = k(T_exp(3,1), DG_hi_Na, kb, h, R);
k_hi_K = k(T_exp(3,1), DG_hi_K, kb, h, R);
k_hi_Ca = k(T_exp(3,1), DG_hi_Ca, kb, h, R);

figure('Position',[100 100 1200 300])
subplot(1,3,1)
hold on
plot(P/101325, k_lo_Na, 'r','LineWidth',1)
plot(P/101325, k_lo_K, 'g','LineWidth',1)
plot(P/101325, k_lo_Ca, 'b','LineWidth',1)
xlabel('Pressure (atm)')
ylabel('Rate Coefficient')
ylim([0 0.2])
title(['10 ' char(176) 'C'])
legend('Na','K','Ca','Location','southwest')
legend boxoff
box off
hold off

subplot(1,3,2)
hold on
plot(P/101325, k_md_Na, 'r','LineWidth',1)
plot(P/101325, k_md_K, 'g','LineWidth',1)
plot(P/101325, k_md_Ca, 'b','LineWidth',1)
xlabel('Pressure (atm)')
ylabel('Rate Coefficient')
ylim([0 1])
title(['20 ' char(176) 'C'])
legend('Na','K','Ca','Location','southwest')
legend boxoff
box off
hold off

subplot(1,3,3)
hold on
plot(P/101325, k_hi_Na, 'r','LineWidth',1)
plot(P/101325, k_hi_K, 'g','LineWidth',1)
plot(P/101325, k_hi_Ca, 'b','LineWidth',1)
xlabel('Pressure (atm)')
ylabel('Rate Coefficient')
ylim([0 15])
title(['52 ' char(176) 'C'])
legend('Na','K','Ca','Location','southwest')
legend boxoff
box off
hold off

%----------Adiabatic Heating---------------------------------------------%
%GENERATES OLD FIGURE

T_pts = transpose([20,52.33,84.67]+273.15);    % Temperatures to plot for
T_a = P/101325/200+T_pts;                      % Adding adiabtic heating
DSo_a = -(R*T_pts(1,1).*log(kb*T_pts(1,1)/h)-DCp*(T_pts(1,1)-To-T_pts(1,1).*log(T_pts(1,1)/To))-DHo)./T_pts(1,1);

% Activation energy at three temperatures for Sodium with NO heating
DG_lo_Na_n = DG(P, T_pts(1,1), Po, To, DCp(1,1), DSo_a(1,1), DHo(1,1), DVo, Da(1,1), Dk(1,1));
DG_md_Na_n = DG(P, T_pts(2,1), Po, To, DCp(1,1), DSo_a(1,1), DHo(1,1), DVo, Da(1,1), Dk(1,1));
DG_hi_Na_n = DG(P, T_pts(3,1), Po, To, DCp(1,1), DSo_a(1,1), DHo(1,1), DVo, Da(1,1), Dk(1,1));

% Rate coefficient at three temperatures for Sodium with NO heating
k_lo_Na_n = k(T_pts(1,1), DG_lo_Na_n, kb, h, R);
k_md_Na_n = k(T_pts(2,1), DG_md_Na_n, kb, h, R);
k_hi_Na_n = k(T_pts(3,1), DG_hi_Na_n, kb, h, R);

% Activation energy at three temperatures for Sodium with heating
DG_lo_Na_a = DG(P, T_a(1,:), Po, To, DCp(1,1), DSo_a(1,1), DHo(1,1), DVo, Da(1,1), Dk(1,1));
DG_md_Na_a = DG(P, T_a(2,:), Po, To, DCp(1,1), DSo_a(1,1), DHo(1,1), DVo, Da(1,1), Dk(1,1));
DG_hi_Na_a = DG(P, T_a(3,:), Po, To, DCp(1,1), DSo_a(1,1), DHo(1,1), DVo, Da(1,1), Dk(1,1));

% Rate coefficient at three temperatures for Sodium with heating
k_lo_Na_a = k(T_a(1,:), DG_lo_Na_a, kb, h, R);
k_md_Na_a = k(T_a(2,:), DG_md_Na_a, kb, h, R);
k_hi_Na_a = k(T_a(3,:), DG_hi_Na_a, kb, h, R);

figure('Position',[100 100 1200 300])
subplot(1,3,1)
hold on
plot(P/101325, k_lo_Na_n, 'r','LineWidth',1)
plot(P/101325, k_lo_Na_a,'--','LineWidth',1)
xlabel('Pressure (atm)')
ylabel('Rate Coefficient')
ylim([0 1])
title(['20 ' char(176) 'C'])
legend('Isothermal','Adiabatic Heating','Location','southwest')
legend boxoff
box off
hold off

subplot(1,3,2)
hold on
plot(P/101325, k_md_Na_n, 'r','LineWidth',1)
plot(P/101325, k_md_Na_a,'--','LineWidth',1)
xlabel('Pressure (atm)')
ylabel('Rate Coefficient')
ylim([0 15])
title(['52 ' char(176) 'C'])
legend('Isothermal','Adiabatic Heating','Location','southwest')
legend boxoff
box off
hold off

subplot(1,3,3)
hold on
plot(P/101325, k_hi_Na_n, 'r','LineWidth',1)
plot(P/101325, k_hi_Na_a,'--','LineWidth',1)
xlabel('Pressure (atm)')
ylabel('Rate Coefficient')
ylim([0 5])
title(['85 ' char(176) 'C'])
legend('Isothermal','Adiabatic Heating','Location','southwest')
legend boxoff
box off
hold off

%----------Optimal Temperature-------------------------------------------%

Topt_MMRT_Na = Topt(P, Po, To, DCp(1,1), DHo(1,1), DVo, Da(1,2), Dk(1,2), R) - 273.15;
Topt_MMRT_K = Topt(P, Po, To, DCp(1,2), DHo(1,2), DVo, Da(1,2), Dk(1,2), R) - 273.15;
Topt_MMRT_Ca = Topt(P, Po, To, DCp(1,3), DHo(1,3), DVo, Da(1,2), Dk(1,2), R) - 273.15;

Topt_MMRT_Na_a = Topt(P, Po, To, DCp(1,1), DHo(1,1), DVo, 1e-6, 0, R) - 273.15;
Topt_MMRT_K_a = Topt(P, Po, To, DCp(1,2), DHo(1,2), DVo, 1e-6, 0, R) - 273.15;
Topt_MMRT_Ca_a = Topt(P, Po, To, DCp(1,3), DHo(1,3), DVo, 1e-6, 0, R) - 273.15;

Topt_MMRT_Na_k = Topt(P, Po, To, DCp(1,1), DHo(1,1), DVo, 0, 1e-11, R) - 273.15;
Topt_MMRT_K_k = Topt(P, Po, To, DCp(1,2), DHo(1,2), DVo, 0, 1e-11, R) - 273.15;
Topt_MMRT_Ca_k = Topt(P, Po, To, DCp(1,3), DHo(1,3), DVo, 0, 1e-11, R) - 273.15;

figure('Position',[100 100 1200 300])

subplot(1,3,1)
hold on
plot(P/101325, Topt_MMRT_Na, 'r','LineWidth',1)
plot(P/101325, Topt_MMRT_K,'g','LineWidth',1)
plot(P/101325, Topt_MMRT_Ca,'b','LineWidth',1)
xlabel('Pressure (atm)')
ylabel(['Optimal Temperature (' char(176) 'C)'])
ylim([40 60])
legend('Na','K','Ca','Location','southwest')
legend boxoff
box off
hold off

subplot(1,3,2)
hold on
plot(P/101325, Topt_MMRT_Na_a, 'r','LineWidth',1)
plot(P/101325, Topt_MMRT_K_a,'g','LineWidth',1)
plot(P/101325, Topt_MMRT_Ca_a,'b','LineWidth',1)
xlabel('Pressure (atm)')
ylabel(['Optimal Temperature (' char(176) 'C)'])
ylim([40 60])
title('\alpha')
legend('Na','K','Ca','Location','southwest')
legend boxoff
box off
hold off

subplot(1,3,3)
hold on
plot(P/101325, Topt_MMRT_Na_k, 'r','LineWidth',1)
plot(P/101325, Topt_MMRT_K_k,'g','LineWidth',1)
plot(P/101325, Topt_MMRT_Ca_k,'b','LineWidth',1)
xlabel('Pressure (atm)')
ylabel(['Optimal Temperature (' char(176) 'C)'])
ylim([40 60])
title('\kappa')
legend('Na','K','Ca','Location','southwest')
legend boxoff
box off
hold off

%----------Potassium Pressure Plot---------------------------------------%

K_data = [1.36 2.02 3.30;  % 5 C      21 MPa  41 MPa  62 MPa
          1.28 1.83 2.67;  % 10 C
          1.27 1.69 2.54;  % 15 C
          1.31 1.91 2.83]; % Mean
P_vals = [21, 41, 62];
T_vals = [5, 10, 15];

figure('Position',[100 100 300 400])
hold on
plot(P_vals, log(K_data(1,:)), 'xk')
plot(P_vals, log(K_data(2,:)), '+k')
plot(P_vals, log(K_data(3,:)), 'ok','MarkerSize',3)
plot(P_vals, log(K_data(4,:)), '.k')
xlim([0 65])
xlabel('Pressure (MPa)')
ylim([0 1.25])
legend('5','10','15','Mean','Location','northwest')
legend boxoff
box off
hold off

%----------Fitting the Equation------------------------------------------%

data = xlsread('Pressure Data.xlsx',2); % Rate Pressure Temperature
rate_vals = data(:,1);
P_vals = data(:,2);
T_vals = data(:,3);

% Convert to SI units
P_vals_SI = P_vals*1e6; 
T_vals_SI = T_vals+273.15;

figure(5)
hold on
scatter3(P_vals(1:9), T_vals(1:9), rate_vals(1:9), 'filled','red')
scatter3(P_vals(10:24), T_vals(10:24), rate_vals(10:24), 'filled','green')
scatter3(P_vals(25:41), T_vals(25:41), rate_vals(25:41), 'filled','blue')
hold off

P_fit = [0:80]*1e6;                         
T_fit = transpose([0:80]+273.15);           

% Fit by hand %
DCp_fit = -2500;                 % Change in heat capacity (J / mol K)
DSo_fit = 60;                    % Change in reference entropy (J / mol K)
DHo_fit = 35000;                 % Change in reference enthalpy (J / mol)
DVo_fit = 40e-6;                 % Change in reference volume (m^3 / mol)
To_fit = 283.15;                 % Reference temperature (K)
Po_fit = 101325;                 % Reference pressure (Pa)
Da_fit = -1e-7;                  % Change in expansivity (m^3 / mol K)
Dk_fit = 1e-13;                  % Change in compressibility (m^3 / mol Pa)

DSo_fit = DHo_fit/To_fit+R*log(h/(kb*To_fit))
DG_fit = DG(P_fit, T_fit, Po_fit, To_fit, DCp_fit, DSo_fit, DHo_fit, DVo_fit, Da_fit, Dk_fit);
k_fit = k(T_fit, DG_fit, kb, h, R);

DG_hat = DG(P_vals_SI, T_vals_SI, Po_fit, To_fit, DCp_fit, DSo_fit, DHo_fit, DVo_fit, Da_fit, Dk_fit);
k_hat = k(T_vals_SI, DG_hat, kb, h, R);

SSE = sum((k_hat-rate_vals).^2)/mean(rate_vals)

hold on
surf(P_fit/1e6, T_fit-273.15, k_fit)
shading interp
alpha(0.5)
xlabel('Pressure (MPa)')
ylabel('Temperature (C)')
zlabel('Rate Coefficient')
hold off

% Set T=To=10 C and fit across pressure %

fit_pressure = pressure_fit(P_vals_SI(1:24), rate_vals(1:24));
coeffvalues(fit_pressure) % DHo DSo DVo Dk To

% Use pressure fit values and fit across T at 62 MPa

P_62_indices = P_vals == 62;
rate_vals_62 = rate_vals(P_62_indices);
T_vals_62 = T_vals(P_62_indices)+273.15;


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

%----------Optimal Temperature-------------------------------------------%

function Topt = Topt(P, Po, To, DCp, DHo, DVo, Da, Dk, R)
% Calculate the optimal tempertature as a function of pressure
    Topt = (DCp*To-DHo+Da*To*(P-Po)-DVo*(P-Po)+Dk/2*(P-Po).^2)/(DCp+R);
end

%------------------------------------------------------------------------%