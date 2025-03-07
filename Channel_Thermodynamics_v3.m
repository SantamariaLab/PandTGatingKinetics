%% Ion Channel Thermodynamics %%
% Jake Miller %
% 4/21/2023 %

%----------Constants-----------------------------------------------------%

kb = 1.380649e-23;             % Boltzmann constant (m^2 kg / s^2 K)
h = 6.62607015e-34;            % Planck's constant (m^2 kg / s)
R= 8.3145;                     % Ideal gas constant (J / mol K)

%----------Parameters----------------------------------------------------%

DCp = -2860;                   % Change in heat capacity (J / mol K)
DSo = 63.1;                    % Change in reference entropy (J / mol K)
DHo = 90110;                   % Change in reference enthalpy (J / mol)
DVo = 19e-6;                   % Change in reference volume (m^3 / mol)
To = 298.15;                   % Reference temperature (K)
Po = 101325;                   % Reference pressure (Pa)
Da = [1e-8;1e-6;-1e-6];        % Change in expansivity (m^3 / mol K)
Dk = [1e-13;1e-11;-1e-11];     % Change in compressibility (m^3 / mol Pa)

%----------Temperature and Pressure--------------------------------------%

T = transpose((0:100)+273.15); % Temperature from 0 to 100 degrees Celsius
P = (0:100)*10*101325;         % Pressure from 1 to 800 atmospheres

%----------Figure 1------------------------------------------------------%

% Atmosphere, Average Ocean, Mariana Trench Pressures
P_vals = [1;369;1072]*101325;  

%-----Isobars MMRT
DG_T_1 = DG(P_vals(1), T, Po, To, DCp, DSo, DHo, DVo, Da(3), Dk(1));
DG_T_100 = DG(P_vals(2), T, Po, To, DCp, DSo, DHo, DVo, Da(3), Dk(1));
DG_T_1000 = DG(P_vals(3), T, Po, To, DCp, DSo, DHo, DVo, Da(3), Dk(1));
k_T_1 = k(T, DG_T_1, kb, h, R);
k_T_100 = k(T, DG_T_100, kb, h, R);
k_T_1000 = k(T, DG_T_1000, kb, h, R);

%-----Optimal Temperatures on MMRT 
T_opt_1 = Topt(P_vals(1), Po, To, DCp, DHo, DVo, Da(3), Dk(1), R);
T_opt_100 = Topt(P_vals(2), Po, To, DCp, DHo, DVo, Da(3), Dk(1), R);
T_opt_1000 = Topt(P_vals(3), Po, To, DCp, DHo, DVo, Da(3), Dk(1), R);

%-----Optimal Temperatures on Pressure
T_opt_0 = Topt(P, Po, To, DCp, DHo, DVo, Da(1), Dk(1), R);
T_opt_a = Topt(P, Po, To, DCp, DHo, DVo, Da(3), Dk(1), R);
T_opt_k = Topt(P, Po, To, DCp, DHo, DVo, Da(1), Dk(2), R);
T_opt_ak = Topt(P, Po, To, DCp, DHo, DVo, Da(3), Dk(2), R);

%-----Create the Figure
%figure('Position',[100 100 700 400])
figure(8)
clf
tiledlayout(1, 2, 'TileSpacing', 'compact','Padding','compact')

plt1 = nexttile();
hold on
title('A')
plt1.TitleHorizontalAlignment = 'left';
plot(T-273.15, k_T_1000, 'k','LineWidth',1)
plot(T_opt_1000-273.15,max(k_T_1000),'ko', 'MarkerSize', 5, 'LineWidth',1)
plot(T-273.15, k_T_100, 'k','LineWidth',1)
plot(T_opt_100-273.15, max(k_T_100), 'ksquare', 'MarkerSize', 5, 'LineWidth',1)
plot(T-273.15, k_T_1, 'k','LineWidth',1)
plot(T_opt_1-273.15, max(k_T_1), 'kdiamond', 'MarkerSize', 5, 'LineWidth',1)
xlabel(['Temperature (' char(176) 'C)'])
ylabel('Rate Coefficient')
legend('','Mariana Trench','','Average Ocean Depth','','Atmosphere', 'Location', 'northwest')
legend boxoff
box off
hold off

plt2 = nexttile();
hold on
title('B')
plt2.TitleHorizontalAlignment = 'left';
plot(P/101325, T_opt_0-273.15,'-k','LineWidth',1,'color',[0 0 0])
plot(P/101325, T_opt_a-273.15,'-k','LineWidth',1,'Color',[0 0 1])
plot(P/101325, T_opt_k-273.15, '-k','LineWidth',1,'Color',[0 1 0])
plot(P/101325, T_opt_ak-273.15, '-k','LineWidth',1,'Color',[0 1 1])
xlabel('Pressure (atm)')
ylabel(['Optimal Temperature (' char(176) 'C)'])
legend('\alpha = \kappa = Low',...
    '\alpha = High \kappa = Low',...
    '\alpha = Low \kappa = High',...
    '\alpha = \kappa = High','Location','northwest')
legend box off
box off
hold off

%----------Figure 4------------------------------------------------------%

% Low, Parameter, Optimal, Supraoptimal Temperatures
T_vals = [10;20;58;85]+273.15; 
T_a = P/101325/200+T_vals;

%-----Without adiabatic heating
DG_10 = DG(P, T_vals(1), Po, To, DCp, DSo, DHo, DVo, Da(1), Dk(1));
DG_20 = DG(P, T_vals(2), Po, To, DCp, DSo, DHo, DVo, Da(1), Dk(1));
DG_58 = DG(P, T_vals(3), Po, To, DCp, DSo, DHo, DVo, Da(1), Dk(1));
DG_85 = DG(P, T_vals(4), Po, To, DCp, DSo, DHo, DVo, Da(1), Dk(1));
k_10 = k(T_vals(1), DG_10, kb, h, R);
k_20 = k(T_vals(2), DG_20, kb, h, R);
k_58 = k(T_vals(3), DG_58, kb, h, R);
k_85 = k(T_vals(4), DG_85, kb, h, R);

%-----With adiabatic heating
DG_10_a = DG(P, T_a(1,:), Po, To, DCp, DSo, DHo, DVo, Da(1), Dk(1));
DG_20_a = DG(P, T_a(2,:), Po, To, DCp, DSo, DHo, DVo, Da(1), Dk(1));
DG_58_a = DG(P, T_a(3,:), Po, To, DCp, DSo, DHo, DVo, Da(1), Dk(1));
DG_85_a = DG(P, T_a(4,:), Po, To, DCp, DSo, DHo, DVo, Da(1), Dk(1));
k_10_a = k(T_a(1,:), DG_10_a, kb, h, R);
k_20_a = k(T_a(2,:), DG_20_a, kb, h, R);
k_58_a = k(T_a(3,:), DG_58_a, kb, h, R);
k_85_a = k(T_a(4,:), DG_85_a, kb, h, R);

%-----Create the Figure
figure('Position',[100 100 400 400])\data
tiledlayout(2, 2, 'TileSpacing', 'compact','Padding','compact')

nexttile
hold on
title(['10 ' char(176) 'C'])
plot(P/101325,k_10,'k','LineWidth',1)
plot(P/101325,k_10_a,'--k','LineWidth',1)
xlabel('Pressure (atm)')
ylabel('Rate Coefficient')
xlim([0 1000])
ylim([k_10(1)-0.5*k_10(1) k_10(1)+0.5*k_10(1)])
hold off

nexttile
hold on
title(['20 ' char(176) 'C'])
plot(P/101325,k_20,'k','LineWidth',1)
plot(P/101325,k_20_a,'--k','LineWidth',1)
xlabel('Pressure (atm)')
ylabel('Rate Coefficient')
xlim([0 1000])
ylim([k_20(1)-0.5*k_20(1) k_20(1)+0.5*k_20(1)])
hold off

nexttile
hold on
title(['58 ' char(176) 'C'])
plot(P/101325,k_58,'k','LineWidth',1)
plot(P/101325,k_58_a,'--k','LineWidth',1)
xlabel('Pressure (atm)')
ylabel('Rate Coefficient')
xlim([0 1000])
ylim([k_58(1)-0.5*k_58(1) k_58(1)+0.5*k_58(1)])
legend('Isothermal','Real Heating','Location','northeast')
legend boxoff
hold off

nexttile
hold on
title(['85 ' char(176) 'C'])
plot(P/101325,k_85,'k','LineWidth',1)
plot(P/101325,k_85_a,'--k','LineWidth',1)
xlabel('Pressure (atm)')
ylabel('Rate Coefficient')
xlim([0 1000])
ylim([k_85(1)-0.5*k_85(1) k_85(1)+0.5*k_85(1)])
hold off

%----------Figure 5------------------------------------------------------%

DG_PT = DG(P, T, Po, To, DCp, DSo, DHo, DVo, Da(1), Dk(1));
k_PT = k(T, DG_PT, kb, h, R);

figure('Position',[100 100 400 400])
surf(P/101325,T-273.15,k_PT)
shading interp
colormap('turbo')
grid minor
xlabel('Pressure (atm)')
ylabel(['Temperature (' char(176) 'C)'])
zlabel('Rate Coefficient')
zlim([0 15])
caxis([0 15])

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

%----------Optimal Temperature-------------------------------------------%

function Topt = Topt(P, Po, To, DCp, DHo, DVo, Da, Dk, R)
% Calculate the optimal tempertature as a function of pressure
    Topt = (DCp*To-DHo+Da*To*(P-Po)-DVo*(P-Po)+Dk/2*(P-Po).^2)/(DCp+R);
end

%------------------------------------------------------------------------%