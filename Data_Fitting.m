%% Fitting to Data %%
% Jake Miller %
% 4/10/2023 %

%----------Loading Data--------------------------------------------------%

data = xlsread('pressure_data/Pressure Data.xlsx',2); 
k_data = data(:,1);                                 % Rate coefficient
P_data = data(:,2);                                 % Pressure (MPa)
T_data = data(:,3);                                 % Temperature (C)

% Data from Conti 1982 Potassium Experiment %
k_ContiK = k_data(1:9);
%k_ContiK = k_ContiK/min(k_ContiK);
P_ContiK = P_data(1:9);
T_ContiK = T_data(1:9);

% Data from Conti 1982 Sodium Experiment %
k_ContiNa = k_data(10:24);
%k_ContiNa = k_ContiNa/min(k_ContiNa);
P_ContiNa = P_data(10:24);
T_ContiNa = T_data(10:24);

% Data from Heinemann 1987 Acetylcholine Experiment %
k_Heine = k_data(25:41);
%k_Heine = k_Heine/min(k_Heine);
P_Heine = P_data(25:41);
T_Heine = T_data(25:41);

% Data from Spyropoulos 1957 Toad Nerve Experiment %
k_Spyro = k_data(42:64);
%k_Spyro = k_Spyro/min(k_Spyro);
P_Spyro = P_data(42:64);
T_Spyro = T_data(42:64);

% Data from Conti 1984 Sodium Experiment %
k_Conti84 = k_data(65:71);
P_Conti84 = P_data(65:71);
T_Conti84 = T_data(65:71);

% Data from Meyer 1997 Potassium Experiment %
k_Meyer = k_data(72:77);
P_Meyer = P_data(72:77);
T_Meyer = T_data(72:77);

% Plot the data %
figure(1)
hold on
scatter3(P_ContiK, T_ContiK, k_ContiK, 80, 'ok')
scatter3(P_ContiNa, T_ContiNa, k_ContiNa, 80, '+k')
scatter3(P_Heine, T_Heine, k_Heine, 80, 'diamondk')
scatter3(P_Spyro, T_Spyro, k_Spyro, 100, '.k')
scatter3(P_Conti84, T_Conti84, k_Conti84, 80, 'xk')
scatter3(P_Meyer, T_Meyer, k_Meyer, 80, 'squarek')
legend('Conti K','Conti Na','Heinemann','Spyropoulos','Conti 1984','Meyer')
legend boxoff
hold off

%----------Pressure Fits-------------------------------------------------%
%GENERATES FIGURE 2

[ContiK_fit_P, ContiK_gof] = fit_to_pressure(P_ContiK, k_ContiK);
[ContiNa_fit_P, ContiNa_gof] = fit_to_pressure(P_ContiNa, k_ContiNa);
[Heine_fit_P, Heine_gof] = fit_to_pressure(P_Heine, k_Heine);
[Spyro_fit_P, Spyro_gof] = fit_to_pressure(P_Spyro, k_Spyro);
[Conti84_fit_P, Conti84_gof] = fit_to_pressure(P_Conti84, k_Conti84);
[Meyer_fit_P, Meyer_gof] = fit_to_pressure(P_Meyer, k_Meyer);

Rsq = [ContiK_gof.rsquare ContiNa_gof.rsquare Heine_gof.rsquare ...
    Spyro_gof.rsquare Conti84_gof.rsquare Meyer_gof.rsquare];
average_Rsq = mean(Rsq)
std_Rsq = std(Rsq)

Rsq_no_Heine = [Rsq(1:2) Rsq(4:6)];
average_Rsq_no_Heine = mean(Rsq_no_Heine)
std_Rsq_no_Heine = std(Rsq_no_Heine)

figure('Position',[100 100 400 300])
hold on

ContiK_plot = plot(ContiK_fit_P, 'k', P_ContiK, k_ContiK, 'ok');
set(ContiK_plot, 'LineWidth', 1, 'MarkerSize', 5)

ContiNa_plot = plot(ContiNa_fit_P, 'k', P_ContiNa, k_ContiNa, '+k');
set(ContiNa_plot, 'LineWidth', 1, 'MarkerSize', 5)

Heine_plot = plot(Heine_fit_P, 'k', P_Heine, k_Heine, 'diamondk');
set(Heine_plot, 'LineWidth', 1, 'MarkerSize', 5)

Spyro_plot = plot(Spyro_fit_P, 'k', P_Spyro, k_Spyro, '.k');
set(Spyro_plot, 'LineWidth', 1, 'MarkerSize', 20)

Conti84_plot = plot(Conti84_fit_P, 'k', P_Conti84, k_Conti84, 'xk');
set(Conti84_plot, 'LineWidth', 1, 'MarkerSize', 5)

Meyer_plot = plot(Meyer_fit_P, 'k', P_Meyer, k_Meyer, 'squarek');
set(Meyer_plot, 'LineWidth', 1, 'MarkerSize', 5)

xlabel('Pressure (MPa)')
xlim([0 72])
ylabel('Rate Coefficient')
ylim([0.2 1.1])
legend off
hold off

%----------Temperature Fits----------------------------------------------%

% Fit to Combined Conti Data
k_Conti = k_data(1:24);
P_Conti = P_data(1:24);
T_Conti = T_data(1:24);

P_21_i = P_Conti == 21;                         % Indices where P = 21 MPa
k_Conti_21 = k_Conti(P_21_i);
T_Conti_21 = T_Conti(P_21_i);

P_41_i = (P_Conti == 41) | (P_Conti == 42);     % Indices where P = 41 MPa
k_Conti_41 = k_Conti(P_41_i);
T_Conti_41 = T_Conti(P_41_i);

P_62_i = P_Conti == 62;                         % Indices where P = 62 MPa
k_Conti_62 = k_Conti(P_62_i);
T_Conti_62 = T_Conti(P_62_i);

[Conti_fit_T_21, Conti_gof_T_21] = fit_to_temperature(T_Conti_21, k_Conti_21);
[Conti_fit_T_41, Conti_gof_T_41] = fit_to_temperature(T_Conti_41, k_Conti_41);
[Conti_fit_T_62, Conti_gof_T_62] = fit_to_temperature(T_Conti_62, k_Conti_62);

% Fit to Heinemann Data
P_40_i = (P_Heine >= 38) & (P_Heine <= 42);     % Indices where P = 40 MPa
k_Heine_40 = k_Heine(P_40_i);
T_Heine_40 = T_Heine(P_40_i);
[Heine_fit_T_40, Heine_gof_T_40] = fit_to_temperature(T_Heine_40, k_Heine_40);

%----------Conti K Pressure Plot---------------------------------------%
% GENERATES FIGURE 3

t_data = [1.36 2.02 3.30;  % 5 C      21 MPa  41 MPa  62 MPa
          1.28 1.83 2.67;  % 10 C
          1.27 1.69 2.54;  % 15 C
          1.31 1.91 2.83]; % Mean
k_data = 1./t_data; % Convert to rate
P_vals = [21, 41, 62];
T_vals = [5, 10, 15];

k_5 = k_data(1,:);        % Rate at 5 C
k_10 = k_data(2,:);       % Rate at 10 C
k_15 = k_data(3,:);       % Rate at 15 C

% First fit to rate at 10 C using DVo, DHo, DSo with Da=Dk=0 and T=To=10 C
% Then use that DVo to fit at 5 and 15 C with DHo, DSo, DCp
% Then fix DVo, DCp, DSo and fit with DHo and Da

[Conti_K_fits, Conti_K_gof] = conti_K_fits(P_vals,k_5,k_10,k_15);
fit_10 = Conti_K_fits{1};
fit_5 = Conti_K_fits{2};
fit_5_alpha = Conti_K_fits{3};
fit_15 = Conti_K_fits{4};
fit_15_alpha = Conti_K_fits{5};

figure('Position',[100 100 400 600])
tiledlayout(2,1, 'TileSpacing', 'compact','Padding','compact')

nexttile
hold on
fit15 = plot(fit_15, 'k', P_vals, k_15, 'ok');
fit10 = plot(fit_10, 'k', P_vals, k_10, '+k');
fit5 = plot(fit_5, 'k', P_vals, k_5, 'xk');
set(fit15, 'LineWidth', 1, 'MarkerSize', 5)
set(fit10, 'LineWidth', 1, 'MarkerSize', 5)
set(fit5, 'LineWidth', 1, 'MarkerSize', 5)
xlim([18 65])
ylim([0.25 0.85])
xlabel('')
ylabel('Rate Coefficient')
title('MMRT Fit')
legend(['15 ' char(176) 'C'],'',['10 ' char(176) 'C'],'',['5 ' char(176) 'C'],'Location','northeast')
legend boxoff
box off
hold off

nexttile
hold on
fit15a = plot(fit_15_alpha, 'k', P_vals, k_15, 'ok');
fit10a = plot(fit_10, 'k', P_vals, k_10, '+k');
fit5a = plot(fit_5_alpha, 'k', P_vals, k_5, 'xk');
set(fit15a, 'LineWidth', 1, 'MarkerSize', 5)
set(fit10a, 'LineWidth', 1, 'MarkerSize', 5)
set(fit5a, 'LineWidth', 1, 'MarkerSize', 5)
xlim([18 65])
xlabel('Pressure (MPa)')
ylim([0.25 0.85])
ylabel('Rate Coefficient')
title('Expansivity Fit')
legend off
box off
hold off




