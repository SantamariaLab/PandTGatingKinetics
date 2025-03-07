%% Hodgkin-Huxley Model %%
% Jake Miller %
% 3/30/2023 %

%----------Constants-----------------------------------------------------%

kb = 1.380649e-23;             % Boltzmann constant (m^2 kg / s^2 K)
hplanck = 6.62607015e-34;      % Planck's constant (m^2 kg / s)
R= 8.3145;                     % Ideal gas constant (J / mol K)

%----------Thermodynamics Parameters-------------------------------------%

DCp = [-890, -2820];           % Change in heat capacity (J / mol K)
DSo = [-122.2, -51.44];        % Change in reference entropy (J / mol K)
DHo = [33050, 51510];          % Change in reference enthalpy (J / mol)
DVo = 36e-6;                   % Change in reference volume (m^3 / mol)
To = 298.15;                   % Reference temperature (K)
Po = 101325;                   % Reference pressure (Pa)
Da = 1e-8;                     % Change in expansivity (m^3 / mol K)
Dk = 1e-13;                    % Change in compressibility (m^3 / mol Pa)

%----------Hodgkin-Huxley Parameters-------------------------------------%

Vo = -65;                      % Equilibrium potential (mV)
Vq = -20;                      % Qualifying potential (mV)
n = 0.3177;                    % Potassium activation variable
m = 0.0529;                    % Sodium activation variable
h = 0.5960;                    % Sodium inactivation variable
gl = 0.3;                      % Leak conductance (mS / cm^2)
gK = 36;                       % Potassium conductance (mS / cm^2)
gNa = 120;                     % Sodium conductance (mS / cm^2)
El = -54.4;                    % Leak reversal potential (mV)
EK = -77;                      % Potassium reversal potential (mV)
ENa = 50;                      % Sodium reversal potential (mV)
I = 0;                         % Initial input current (uA / cm^2)
Cm = 1;                        % Membrane capacitance (uF / cm^2)

%----------Rate Calculation----------------------------------------------%

P = [1,25,625]*Po;       % Experimental pressures
T = [0,6.3,15]+273.15;   % Experimental temperatures

DG_11 = DG(P(1), T(1), Po, To, DCp, DSo, DHo, DVo, Da, Dk);
k_11 = k(T(1), DG_11, kb, hplanck, R);
DG_21 = DG(P(2), T(1), Po, To, DCp, DSo, DHo, DVo, Da, Dk);
k_21 = k(T(1), DG_21, kb, hplanck, R);
DG_31 = DG(P(3), T(1), Po, To, DCp, DSo, DHo, DVo, Da, Dk);
k_31 = k(T(1), DG_31, kb, hplanck, R);

DG_12 = DG(P(1), T(2), Po, To, DCp, DSo, DHo, DVo, Da, Dk);
k_12 = k(T(2), DG_12, kb, hplanck, R);
DG_22 = DG(P(2), T(2), Po, To, DCp, DSo, DHo, DVo, Da, Dk);
k_22 = k(T(2), DG_22, kb, hplanck, R);
DG_32 = DG(P(3), T(2), Po, To, DCp, DSo, DHo, DVo, Da, Dk);
k_32 = k(T(2), DG_32, kb, hplanck, R);

DG_13 = DG(P(1), T(3), Po, To, DCp, DSo, DHo, DVo, Da, Dk);
k_13 = k(T(3), DG_13, kb, hplanck, R);
DG_23 = DG(P(2), T(3), Po, To, DCp, DSo, DHo, DVo, Da, Dk);
k_23 = k(T(3), DG_23, kb, hplanck, R);
DG_33 = DG(P(3), T(3), Po, To, DCp, DSo, DHo, DVo, Da, Dk);
k_33 = k(T(3), DG_33, kb, hplanck, R);

%----------Simulation----------------------------------------------------%

% Single Action Potential
Time = 30;                         % Run time (ms)

AP_11 = HH(k_11, n, m, h, gl, gK, gNa, Vo, El, EK, ENa, I, Cm, Vq, Time, true);
AP_21 = HH(k_21, n, m, h, gl, gK, gNa, Vo, El, EK, ENa, I, Cm, Vq, Time, true);
AP_31 = HH(k_31, n, m, h, gl, gK, gNa, Vo, El, EK, ENa, I, Cm, Vq, Time, true);
AP_12 = HH(k_12, n, m, h, gl, gK, gNa, Vo, El, EK, ENa, I, Cm, Vq, Time, true);
AP_22 = HH(k_22, n, m, h, gl, gK, gNa, Vo, El, EK, ENa, I, Cm, Vq, Time, true);
AP_32 = HH(k_32, n, m, h, gl, gK, gNa, Vo, El, EK, ENa, I, Cm, Vq, Time, true);
AP_13 = HH(k_13, n, m, h, gl, gK, gNa, Vo, El, EK, ENa, I, Cm, Vq, Time, true);
AP_23 = HH(k_23, n, m, h, gl, gK, gNa, Vo, El, EK, ENa, I, Cm, Vq, Time, true);
AP_33 = HH(k_33, n, m, h, gl, gK, gNa, Vo, El, EK, ENa, I, Cm, Vq, Time, true);

figure('Position',[100 100 900 600])
tiledlayout(3, 5, 'TileSpacing', 'compact','Padding','compact')

plt1ap = nexttile(1);
hold on
title(['0 ' char(176) 'C'])
plot(AP_11(1,:),AP_11(2,:),'k','LineWidth',1)
plot(AP_21(1,:),AP_21(2,:),':k','LineWidth',1)
plot(AP_31(1,:),AP_31(2,:),'--k','LineWidth',1)
xlabel('Time (ms)')
xlim([0 30])
ylabel('Membrane Potential (mV)')
ylim([-100 50])
hold off

plt2ap = nexttile(6);
hold on
title(['6.3 ' char(176) 'C'])
plot(AP_12(1,:),AP_12(2,:),'k','LineWidth',1)
plot(AP_22(1,:),AP_22(2,:),':k','LineWidth',1)
plot(AP_32(1,:),AP_32(2,:),'--k','LineWidth',1)
xlabel('Time (ms)')
xlim([0 30])
ylabel('Membrane Potential (mV)')
ylim([-100 50])
hold off

plt3ap = nexttile(11);
hold on
title(['15 ' char(176) 'C'])
plot(AP_13(1,:),AP_13(2,:),'k','LineWidth',1)
plot(AP_23(1,:),AP_23(2,:),':k','LineWidth',1)
plot(AP_33(1,:),AP_33(2,:),'--k','LineWidth',1)
xlabel('Time (ms)')
xlim([0 30])
ylabel('Membrane Potential (mV)')
ylim([-100 50])
hold off

% Continuous Spiking
I = 16;
Time = 600;                         % Run time (ms)

Spikes_11 = HH(k_11, n, m, h, gl, gK, gNa, Vo, El, EK, ENa, I, Cm, Vq, Time, false);
Spikes_21 = HH(k_21, n, m, h, gl, gK, gNa, Vo, El, EK, ENa, I, Cm, Vq, Time, false);
Spikes_31 = HH(k_31, n, m, h, gl, gK, gNa, Vo, El, EK, ENa, I, Cm, Vq, Time, false);
Spikes_12 = HH(k_12, n, m, h, gl, gK, gNa, Vo, El, EK, ENa, I, Cm, Vq, Time, false);
Spikes_22 = HH(k_22, n, m, h, gl, gK, gNa, Vo, El, EK, ENa, I, Cm, Vq, Time, false);
Spikes_32 = HH(k_32, n, m, h, gl, gK, gNa, Vo, El, EK, ENa, I, Cm, Vq, Time, false);
Spikes_13 = HH(k_13, n, m, h, gl, gK, gNa, Vo, El, EK, ENa, I, Cm, Vq, Time, false);
Spikes_23 = HH(k_23, n, m, h, gl, gK, gNa, Vo, El, EK, ENa, I, Cm, Vq, Time, false);
Spikes_33 = HH(k_33, n, m, h, gl, gK, gNa, Vo, El, EK, ENa, I, Cm, Vq, Time, false);

plt1 = nexttile(2, [1 4]);
hold on
plot(Spikes_11(1,:),Spikes_11(2,:),'k','LineWidth',0.75)
plot(Spikes_21(1,:),Spikes_21(2,:),':k','LineWidth',0.75)
plot(Spikes_31(1,:),Spikes_31(2,:),'--k','LineWidth',0.75)
xlabel('Time (ms)')
xlim([0 150])
ylim([-100 50])
hold off

plt2 = nexttile(7, [1 4]);
hold on
plot(Spikes_12(1,:),Spikes_12(2,:),'k','LineWidth',0.75)
plot(Spikes_22(1,:),Spikes_22(2,:),':k','LineWidth',0.75)
plot(Spikes_32(1,:),Spikes_32(2,:),'--k','LineWidth',0.75)
xlabel('Time (ms)')
xlim([0 150])
ylim([-100 50])
legend('1 atm','25 atm','625 atm','Location','eastoutside')
legend boxoff
hold off

plt3 = nexttile(12, [1 4]);
hold on
plot(Spikes_13(1,:),Spikes_13(2,:),'k','LineWidth',0.75)
plot(Spikes_23(1,:),Spikes_23(2,:),':k','LineWidth',0.75)
plot(Spikes_33(1,:),Spikes_33(2,:),'--k','LineWidth',0.75)
xlabel('Time (ms)')
xlim([0 150])
ylim([-100 50])
hold off

%---------Firing Rate Calculation----------------------------------------%

fr_11 = fr(Spikes_11(2,:),10000);     % 10000 --> start at 100 ms                     
fr_21 = fr(Spikes_21(2,:),10000);
fr_31 = fr(Spikes_31(2,:),10000);
fr_12 = fr(Spikes_12(2,:),10000);
fr_22 = fr(Spikes_22(2,:),10000);
fr_32 = fr(Spikes_32(2,:),10000);
fr_13 = fr(Spikes_13(2,:),10000);
fr_23 = fr(Spikes_23(2,:),10000);
fr_33 = fr(Spikes_33(2,:),10000);

frs = [fr_11 fr_12 fr_13;    % Columns = 0 6.3 15 Celsius
       fr_21 fr_22 fr_23;    % Rows = 1 25 625 atm
       fr_31 fr_32 fr_33]

%----------Firing Rate Curves--------------------------------------------%

I = 100;
P_range = transpose([1:160]*5*Po);
T_range = transpose([-20:80]+273.15);

fr_P = ones(length(P_range),1);
fr_T = ones(length(T_range),1);

figure('Position',[100 100 900 300])
tiledlayout(1, 3, 'TileSpacing', 'compact','Padding','compact')

% Pressure Curve
nexttile()
for j = (1:3)
    DG_P = DG(P_range, T(j), Po, To, DCp, DSo, DHo, DVo, Da, Dk);
    k_P = k(T(j), DG_P, kb, hplanck, R);
    for i = (1:length(k_P))
        spikes = HH(k_P(i,:), n, m, h, gl, gK, gNa, Vo, El, EK, ENa, I, Cm, Vq, Time, false);
        fr_P(i) = fr(spikes(2,:),10000);
    end
    fr_Ps(:,j) = fr_P;
end

hold on
%plot(P_range/Po,fr_Ps(:,1), 'Color', [0.4 0.4 0.4], 'LineWidth', 2)
%plot(P_range/Po,fr_Ps(:,2), 'k', 'LineWidth', 1)    % T = 6.3 C
plot(P_range/Po,fr_Ps(:,3), 'k', 'LineWidth', 1)     % T = 15 C
xlabel('Pressure (atm)')
xlim([0 800])
ylabel('Firing Rate (Hz)')
ylim([0 605])
box off
hold off

% Temperature Curve
nexttile()
for j = (1:3)
    DG_T = DG(P(j), T_range, Po, To, DCp, DSo, DHo, DVo, Da, Dk);
    k_T = k(T_range, DG_T, kb, hplanck, R);
    for i = (1:length(k_T))
        spikes = HH(k_T(i,:), n, m, h, gl, gK, gNa, Vo, El, EK, ENa, I, Cm, Vq, Time, false);
        fr_T(i) = fr(spikes(2,:),10000);
    end
    fr_Ts(:,j) = fr_T;
end

hold on
%plot(T_range-273.15,fr_Ts(:,1), 'k', 'LineWidth', 1)
%plot(T_range-273.15,fr_Ts(:,2), 'k', 'LineWidth', 1)    % P = 25 atm
plot(T_range-273.15,fr_Ts(:,3), 'k', 'LineWidth', 1)     % P = 625 atm
xlabel(['Temperature (' char(176) 'C)'])
xlim([-20 80])
ylabel('Firing Rate (Hz)')
ylim([0 605])
box off
hold off

% Current Curve
I_range = 0:225;
fr_I = ones(length(I_range),1);

for i = (1:(length(I_range)))
    spikes = HH(k_33, n, m, h, gl, gK, gNa, Vo, El, EK, ENa, I_range(i), Cm, Vq, Time, false);
    fr_I(i) = fr(spikes(2,:),10000);
end

nexttile()
plot(I_range,fr_I, 'k', 'LineWidth', 1)
xlabel('Input Current \muA/cm^2')
xlim([0 length(I_range)-1])
ylabel('Firing Rate (Hz)')
ylim([0 605])
box off

%---------Steady-State Sodium Inactivation-------------------------------%

V = [-100:-20];
[~, ~, ~, ~, ah_p, bh_p] = rates(V, Vo);
hinf = ah_p./(ah_p + bh_p);
sig = 1./(1+exp((V+52)/5));
sig_62 = 1./(1+exp((V+54)/7));

% figure(3)
% hold on
% plot(V,hinf, 'LineWidth', 1)
% plot(V,sig,'LineWidth', 1)
% plot(V,sig_62,'LineWidth', 1)
% hold off

%---------Hodgkin-Huxley Model-------------------------------------------%

function AP = HH(k, n0, m0, h0, gl, gK, gNa, Vo, El, EK, ENa, I, Cm, Vq, T, useMinCurrent)
% Model a single action potential using the HH model
    iters = 0;
    if useMinCurrent == true
        I = 0;
    end
    while iters < 500
        n = n0;
        m = m0;
        h = h0;
        dt = 0.01;               % Time step (ms)
        t = [1:100*T];           % Time indices
        time = t*dt;             % Time (ms)
        V = 0*ones(1,length(t)); 
        V(1) = Vo;               % Start potential (mV)
        for t = t(1:length(t)-1)
            % Potential step (mV)
            dV = -(gl*(V(t)-El)+gK*n^4*(V(t)-EK)+gNa*m^3*h*(V(t)-ENa)-I)*dt/Cm;
            V(1,t+1) = V(1,t) + dV;     
     
            % Rates
            [an, bn, am, bm, ah, bh] = rates(V(1,t), Vo);
            n = dt*k(1,2)*(an*(1-n)-bn*n) + n;
            m = dt*k(1,1)*(am*(1-m)-bm*m) + m;
            h = dt*mean(k)*(ah*(1-h)-bh*h) + h;
        end
        if useMinCurrent == true
            if max(V) > Vq
                AP = [time; V];          % Return action potential
                I
                return
            end
            I = I + 0.1;
            iters = iters + 1;
        else
            AP = [time; V];
            return
        end
    end
end

%----------Forward and Backward Rates------------------------------------%

function [an, bn, am, bm, ah, bh] = rates(V, Vo)
% Calculate the n, m, and h forward and backward rates

% K activation
an = (0.1-0.01*(V-Vo))./(exp(1-0.1*(V-Vo))-1);
bn = 0.125*exp(-(V-Vo)./80);

% Na activation
am = (2.5-0.1*(V-Vo))./(exp(2.5-0.1*(V-Vo))-1);
bm = 4*exp(-(V-Vo)./18);

% Na inactivation
ah = 0.07*exp(-(V-Vo)./20);
bh = 1./(1+exp(3-0.1*(V-Vo)));
end

%----------Rate Coefficient----------------------------------------------%

function k = k(T, G, kb, h, R)
% Calculate the rate coefficient
    k = (kb*T/h).*exp(-G./(R*T));
end

%----------Energy Barrier------------------------------------------------%

function DG = DG(P, T, Po, To, DCp, DSo, DHo, DVo, Da, Dk)
% Calculate the change in Gibbs free energy across activation barrier
    DG = DCp.*(T-To)-DCp.*T.*log(T/To)-T.*DSo+Da*(T-To).*(P-Po)+DVo*(P-Po)-Dk/2*(P-Po).^2+DHo;
end

%----------Firing Rate---------------------------------------------------%

function fr = fr(spikes, i)
    time = (length(spikes)-i)/100000; % indices to seconds (s)
    spike_number = length(findpeaks(spikes(i:length(spikes))));
    fr = spike_number/time;
end

%------------------------------------------------------------------------%