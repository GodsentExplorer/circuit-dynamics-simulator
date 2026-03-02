%% RC/RL/LC/RLC Circuit Dynamics Simulator
% Intermediate Electromagnetism
% Tools: MATLAB, ode45
% Mirrors: python/circuit_dynamics.ipynb

clear; clc; close all;

%% 1. RC Circuit — Charging & Discharging
% Governing ODE: dVc/dt = (Vs - Vc) / (R*C)
% Analytical:    Vc(t) = Vs * (1 - exp(-t/tau))

% Parameters
R   = 1e3;       % Resistance (Ohms) — 1 kΩ
C   = 1e-6;      % Capacitance (F)   — 1 µF
Vs  = 5.0;       % Source voltage (V)
tau = R * C;     % Time constant (s)

t_span = [0, 5*tau];
t_eval = linspace(0, 5*tau, 500);

% ODE function
rc_ode = @(t, y) (Vs - y(1)) / (R * C);

% Solve numerically
[t_num, Vc_num] = ode45(rc_ode, t_span, 0);

% Analytical solution
Vc_analytical = Vs * (1 - exp(-t_eval / tau));

% Current
I_num = (Vs - Vc_num(:,1)) / R;

% ── Plot ─────────────────────────────────────────────────────────────────────
figure('Position', [100 100 1200 400]);

subplot(1,2,1)
plot(t_num*1e3, Vc_num(:,1), 'b-', 'LineWidth', 2, 'DisplayName', 'Numerical (ode45)'); hold on;
plot(t_eval*1e3, Vc_analytical, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Analytical');
yline(Vs, 'k:', 'DisplayName', sprintf('Vs = %.1fV', Vs));
xline(tau*1e3, 'g--', 'DisplayName', sprintf('\\tau = %.2f ms', tau*1e3));
xlabel('Time (ms)'); ylabel('Voltage (V)');
title('RC Circuit — Capacitor Voltage');
legend; grid on;

subplot(1,2,2)
plot(t_num*1e3, I_num*1e3, 'Color', [1 0.5 0], 'LineWidth', 2, 'DisplayName', 'Current');
xline(tau*1e3, 'g--', 'DisplayName', sprintf('\\tau = %.2f ms', tau*1e3));
xlabel('Time (ms)'); ylabel('Current (mA)');
title('RC Circuit — Charging Current');
legend; grid on;

sgtitle(sprintf('RC Circuit | R=%.0fk\\Omega  C=%.0f\\muF  \\tau=%.2fms', ...
    R/1e3, C*1e6, tau*1e3), 'FontSize', 13);

fprintf('τ = %.2f ms | At t=τ, Vc = %.3fV (should be ~%.3fV)\n', ...
    tau*1e3, Vs*(1-exp(-1)), 0.632*Vs);

%% 2. RL Circuit — Current Buildup & Decay
% Governing ODE: dI/dt = (Vs - I*R) / L
% Analytical:    I(t) = (Vs/R) * (1 - exp(-t/tau))

% Parameters
R_rl  = 1e3;          % Resistance (Ohms) — 1 kΩ
L     = 1.0;          % Inductance (H)    — 1 H
Vs_rl = 5.0;          % Source voltage (V)
tau_rl = L / R_rl;    % Time constant (s)
I_ss   = Vs_rl / R_rl; % Steady-state current

t_span_rl = [0, 5*tau_rl];
t_eval_rl = linspace(0, 5*tau_rl, 500);

% ODE function
rl_ode = @(t, y) (Vs_rl - y(1) * R_rl) / L;

% Solve numerically
[t_num_rl, I_num] = ode45(rl_ode, t_span_rl, 0);

% Analytical solution
I_analytical = (Vs_rl / R_rl) * (1 - exp(-t_eval_rl / tau_rl));

% Inductor voltage
VL = Vs_rl * exp(-t_eval_rl / tau_rl);

% ── Plot ─────────────────────────────────────────────────────────────────────
figure('Position', [100 100 1200 400]);

subplot(1,2,1)
plot(t_num_rl*1e3, I_num(:,1)*1e3, 'b-', 'LineWidth', 2, 'DisplayName', 'Numerical (ode45)'); hold on;
plot(t_eval_rl*1e3, I_analytical*1e3, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Analytical');
yline(I_ss*1e3, 'k:', 'DisplayName', sprintf('I_{ss} = %.2f mA', I_ss*1e3));
xline(tau_rl*1e3, 'g--', 'DisplayName', sprintf('\\tau = %.3f ms', tau_rl*1e3));
xlabel('Time (ms)'); ylabel('Current (mA)');
title('RL Circuit — Inductor Current');
legend; grid on;

subplot(1,2,2)
plot(t_eval_rl*1e3, VL, 'Color', [1 0.5 0], 'LineWidth', 2, 'DisplayName', 'V_L');
xline(tau_rl*1e3, 'g--', 'DisplayName', sprintf('\\tau = %.3f ms', tau_rl*1e3));
xlabel('Time (ms)'); ylabel('Voltage (V)');
title('RL Circuit — Inductor Voltage Decay');
legend; grid on;

sgtitle(sprintf('RL Circuit | R=%.0fk\\Omega  L=%.1fH  \\tau=%.3fms', ...
    R_rl/1e3, L, tau_rl*1e3), 'FontSize', 13);

fprintf('τ = %.3f ms | I_ss = %.2f mA | At t=τ, I = %.3f mA\n', ...
    tau_rl*1e3, I_ss*1e3, I_ss*(1-exp(-1))*1e3);

%% 3. LC Circuit — Undamped Oscillations
% System:  dVc/dt = I/C
%          dI/dt  = -Vc/L
% Analytical: Vc(t) = Vc0*cos(w0*t), I(t) = -Vc0*sqrt(C/L)*sin(w0*t)

% Parameters
L_lc  = 1.0;     % Inductance (H)
C_lc  = 1e-6;    % Capacitance (F)
Vc0   = 5.0;     % Initial voltage (V)
I0    = 0.0;     % Initial current (A)

omega0 = 1 / sqrt(L_lc * C_lc);   % Natural frequency (rad/s)
f0     = omega0 / (2*pi);          % Natural frequency (Hz)
T      = 1 / f0;                   % Period (s)

t_span_lc = [0, 3*T];
t_eval_lc = linspace(0, 3*T, 1000);

% ODE system
lc_ode = @(t, y) [y(2)/C_lc; -y(1)/L_lc];

% Solve numerically
[t_num_lc, sol_lc] = ode45(lc_ode, t_span_lc, [Vc0; I0]);
Vc_lc = sol_lc(:,1);
I_lc  = sol_lc(:,2);

% Analytical solution
Vc_analytical_lc = Vc0 * cos(omega0 * t_eval_lc);
I_analytical_lc  = -Vc0 * sqrt(C_lc/L_lc) * sin(omega0 * t_eval_lc);

% ── Plot ─────────────────────────────────────────────────────────────────────
figure('Position', [100 100 1600 400]);

subplot(1,3,1)
plot(t_num_lc*1e3, Vc_lc, 'b-', 'LineWidth', 2, 'DisplayName', 'Numerical (ode45)'); hold on;
plot(t_eval_lc*1e3, Vc_analytical_lc, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Analytical');
xlabel('Time (ms)'); ylabel('Voltage (V)');
title('LC — Capacitor Voltage');
legend; grid on;

subplot(1,3,2)
plot(t_num_lc*1e3, I_lc*1e3, 'Color', [1 0.5 0], 'LineWidth', 2, 'DisplayName', 'Numerical (ode45)'); hold on;
plot(t_eval_lc*1e3, I_analytical_lc*1e3, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Analytical');
xlabel('Time (ms)'); ylabel('Current (mA)');
title('LC — Inductor Current');
legend; grid on;

subplot(1,3,3)
plot(Vc_lc, I_lc*1e3, 'Color', [0.5 0 0.8], 'LineWidth', 2);
xlabel('Voltage (V)'); ylabel('Current (mA)');
title('LC — Phase Portrait (V vs I)');
axis equal; grid on;

sgtitle(sprintf('LC Circuit | L=%.1fH  C=%.0f\\muF  f_0=%.2fHz  T=%.2fms', ...
    L_lc, C_lc*1e6, f0, T*1e3), 'FontSize', 13);

fprintf('ω₀ = %.2f rad/s | f₀ = %.2f Hz | T = %.2f ms\n', omega0, f0, T*1e3);

%% 4. RLC Circuit — Damped Oscillations
% System:  dVc/dt = I/C
%          dI/dt  = -Vc/L - (R/L)*I
% Damping ratio: zeta = (R/2)*sqrt(C/L)

% Fixed parameters
L_rlc  = 1.0;    % Inductance (H)
C_rlc  = 1e-6;   % Capacitance (F)
Vc0_rlc = 5.0;   % Initial voltage (V)

omega0_rlc = 1 / sqrt(L_rlc * C_rlc);
R_critical = 2 * sqrt(L_rlc / C_rlc);
T_rlc = 2*pi / omega0_rlc;

t_span_rlc = [0, 6*T_rlc];
t_eval_rlc = linspace(0, 6*T_rlc, 2000);

% Three damping cases
R_values = [R_critical*0.2, R_critical, R_critical*3.0];
labels   = {'Underdamped (\zeta < 1)', 'Critically Damped (\zeta = 1)', 'Overdamped (\zeta > 1)'};
colors   = {[0.25 0.45 0.90], [0.18 0.65 0.18], [0.85 0.20 0.20]};

figure('Position', [100 100 1600 400]);

for k = 1:3
    R_rlc = R_values(k);
    zeta  = (R_rlc/2) * sqrt(C_rlc/L_rlc);

    % ODE system
    rlc_ode = @(t, y) [y(2)/C_rlc; -y(1)/L_rlc - (R_rlc/L_rlc)*y(2)];

    % Solve
    [t_num_rlc, sol_rlc] = ode45(rlc_ode, t_span_rlc, [Vc0_rlc; 0]);
    Vc_rlc = sol_rlc(:,1);
    I_rlc  = sol_rlc(:,2);

    leg = sprintf('%s | \\zeta=%.2f', labels{k}, zeta);

    subplot(1,3,1)
    plot(t_num_rlc*1e3, Vc_rlc, 'Color', colors{k}, 'LineWidth', 2, 'DisplayName', leg); hold on;

    subplot(1,3,2)
    plot(t_num_rlc*1e3, I_rlc*1e3, 'Color', colors{k}, 'LineWidth', 2, 'DisplayName', leg); hold on;

    subplot(1,3,3)
    plot(Vc_rlc, I_rlc*1e3, 'Color', colors{k}, 'LineWidth', 2, 'DisplayName', leg); hold on;

    fprintf('%-30s R=%.1f Ohm  zeta=%.3f\n', labels{k}, R_rlc, zeta);
end

subplot(1,3,1)
yline(0, 'k:', 'LineWidth', 1); grid on;
xlabel('Time (ms)'); ylabel('Voltage (V)');
title('RLC — Capacitor Voltage'); legend('FontSize', 8);

subplot(1,3,2)
yline(0, 'k:', 'LineWidth', 1); grid on;
xlabel('Time (ms)'); ylabel('Current (mA)');
title('RLC — Current'); legend('FontSize', 8);

subplot(1,3,3)
grid on;
xlabel('Voltage (V)'); ylabel('Current (mA)');
title('RLC — Phase Portraits'); legend('FontSize', 8);

sgtitle(sprintf('RLC Circuit | L=%.1fH  C=%.0f\\muF  R_{crit}=%.1f\\Omega', ...
    L_rlc, C_rlc*1e6, R_critical), 'FontSize', 13);

fprintf('\nCritical resistance: R_crit = %.2f Ohm\n', R_critical);
%ζωητσπ