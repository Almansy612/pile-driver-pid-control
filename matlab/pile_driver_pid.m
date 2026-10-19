%% Pile-driver rig: active vibration damping with a hydraulic actuator + PID
% Corrected version of the Control Engineering coursework model.
% Requires Control System Toolbox. Mirrors analysis/model.py (linear part).

clear; clc; close all;

%% 1. Rig / pile structure (single degree of freedom)
m    = 1000;                  % kg
k    = 2e6;                   % N/m
zeta = 0.05;                  % structural damping ratio
wn   = sqrt(k/m);             % 44.72 rad/s
fn   = wn/(2*pi);             % 7.12 Hz   (wn/(2*pi), not wn/pi)
c    = 2*zeta*sqrt(k*m);      % 4472 N*s/m

%% 2. Hydraulic actuator (valve-controlled double-acting cylinder)
D = 0.050; d = 0.020; L = 0.200;           % bore, rod, stroke (m)
Pmax = 200e5;                              % 200 bar
beta = 1.5e9;                              % oil bulk modulus, Pa (1.5 GPa)
Qmax = 20/60000;                           % 20 L/min -> m^3/s
A    = pi/4*(D^2 - d^2);                   % annulus area 1.649e-3 m^2
Fmax = Pmax*A;                             % 32.99 kN
V    = A*L/2;                              % chamber volume at mid-stroke
tau_a = V*Pmax/(beta*Qmax);                % 6.6 ms  (pressure build-up)
% (dP/dt = beta/V*(Qmax*u - (Qmax/Pmax)*P)  ->  first order, gain Pmax)

s     = tf('s');
G_act = Fmax/(tau_a*s + 1);                % force per unit valve command u
G_rig = 1/(m*s^2 + c*s + k);               % displacement per unit force

%% 3. PID on rig displacement (setpoint 0), filtered derivative
Kp = 20; Ki = 1200; Kd = 1.8; N = 600;
C  = pid(Kp, Ki, Kd, 1/N);                 % Kp + Ki/s + Kd*s/(s/N + 1)

Lop = C*G_act*G_rig;                       % open loop
T_d = feedback(G_rig, C*G_act);            % disturbance force -> x
T_u = -C*T_d;                              % disturbance force -> valve command

[Gm, Pm, ~, Wcp] = margin(Lop);
fprintf('PM = %.1f deg, GM = %.1f dB, crossover = %.1f Hz\n', Pm, 20*log10(Gm), Wcp/(2*pi));

%% 4. Test 1: ram-impact pulse (20 kN half-sine, 25 ms)
t  = (0:1e-4:2)';
Fd = 20e3*sin(pi*t/0.025).*(t < 0.025);
x_open = lsim(G_rig, Fd, t);
x_pid  = lsim(T_d,   Fd, t);
u_pid  = lsim(T_u,   Fd, t);
fprintf('Impact: peak %.2f -> %.2f mm, max |u| = %.2f\n', ...
        1e3*max(abs(x_open)), 1e3*max(abs(x_pid)), max(abs(u_pid)));

figure('Name','Impact');
subplot(2,1,1); plot(t, 1e3*x_open, t, 1e3*x_pid, 'LineWidth', 1.4); grid on
ylabel('displacement (mm)'); legend('Uncontrolled','PID'); title('Ram-impact pulse')
subplot(2,1,2); plot(t, Fmax*u_pid/1e3, 'LineWidth', 1.4); grid on
yline([Fmax -Fmax]/1e3, '--'); ylabel('actuator force (kN)'); xlabel('time (s)')

%% 5. Test 2: ram weight (9.81 kN) as a static load step
Fs = 9.81e3*ones(size(t));
xs_open = lsim(G_rig, Fs, t);
xs_pid  = lsim(T_d,   Fs, t);
fprintf('Static load: final %.2f mm -> %.3f mm\n', 1e3*xs_open(end), 1e3*xs_pid(end));

figure('Name','Static load');
plot(t, 1e3*xs_open, t, 1e3*xs_pid, 'LineWidth', 1.4); grid on
ylabel('displacement (mm)'); xlabel('time (s)'); legend('Uncontrolled','PID')
title('Ram weight transferred to the rig')

%% 6. Frequency domain
figure('Name','Margins'); margin(Lop); grid on
figure('Name','Disturbance rejection'); bodemag(G_rig*k, T_d*k); grid on
legend('Uncontrolled','PID'); title('|x/F_d|·k')
