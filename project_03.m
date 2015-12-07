%% ECE 4712 - Project 3
% Problem 1
% Consider a simple circuit consisting of a 110V (single phase) source 
% connected to a load through a transmission line of impedance 0.1+j1 ohms.
% There are 500 load impedances that are switched ON one at a time. Express
% the circuit in per unit using a base of 110V and 1 KVA. Calculate the 
% load voltage, load power and load reactive power as the loads are 
% switched ON. Plot the load voltage as a function of load (real) power.
%   Case I: Each load impedance is a standard light bulb rated 100W, 120V
%   Case II: Each load impedance is a fluorescent tube light rated 100W, 
%           120V. Assume power factor of 0.8 lagging.
clear; clc;
% base units
V_base = 110; S_base = 1000;
Z_base = V_base^2 / S_base;

% single phase source
Vs = 110; %Volts 
Vs = Vs/V_base; %pu

% impedance of transmission line
Z_line = 0.1 + 1j; % ohms
Z_line = Z_line/Z_base; %pu

% number of load impedences
num_load = 0:500;

%% Case I: standard light bulb rated 100W, 120V

P_bulb = 100; % Watts
V_bulb = 120; % Volts

P_bulb = P_bulb/S_base; % pu
V_bulb = V_bulb/V_base; % pu 


Z_bulb = V_bulb^2 / P_bulb; 
Z_load = Z_bulb./num_load;
Z_total = Z_load + Z_line; 

% Current from source
I = Vs./Z_total;

% Voltage load
V_load = I.*Z_load

% Real and Reactive power of load
P_load = real(V_load.*conj(I))
Q_load = imag(V_load.*conj(I))

% plot load voltage as a function of real load power
figure(1)
plot(P_load,V_load)
xlabel('Real Power of load (pu)')
ylabel('Load Voltage (pu)')
title('Case I: PV Curve')
grid on;
parkit(0); % function to dock all figures --written by Dr. Obeid

%% Case II: fluorescent tube light rated 100W, 120V, and pf 0.8 lagging
P_tube = 100; % Watts
V_tube = 120; % Volts

% convert to power and voltage of tube light to pu
P_tube = P_tube/S_base;
V_tube = V_tube/V_base;

% power factor of the fluorescent tube light
pf_tube = 0.8;
angle = acos(0.8)
Z_tube = V_tube^2 / P_tube;
Z_tube = abs(Z_tube)*cos(angle) + 1i*abs(Z_tube)*sin(angle)
Zl_tube = Z_tube./num_load

Zt_tube = Z_line + Zl_tube

I2 = Vs./Zt_tube
Vload_t = I2.*Zl_tube
Pload_t = real(Vload_t.*conj(I2))
Qload_t = imag(Vload_t.*conj(I2))
figure(2)
plot(Pload_t,Vload_t)
xlabel('Real Power of load (pu)')
ylabel('Load Voltage (pu)')
title('Case II: PV Curve')
grid on;
parkit(0);
%% Comparison
figure(3)
plot(P_load,V_load)
xlabel('Real Power of load (pu)')
ylabel('Load Voltage (pu)')
title('Comparison')
hold on;
grid on;
plot(Pload_t,Vload_t)
hold off
legend('Standard light bulb', 'Fluorescent light tube, pf of 0.8 lagging')
parkit(0)




