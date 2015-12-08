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
%% initialize given
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
V_load = I.*Z_load;

% Real and Reactive power of load
P_load = real(V_load.*conj(I));
Q_load = imag(V_load.*conj(I));

% plot load voltage as a function of real load power
figure(1)
plot(P_load,V_load)
xlabel('Real Power of load (pu)')
ylabel('Load Voltage (pu)')
title('Case I: PV Curve')
grid on;
parkit(0); % function to dock all figures --written by Dr. Obeid

% display calculations
V_load(end-1)


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
%% Comparison between Case I and Case II
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
%% 2 (a)
clear; clc;
kV_LL = 345; %line voltage in kV
S_base = 100; % MVAthree phase
Z_base = kV_LL^2/S_base;

V1 = ones(1,11).*0.6;
V2 = [0.7812 0.7733 0.7652 0.7568 0.748 0.7389 0.7294 0.7194 0.7089 0.6979 0.6863];
V3 = ones(1,11);
V4 = [0.8712 0.8712 0.8712 0.8712 0.8711 0.8711 0.871 0.871 0.8709 0.8708 0.8707];
S = [abs(0+0i) abs(40+30i) abs(80+60i) abs(120+90i) abs(160+120i) ...
  abs(200+150i) abs(240+180i) abs(280+210i) abs(320+240i) abs(360+270i) abs(400+300i)];

plot(S./S_base, V1)
hold on;
plot(S./S_base, V2)
plot(S./S_base, V3)
plot(S./S_base, V4)
hold off;
grid on;
ylim([0.4,1.05])
title('V(S)')
xlabel('S (pu)')
ylabel('V_{bus} (pu)')
legend('V_1','V_2','V_3','V_4')
parkit(0);

%% 2 (b)
I_base = (100e6)/(345e3);
I12 = [972 946 924 907 895 890 890 898 912 934 964];
I13 = [1069 1073 1079 1085 1091 1099 1107 1117 1127 1138 1151];
I14 = [105 107 108 110 111 113 115 117 119 121 123];
I23 = [972 1003 1036 1071 1109 1148 1189 1233 1280 1329 1382];
I34 = [536 534 533 532 531 529 528 527 526 525 524];
S = [abs(0+0i) abs(40+30i) abs(80+60i) abs(120+90i) abs(160+120i) ...
  abs(200+150i) abs(240+180i) abs(280+210i) abs(320+240i) abs(360+270i) abs(400+300i)];
plot(S./S_base, I12./I_base)
hold on;
plot(S./S_base, I13./I_base)
plot(S./S_base, I14./I_base)
plot(S./S_base, I23./I_base)
plot(S./S_base, I34./I_base)
hold off;
grid on;
% ylim([0.4,1.05])
title('I(S)')
xlabel('S (pu)')
ylabel('Line flow, I (pu)')
legend('I_{12}','I_{13}','I_{14}','I_{23}', 'I_{34}')
parkit(0);
% %%
% P1 = [103 144 184 225 266 308 350 392 434 477 519];
% P2 = [  0  40  80 120 160 200 240 280 320 360 400];
% P3 = [150 150 150 150 150 150 150 150 150 150 150];
% P4 = [200 200 200 200 200 200 200 200 200 200 200];





