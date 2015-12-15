%% ECE 4712 - Project 4: Short Circuit Analysis
% Devin Jiang, John Snyder, Abdulmagid Dahbali
% Analyze the short circuit condition for each of the following three cases
%   (a) Single line to ground fault with fault impedance
%   (b) Double line fault (no ground) with impedance
%   (c) Double line to ground fault
%% Initialize parameters
clear; clc;

% generator per units
gx1 = 0.12*1i; %pu
gx2 = 0.12*1i; %pu
gx0 = 0.06*i1; %pu 

% transformer per units
tx1 = 0.10*1i; %pu
tx2 = 0.10*1i; %pu
tx0 = 0.08*1i; %pu

% transmission line per units
trx1 = 0.15*1i; %pu
trx2 = 0.15*1i; %pu
trx0 = 0.20*1i; %pu

% load bus in rect
S = 1.0+(1i*0.5); %pu

% bus 1 voltage
bv1 = 1.0 + 0*1i; 

% fault line impedance 
zf = 0.05 ; %pu

z_1th = gx1;
z_2th = gx2;
z_0th = gx0;

%% (a) Single line to ground fault with fault impedance
fprintf('(a) Single line to ground fault with fault impedance')
% 1<120 and 1<240 as rect, a_1 = a, a_2 = a^2
a_1 = -0.5+ (-0.866*1i);
a_2 = -0.5+ (0.866*1i);

i_1 = bv1/(z_1th + z_2th + z_0th + 3*zf);
i_2 = i_1;
i_0 = i_1;

% %sequence voltage components
v_1 = bv1-(i_1*(gx1));
v_2 = -(i_2)*(gx1);
v_0 = -(i_0)*(gx0);

% 1<120 and 1<240 as rect, a_1 = a, a_2 = a^2
a_1 = -0.5 + (-0.866*1i);
a_2 = -0.5 + (0.866*1i);

% solve for i_a, i_b, and i_c
I_seq = [i_1; i_2; i_0];
I = [1 1 1; a_2 a_1 1; a_1 a_2 1]*I_seq;
I_a = I(1) 
I_b = I(2) 
I_c = I(3)

% %solving for v_a, v_b, and v_c
V_a = v_0+v_1+v_2
V_b = v_0+a_2*v_1+a_1*v_2
V_c = v_0+a_1*v_1+a_2*v_2

%% (b) Double Line Fault (no ground) with Impedance
fprintf('(b) Double Line Fault (no ground) with Impedance')
%currents
i_tot = bv1/(z_1th+z_2th);
i_1 = i_tot*((tx1+trx1+tx2+gx2)/((tx1+trx1+tx2+gx2)+gx1));
i_2 = -i_1;

%1<120 a_1 = a, a_2 = a^2
a_1=-0.5+ (0.87*1i);
a_2=-0.5+ (-0.87*1i);

%solve for i_a, i_b, i_c
i_a = 0
i_b = (a_2*i_1)+(a_1*i_2)
i_c = (a_1*i_1)+(a_2*i_2)

%solve for v_a, v_b, and v_c
v_a1 = bv1-(i_1*gx1);
v_a2 = -(i_2*gx1);
v_a = v_a1 + v_a2
v_b = (i_2*gx1)
v_c = -(i_2*gx1)


%% (c) Double line to ground fault
fprintf('(c) Double line to ground fault')
i_a = 0;

