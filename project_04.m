%% ECE 4712 - Project 4: Short Circuit Analysis
% Devin Jiang, John Snyder, Abdulmagid Dahbali
% Analyze the short circuit condition for each of the following three cases
%   (a) Single line to ground fault with fault impedance
%   (b) Double line fault (no ground) with impedance
%   (c) Double line to ground fault
%% Initialize parameters
clear; clc;
% generator (pu)
X_g1 = 0.12; X_g2 = 0.12; X_g0 = 0.06;

% transformer (pu)
X_t1 = 0.10; X_t2 = 0.10; X_t0 = 0.08;

% transmission line (pu)
X_l1 = 0.15; X_l2 = 0.10; X_l0 = 0.20;

% load
S = 1.0 + 1i*0.5;

% bus 1 voltage (pu)
V1 = 1.0;

% fault impedance (pu)
Zf = 0.05;

%% (a) Single line to ground fault with fault impedance
