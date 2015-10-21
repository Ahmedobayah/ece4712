% ece 4712
% project 1- load flow analysis

clear; clc;
% initialize
j = sqrt(-1);
num_bus = 4;
XL_km=0.376; % ohm/km at 60 Hz
RL_km= 0.037;
B_km=4.5e-6; % B mho/km
kV_LL = 345;
MVA3Ph=100;
Z_base = kV_LL^2/MVA3Ph;
B_base = 1/Z_base;
Y = zeros(4,4);
G = zeros(4,4);
B = zeros(4,4);

% transmission lines in km
l12 = 100;l13 = 200;l14 = 150;l23 = 120;l34 = 150;
l = [l12 l13 l14 l23 l34];

% impedance
Z12 = l12*RL_km + j*l12*XL_km; %?
Z13 = l13*RL_km + j*l13*XL_km;
Z14 = l14*RL_km + j*l14*XL_km;
Z23 = l23*RL_km + j*l23*XL_km;
Z34 = l34*RL_km + j*l34*XL_km;
Z_ohm = RL_km.*l + j*XL_km.*l;
Z(1,2) = Z12/Z_base; Z(2,1) = Z(1,2);
Z(1,3) = Z13/Z_base; Z(3,1) = Z(1,3);
Z(1,4) = Z14/Z_base; Z(4,1) = Z(1,4);
Z(2,3) = Z23/Z_base; Z(3,2) = Z(2,3);
Z(3,4) = Z34/Z_base; Z(4,3) = Z(3,4);
Z

% shunt susceptance B
B12 = B_km*l12/B_base;
B13 = B_km*l13/B_base;
B14 = B_km*l14/B_base;
B23 = B_km*l23/B_base;
B34 = B_km*l34/B_base;
% B = B_km.*l;
% B_pu = B/B_base;

% admittance matrix
for r=1:num_bus
    for c=1:num_bus
        Y(r,c) = -1/Z(r,c);
    end
end
Y(1,1) = 1/Z(1,2) + 1/Z(1,3) + Z(1,4);
Y(2,2) = 1/Z(1,2) + Z(2,3);
Y(3,3) = 1/Z(2,3) + 1/Z(1,3) + 1/Z(3,4);
Y(4,4) = 1/Z(1,4) + 1/Z(3,4);
Y(4,2) = 0;
Y(2,4) = 0;
Y
%
for r = 1:num_bus
    for c = 1:num_bus
        G(r,c) = real(Y(r,c));
        B(r,c) = imag(Y(r,c));
    end
end
G
B

% Solution Parameters
tolerance= 1e-07; iter_max=10;

% Given Specifications
V1MAG = 1.0; theta1 = 0;
P2sp = -1.2; Q2sp = -0.8;
P3sp = 1.5; V3MAG=1.0;
P4sp = -2.0; Q4sp = -1.6;

% Solve for:
P1sp = 0; Q1sp = 0;
V2MAG = 1.0; theta2 = 0;
Q3sp = 0; theta3 = 0;
V4MAG = 1.0; theta4 =0;

% initialize
iter=0; convFlag=1;
d_theta2 = 0; d_theta3 = 0; d_theta4 = 0;
d_v2_mag = 0; d_v4_mag = 0;

% Jacobian matrix
% Start Iteration Process for N-R
while( convFlag==1 && iter < iter_max)
    iter = iter + 1
    theta2 = theta2 + d_theta2;
    theta3 = theta3 + d_theta3;
    theta4 = theta4 + d_theta4;
    V2MAG = V2MAG + d_v2_mag;
    V4MAG = V4MAG + d_v4_mag;
    
    % Creation of Jacobian J
    % J(1,1) = dP2/dtheta2; k=2, m=1,3,4
    J(1,1) = V2MAG*( V1MAG*(-G(2,1)*sin(theta2-theta1) + B(2,1)*cos(theta2-theta1)) + V3MAG*(-G(2,3)*sin(theta2-theta3) + B(2,3)*cos(theta2-theta3)) + V4MAG*(-G(2,4)*sin(theta2-theta4) + B(2,4)*cos(theta2-theta4)) );
    % J(1,2) = dP2/dtheta3
    J(1,2) = V2MAG*V3MAG*(G(2,3)*sin(theta2-theta3) - B(2,3)*cos(theta2-theta3));
    % J(1,3) = dP2/d_theta4
    J(1,3) = V2MAG*V4MAG*(G(2,4)*sin(theta2-theta4) - B(2,4)*cos(theta2-theta4));
    % J(1,4) = dP2/d_v2_mag
    J(1,4) = 2*G(2,2)*V2MAG + V1MAG*(G(2,1)*cos(theta2-theta1) + B(2,1)*sin(theta2-theta1)) + V3MAG*(G(2,3)*cos(theta2-theta3) + B(2,3)*sin(theta2-theta3)) + V4MAG*(G(2,4)*cos(theta2-theta4) + B(2,4)*sin(theta2-theta4)) ;
    % J(1,5) = dP2/d_v4_mag
    J(1,5) = V2MAG*( G(2,4)*cos(theta2-theta4) + B(2,4)*sin(theta2-theta4) );

    %pause()
    % J(2,1) = dP3/d_theta2
    J(2,1) = V3MAG*V2MAG*(G(3,2)*sin(theta3-theta2) - B(3,2)*cos(theta3-theta2));
    % J(2,2) = dP3/d_theta3
    J(2,2) = V3MAG*( V1MAG*(-G(3,1)*sin(theta3-theta1) + B(3,1)*cos(theta3-theta1)) + V2MAG*(-G(3,2)*sin(theta3-theta2) + B(2,3)*cos(theta3-theta2)) + V4MAG*(-G(3,4)*sin(theta3-theta4) + B(3,4)*cos(theta3-theta4)) );
    % J(2,3) = dP3/d_theta4
    J(2,3) = V3MAG*V4MAG*( G(3,4)*sin(theta3-theta4) - B(3,4)*cos(theta3-theta4));
    % J(2,4) = dP3/d_v2_mag
    J(2,4) = V3MAG*( G(3,2)*cos(theta3-theta2) + B(3,2)*sin(theta3-theta2) );
    % J(2,5) = dP3/d_v4_mag
    J(2,5) = V3MAG*( G(3,4)*cos(theta3-theta4) + B(3,4)*sin(theta3-theta4) );

    % J(3,1) = dP4/d_theta2
    J(3,1) = V4MAG*V2MAG*( G(4,2)*sin(theta4-theta2) - B(4,2)*cos(theta4-theta2));
    % J(3,2) = dP4/d_theta3
    J(3,2) = V4MAG*V3MAG*( G(4,3)*sin(theta4-theta3) - B(4,3)*cos(theta4-theta3));
    % J(3,3) = dP4/d_theta4
    J(3,3) = V4MAG*( V1MAG*(-G(4,1)*sin(theta4-theta1) + B(4,1)*cos(theta4-theta1)) + V2MAG*(-G(4,2)*sin(theta4-theta2) + B(4,2)*cos(theta4-theta2)) + V3MAG*(-G(4,3)*sin(theta4-theta3) + B(4,3)*cos(theta4-theta3)) );
    % J(3,4) = dP4/d_v2_mag
    J(3,4) = V4MAG*( G(4,2)*cos(theta4-theta2) + B(4,2)*sin(theta4-theta2) );
    % J(3,5) = dP4/d_v4_mag
    J(3,5) = 2*G(4,4)*V4MAG + V1MAG*(G(4,1)*cos(theta4-theta1) + B(4,1)*sin(theta4-theta1)) + V2MAG*(G(4,2)*cos(theta4-theta2) + B(4,2)*sin(theta4-theta2)) + V3MAG*(G(4,3)*cos(theta4-theta3) + B(4,3)*sin(theta4-theta3)) ;

    % J(4,1) = dQ2/d_theta2
    J(4,1) = V2MAG*( V1MAG*(G(2,1)*cos(theta2-theta1) + B(2,1)*sin(theta2-theta1)) + V3MAG*(G(2,3)*cos(theta2-theta3) + B(2,3)*sin(theta2-theta3)) + V4MAG*( G(2,4)*cos(theta2-theta4) + B(2,4)*sin(theta2-theta4)) );
    % J(4,2) = dQ2/d_theta3
    J(4,2) = V2MAG*V3MAG*( -G(2,3)*cos(theta2-theta3) - B(2,3)*sin(theta2-theta3) );
    % J(4,3) = dQ2/d_theta4
    J(4,3) = V2MAG*V4MAG*( -G(2,4)*cos(theta2-theta4) - B(2,4)*sin(theta2-theta4) );
    % J(4,4) = dQ2/d_v2_mag
    J(4,4) = -2*B(2,2)*V2MAG + V1MAG*( G(2,1)*sin(theta2-theta1) - B(2,1)*cos(theta2-theta1) ) + V3MAG*( G(2,3)*sin(theta2-theta3) - B(2,3)*cos(theta2-theta3) ) + V4MAG*( G(2,4)*sin(theta2-theta4) - B(2,4)*cos(theta2-theta4) );
    % J(4,5) = dQ2/d_v4_mag
    J(4,5) = V2MAG*( G(2,4)*sin(theta2-theta4) - B(2,4)*cos(theta2-theta4) );

    % J(5,1) = dQ4/d_theta2
    J(5,1) = V4MAG*V2MAG*( -G(4,2)*cos(theta4-theta2) - B(4,2)*sin(theta4-theta2) );
    % J(5,2) = dQ4/d_theta3
    J(5,2) = V4MAG*V3MAG*( -G(4,3)*cos(theta4-theta3) - B(4,3)*sin(theta4-theta3) );
    % J(5,3) = dQ4/d_theta4
    J(5,3) = V4MAG*( V1MAG*(G(4,1)*cos(theta4-theta1) + B(4,1)*sin(theta4-theta1)) + V2MAG*(G(4,2)*cos(theta4-theta2) + B(4,2)*sin(theta4-theta2)) + V3MAG*( G(4,3)*cos(theta4-theta3) + B(4,3)*sin(theta4-theta3)) );
    % J(5,4) = dQ4/d_v2_mag
    J(5,4) = V4MAG*( G(4,2)*sin(theta4-theta2) - B(4,2)*cos(theta4-theta2) );
    % J(5,5) = dQ4/d_v4_mag
    J(5,5) = -2*B(4,4)*V4MAG + V1MAG*( G(4,1)*sin(theta4-theta1) - B(4,1)*cos(theta4-theta1) ) + V2MAG*( G(4,2)*sin(theta4-theta2) - B(4,2)*cos(theta4-theta2) ) + V3MAG*( G(4,3)*sin(theta4-theta3) - B(4,3)*cos(theta4-theta3) );
    J
    
    % Bus Voltages
    V(1,1)=V1MAG*exp(j*theta1); 
    V(2,1)=V2MAG*exp(j*theta2);
    V(3,1)=V3MAG*exp(j*theta3);
    V(4,1)=V4MAG*exp(j*theta4);
    
    % Injected currents into Buses
    Iinj=Y*V;
    % P and Q Injected into Buses
    S(1,1)=V(1,1)*conj(Iinj(1)); S(2,1)=V(2,1)*conj(Iinj(2)); S(3,1)=V(3,1)*conj(Iinj(3)); S(4,1)=V(4,1)*conj(Iinj(4));
    
    % Mismatch at PQ and PV buses
    Mismatch(1,1)=P2sp-real(S(2,1));
    Mismatch(2,1)=P3sp-real(S(3,1));
    Mismatch(3,1)=P4sp-real(S(4,1));
    Mismatch(4,1)=Q2sp-imag(S(2,1));
    Mismatch(5,1)=Q4sp-imag(S(4,1));
    Mismatch
    
    % calculate new delta values for theta2, theta3, and MAG3
    d = inv(J)*Mismatch; 
    d_theta2 = d(1); 
    d_theta3 = d(2); 
    d_theta4 = d(3); 
    d_v2_mag = d(4); 
    d_v4_mag = d(5);
    
    if max(abs(Mismatch)) > tolerance,
        convFlag=1;
    else
        convFlag=0;
    end
end

J % Final Jacobian Matrix

ANG2DEG=theta2*180/pi, ANG3DEG=theta3*180/pi, ANG4DEG=theta4*180/pi, V2MAG, V4MAG,
% Calculate Power Flow on the Transmission Lines
disp('bus voltages');
fprintf('V1 = %f pu , theta1 = %f° \n', V1MAG,theta1*180/pi);
fprintf('V2 = %f pu , theta2 = %f° \n', V2MAG,theta2*180/pi);
fprintf('V3 = %f pu , theta3 = %f° \n', V3MAG,theta3*180/pi);
fprintf('V4 = %f pu , theta4 = %f° \n', V4MAG,theta4*180/pi);
% V1MAG,theta1*180/pi
% V2MAG,theta2*180/pi
% V3MAG,theta3*180/pi
% V4MAG,theta4*180/pi


P12=real(V(1,1)*conj((V(1,1)-V(2,1))/Z12)), Q12=imag(V(1,1)*conj((V(1,1)-V(2,1))/Z12)), % at Bus 1
P13=real(V(1,1)*conj((V(1,1)-V(3,1))/Z13)), Q13=imag(V(1,1)*conj((V(1,1)-V(3,1))/Z13)), % at Bus 1
P14=real(V(1,1)*conj((V(1,1)-V(4,1))/Z14)), Q14=imag(V(1,1)*conj((V(1,1)-V(4,1))/Z14)), % at Bus 1
P23=real(V(2,1)*conj((V(2,1)-V(3,1))/Z23)), Q23=imag(V(2,1)*conj((V(2,1)-V(3,1))/Z23)), % at Bus 2
P34=real(V(3,1)*conj((V(3,1)-V(4,1))/Z34)), Q34=imag(V(3,1)*conj((V(3,1)-V(4,1))/Z34)), % at Bus 3
fprintf('Power at each bus.\n');
fprintf('P12 = %f pu, Q12 = %f pu\n', P12, Q12);
fprintf('P13 = %f pu, Q13 = %f pu\n', P13, Q13);
fprintf('P14 = %f pu, Q14 = %f pu\n', P14, Q14);
fprintf('P23 = %f pu, Q23 = %f pu\n', P23, Q23);
fprintf('P34 = %f pu, Q34 = %f pu\n', P34, Q34);

S(1,1), S(2,1), S(3,1), S(4,1), 
fprintf('S1 = %f pu, Q34 = %f pu\n', P34, Q34);