close all
clc
clear all
R   = 5.65;     Re  = 0.445;    L   = 2.90e-3;
Ka  = 2.39e-2;  J   = 0.481;    B1  = 980;
K1  = 145000;   mc  = 80;       Blin= 300;
n1  = 1430;     n2  = 190;      r   = 0.1;
d    = 4120;      Vmax = 500;


A = [ -(Re+R)/L,        -n1/(Ka*L*n2),   0,      0;
       n2/(n1*J*Ka),    -B1/J,           0,      0;
       0,                K1/r,           0,     -K1;
       0,                0,              1/mc,  -Blin/mc ]

B = [ 1/L; 0; 0; 0 ]
C = [ 0, 0, 0, 1 ]
D=0

polesBO = eig(A);
fprintf('Poles en boucle ouverte :\n');
disp(polesBO)

Co = ctrb(A, B);
r_Co = rank(Co);
fprintf('Rang matrice commandabilite : %d / 4\n', r_Co);
if r_Co == 4
    fprintf('Systeme COMMANDABLE -> placement de poles possible\n');
end

wn   = 1;    % Pulsation naturelle [rad/s]
zeta = 0.8;    % Amortissement [-]

% Pôles dominants (complexes conjugués)
sigma = zeta * wn;             % = 2.4
wd    = wn * sqrt(1-zeta^2);   % = 1.8 rad/s

p1 = -sigma + 1j*wd;   % = -2.4 + 1.8j
p2 = -sigma - 1j*wd;   % = -2.4 - 1.8j

% Pôles auxiliaires (modes rapides, sans influence sur la dynamique dominante)
p3 = -6  * wn;   % = -18
p4 = -7.2* wn;   % = -21.6

poles_desires = [p1; p2; p3; p4]
K=place(A,B,poles_desires)
A_BF = A - B*K;
polesBF = eig(A_BF)

Kr = -1 / (C * inv(A - B*K) * B);

sysBF = ss(A - B*K, B*Kr, C, D);

figure
step(sysBF)
grid on
title('Reponse indicielle avec suivi de consigne')

sim('untitled2023a.slx')
t=ans.y.time;
y_data=ans.y.signals.values;
info= stepinfo(y_data,t)

