function [t, y, yref, Overshoot, Tr95] = systeme_ascenseur(wn, zeta)

warning off

% Parametres physiques
R   = 5.65;     Re  = 0.445;    L   = 2.90e-3;
Ka  = 2.39e-2;  J   = 0.481;    B1  = 980;
K1  = 145000;   mc  = 80;       Blin= 300;
n1  = 1430;     n2  = 190;      r   = 0.1;

%Matrices d'etat
A = [ -(Re+R)/L,        -n1/(Ka*L*n2),   0,      0;
       n2/(n1*J*Ka),    -B1/J,           0,      0;
       0,                K1/r,           0,     -K1;
       0,                0,              1/mc,  -Blin/mc ];
B = [ 1/L; 0; 0; 0 ];
C = [ 0, 0, 0, 1 ];
D = 0;

%Poles desires
sigma = zeta * wn;
wd    = wn * sqrt(max(0, 1 - zeta^2));
p1 = -sigma + 1j*wd;
p2 = -sigma - 1j*wd;
p3 = -6   * wn;
p4 = -7.2 * wn;
poles_desires = [p1; p2; p3; p4];

%Placement de poles + simulation.
try
    K  = place(A, B, poles_desires);
    Kr = -1 / (C * inv(A - B*K) * B);

    sysBF = ss(A - B*K, B*Kr, C, D);

    % Verifier stabilitee
    if any(real(eig(A - B*K)) >= 0)
        error('Système instable');
    end

    opt = stepDataOptions('StepAmplitude', 1);
    [y, t] = step(sysBF, 10, opt);
    yref = ones(size(t));

    info      = stepinfo(sysBF, 'SettlingTimeThreshold', 0.05);
    Overshoot = info.Overshoot;
    Tr95      = info.SettlingTime;

    % Si valeurs invalides → pénalité
    if isnan(Overshoot) || isinf(Overshoot)
        Overshoot = 100;
    end
    if isnan(Tr95) || isinf(Tr95)
        Tr95 = 100;
    end

catch
   
    t         = (0:0.1:10)';
    y         = zeros(size(t));
    yref      = ones(size(t));
    Overshoot = 100;
    Tr95      = 100;
end

end
