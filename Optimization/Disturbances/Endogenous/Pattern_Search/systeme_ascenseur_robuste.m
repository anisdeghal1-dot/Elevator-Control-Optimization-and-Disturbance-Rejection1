function [t, y, yref, Overshoot, Tr95, Mg] = systeme_ascenseur_robuste(wn, zeta, mc_val)


warning off all

%Parametres physiques
R=5.65; Re=0.445; L=2.90e-3; Ka=2.39e-2; J=0.481;
B1=980; K1=145000; Blin=300; n1=1430; n2=190; r=0.1;


mc_nom = 80;
if nargin < 3
    mc_val = mc_nom;
end

% Matrice nominale
A_nom = [ -(Re+R)/L,     -n1/(Ka*L*n2), 0,       0;
           n2/(n1*J*Ka), -B1/J,          0,       0;
           0,             K1/r,          0,      -K1;
           0,             0,             1/mc_nom, -Blin/mc_nom ];

% Matrice perturbe
A_pert = [ -(Re+R)/L,     -n1/(Ka*L*n2), 0,        0;
            n2/(n1*J*Ka), -B1/J,          0,        0;
            0,             K1/r,          0,       -K1;
            0,             0,             1/mc_val, -Blin/mc_val ];

B = [1/L; 0; 0; 0];
C = [0, 0, 0, 1];
D = 0;

try
    % Placement de poles sur le modele NOMINAL
    sigma = zeta * wn;
    wd    = wn * sqrt(max(0, 1 - zeta^2));
    poles = [-sigma+1j*wd; -sigma-1j*wd; -6*wn; -7.2*wn];

    K  = place(A_nom, B, poles);

    
    Acl = A_pert - B*K;

    if any(real(eig(Acl)) >= 0)
        error('instable');
    end

    Kr = -1 / (C * (Acl \ B));

    % Reponse indicielle
    sysBF = ss(Acl, B*Kr, C, D);
    [y, t] = step(sysBF, 10);
    yref   = ones(size(t));

    info      = stepinfo(sysBF, 'SettlingTimeThreshold', 0.05);
    Overshoot = info.Overshoot;
    Tr95      = info.SettlingTime;
    if ~isfinite(Overshoot), Overshoot = 100; end
    if ~isfinite(Tr95),      Tr95      = 100; end

    % Marges de stabilite
    sys_ol = ss(Acl, B, C, 0);
    L_tf   = Kr * tf(sys_ol);
    S      = allmargin(L_tf);

    gm_all = S.GainMargin;
    gm_all = gm_all(isfinite(gm_all) & gm_all > 0);
    if isempty(gm_all)
        Mg = 20;
    else
        Mg = 20 * log10(min(gm_all));
    end

catch
    t         = (0:0.1:10)';
    y         = zeros(size(t));
    yref      = ones(size(t));
    Overshoot = 100;
    Tr95      = 100;
    Mg        = 0;
end
end
