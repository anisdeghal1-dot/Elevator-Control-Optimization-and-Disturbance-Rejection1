function [t, y, yref, Overshoot, Tr95, Mg] = systeme_ascenseur_robuste(wn, zeta, ~)
% SYSTEME_ASCENSEUR_ROBUSTE  Retour d'etat avec perturbation exogene en entree.
%   Perturbation : petit echelon +0.5V en entree a t=5s.
%   Simulation via lsim() — evite le probleme de gain Kr*step().


warning off all

% Parametres physiques nominaux
R=5.65; Re=0.445; L=2.90e-3; Ka=2.39e-2; J=0.481;
B1=980; K1=145000; Blin=300; n1=1430; n2=190; r=0.1;
mc=80;

A = [ -(Re+R)/L,     -n1/(Ka*L*n2), 0,    0;
       n2/(n1*J*Ka), -B1/J,          0,    0;
       0,             K1/r,          0,   -K1;
       0,             0,             1/mc, -Blin/mc ];
B = [1/L; 0; 0; 0];
C = [0, 0, 0, 1];
D = 0;

try
    sigma = zeta * wn;
    wd    = wn * sqrt(max(0, 1 - zeta^2));

    %Poles en BF : 2 complexes conjugues + 2 reels
    poles = [-sigma+1j*wd; -sigma-1j*wd; -6*wn; -7.2*wn];

    K   = place(A, B, poles);
    Acl = A - B*K;

    if any(real(eig(Acl)) >= 0)
        error('instable');
    end

    % Gain de pre-compensation (scalaire)
    % Kr tel que la sortie en regime permanent = 1 pour une consigne = 1
    Kr = -1 / (C * (Acl \ B));

    % Vecteur temps
    t_total   = 10;
    t_perturb = 5.0;
    d_amp     = 0.5;   % amplitude perturbation exogene (V)
    Ns        = 2000;
    t_vec     = linspace(0, t_total, Ns)';

 
    %  Simulation correcte via lsim avec 2 entrees separees :
    %    u_consigne(t) = Kr * 1(t)          (echelon de consigne)
    %    u_perturb(t)  = d_amp * 1(t-5s)   (perturbation exogene)
    %  Les deux passent par la meme matrice B (meme canal d'entree).
    %  On utilise la superposition : y = y_c + y_p

    % Reponse a la consigne : u = Kr (echelon unitaire normalise)
    u_c = Kr * ones(Ns, 1);
    y_consigne = lsim(ss(Acl, B, C, D), u_c, t_vec);

    % Reponse a la perturbation : u = d_amp * 1(t >= t_perturb)
    u_p = zeros(Ns, 1);
    u_p(t_vec >= t_perturb) = d_amp;
    y_pert = lsim(ss(Acl, B, C, D), u_p, t_vec);

    % Superposition
    y    = y_consigne + y_pert;
    t    = t_vec;
    yref = ones(size(t));

    
    idx1 = t < t_perturb;
    y1   = y_consigne(idx1);
    yss  = 1.0;

    pk        = max(y1);
    Overshoot = max(0, (pk - yss) / yss * 100);
    if ~isfinite(Overshoot), Overshoot = 100; end

    outside = find(y1 < 0.95*yss | y1 > 1.05*yss, 1, 'last');
    if isempty(outside)
        Tr95 = 0;
    else
        Tr95 = t(outside);
    end
    if ~isfinite(Tr95) || Tr95 <= 0, Tr95 = 0.01; end

    %% Marge de gain (boucle ouverte avec Kr)
    % Boucle ouverte : L(s) = Kr * C * (sI - Acl)^{-1} * B
    sys_ol = ss(A - B*K, B, C, 0);   
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
    t         = (0:0.005:10)';
    y         = zeros(size(t));
    yref      = ones(size(t));
    Overshoot = 100;
    Tr95      = 100;
    Mg        = 0;
end
end
