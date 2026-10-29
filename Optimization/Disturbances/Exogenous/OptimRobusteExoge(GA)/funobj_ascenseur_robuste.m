function fob = funobj_ascenseur_robuste(x)
% FUNOBJ_ASCENSEUR_ROBUSTE  Critere pour optimisation avec perturbation exogene.
%   Minimise erreur de tracking + penalites sur OS, Tr, Mg.

    wn   = x(1);
    zeta = x(2);

    [t, y, yref, OS, Tr, Mg] = systeme_ascenseur_robuste(wn, zeta);

    % Erreur de tracking (IAE)
    erreur = trapz(t, abs(yref - y));

    % Penalites
    pen_Mg = max(0, 6   - Mg)^2;
    pen_OS = max(0, OS  - 5 )^2;
    pen_Tr = max(0, Tr  - 3.0)^2;
    pen_u  = 0.001 * wn^2;  

    fob = erreur + 0.2*pen_Tr + 0.1*pen_OS + 0.1*pen_Mg + pen_u;
end
