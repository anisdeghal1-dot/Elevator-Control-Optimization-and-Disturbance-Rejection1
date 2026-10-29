function fob = funobj_ascenseur_ga(x)
    % x(1) = wn, x(2) = zeta
    wn   = x(1);
    zeta = x(2);

    % Appel du simulateur
    [t, y, yref, OS, Tr, Mg] = systeme_ascenseur(wn, zeta);

    
    erreur = trapz(t, abs(yref - y));

    
    pen_Mg = max(0, 6 - Mg)^2;  % Seuil à 6dB selon cahier des charges
    pen_OS = max(0, OS - 5)^2;   % OS <= 5%
    pen_Tr = max(0, Tr - 1.0)^2; % Tr <= 1s

    %% Effort de commande
    pen_u = 0.01 * wn^2;

    %Fonction obj
    fob = erreur + 0.2 * pen_Tr + 0.1 * pen_OS + 0.1 * pen_Mg + pen_u;
end