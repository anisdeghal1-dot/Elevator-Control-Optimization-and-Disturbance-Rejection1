function fob = funobj_ascenseur(x)

wn   = x(1);
zeta = x(2);

[t, y, yref, OS, Tr, Mg] = systeme_ascenseur(wn, zeta);

% Erreur de tracking
erreur = trapz(t, abs(yref - y));

% Penalitees
pen_Mg = max(0, 10 - Mg)^2;
pen_OS = max(0, OS - 5)^2;
pen_Tr = max(0, Tr - 1.0)^2;

%Effort de commande

pen_u = 0.01 * wn^2;

%Objectif global
fob = erreur ...
    + 0.2 * pen_Tr ...
    + 0.1 * pen_OS ...
    + 0.1 * pen_Mg ...
    + pen_u;

end