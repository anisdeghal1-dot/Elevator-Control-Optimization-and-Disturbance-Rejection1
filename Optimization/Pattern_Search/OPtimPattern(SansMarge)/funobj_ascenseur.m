function fob = funobj_ascenseur(x)
% Fonction objectif ITAE  minimise l'erreur de suivi ponderes par le temps.


wn   = x(1);
zeta = x(2);

[t, y, yref, ~, ~] = systeme_ascenseur(wn, zeta);

N = length(t);

% Critère ITAE : penalise d'autant plus l'erreur qu'elle dure longtemps
fob = real((1/N) * (t' * abs(yref - y)));

end
