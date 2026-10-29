function [c, ceq] = nlcon_ascenseur(x)
%Fonctions des contraintes

wn   = x(1);
zeta = x(2);

[~, ~, ~, Overshoot, Tr95] = systeme_ascenseur(wn, zeta);

Overshoot = real(Overshoot);
Tr95      = real(Tr95);

% --- Contraintes d'inegalite c <= 0 ---
c(1,1) = Overshoot - 5;   % Depassement <= 5%
c(2,1) = 0.5 - Tr95;      % Tr95 >= 0.5 s
c(3,1) = Tr95 - 5;        % Tr95 <= 5 s
c(4,1) = 0.6 - zeta;      % zeta >= 0.6
c(5,1) = 0.5 - wn;        % wn >= 0.5

% Forcer tout en reel
c = real(c);

% --- Contraintes d'egalitees ---
ceq = [];

end
