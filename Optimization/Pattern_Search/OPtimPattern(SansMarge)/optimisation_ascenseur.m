


close all; clc; clear all;


% Chix de Point initial
x0 = [1.0; 0.8];   % wn=1 rad/s, zeta=0.8

fprintf('=== Point initial ===\n');
fprintf('  wn   = %.4f rad/s\n', x0(1));
fprintf('  zeta = %.4f\n',       x0(2));


%Test des points initial

fprintf('\n=== Test des contraintes au point initial ===\n');
[c0, ~] = nlcon_ascenseur(x0);
disp(c0)

f0 = funobj_ascenseur(x0);
fprintf('Valeur initiale de la fonction objectif : f0 = %.4f\n', f0);


%  Bornes sur [wn, zeta]

lb = [0.1; 0.5];  
ub = [10;  1.0];

fprintf('\n=== Bornes ===\n');
fprintf('  wn   ∈ [%.1f, %.1f] rad/s\n', lb(1), ub(1));
fprintf('  zeta ∈ [%.1f, %.1f]\n',        lb(2), ub(2));


%  Lancer l'algorithme d'Optimisation par patternsearch

fprintf('\n=== Optimisation : patternsearch (MADS) ===\n');

options = optimoptions('patternsearch', ...
    'Algorithm',              'classic', ...
    'Display',                'iter', ...
    'MaxFunctionEvaluations', 2000, ...
    'SearchFcn',              @MADSPositiveBasis2N, ...
    'PlotFcn', {@psplotbestf, @psplotmeshsize, @psplotfuncount, @psplotbestx});

[x_opt, fval, exitflag, output] = patternsearch( ...
    @funobj_ascenseur, x0, [], [], [], [], lb, ub, ...
    @nlcon_ascenseur, options);


%Affichage des resultats

fprintf('\n========================================\n');
fprintf('       RÉSULTATS DE L''OPTIMISATION      \n');
fprintf('========================================\n');
fprintf('  wn   optimal = %.4f rad/s\n', x_opt(1));
fprintf('  zeta optimal = %.4f\n',       x_opt(2));
fprintf('  Valeur objectif f* = %.4f\n', fval);
fprintf('  Exit flag = %d\n', exitflag);

% Verification des contraintes
[c_opt, ~] = nlcon_ascenseur(x_opt);
fprintf('\nContraintes au point optimal :\n');
labels = {'Dépassement <= 5%', 'Tr95 >= 0.5 s', 'Tr95 <= 5 s', ...
          'zeta >= 0.6', 'wn >= 0.5'};
for i = 1:length(c_opt)
    if c_opt(i) <= 0
        statut = 'SATISFAITE  ✓';
    else
        statut = 'VIOLÉE      ✗';
    end
    fprintf('  c(%d) = %8.4f  [%s]  (%s)\n', i, c_opt(i), statut, labels{i});
end

% Simulation et Performances avec parametres optimeaux
wn_opt   = x_opt(1);
zeta_opt = x_opt(2);

[t, y, yref, Overshoot, Tr95] = systeme_ascenseur(wn_opt, zeta_opt);

fprintf('\nPerformances obtenues :\n');
fprintf('  Dépassement          = %.2f %%\n', Overshoot);
fprintf('  Temps de réponse 5%% = %.2f s\n',  Tr95);

figure('Name', 'Réponse indicielle optimisée', 'NumberTitle', 'off');
plot(t, yref, 'r--', 'LineWidth', 1.5); hold on;
plot(t, y,    'b-',  'LineWidth', 2);
grid on;
xlabel('Temps (s)');
ylabel('Vitesse cabine (m/s)');
title('Réponse indicielle – Paramètres optimisés par patternsearch');
legend('Consigne (1 m/s)', 'Réponse optimisée', 'Location', 'southeast');
ylim([-0.1, 1.3]);


% Calcul des gains K et Kr finaux a�copier dans Simulink

R_   = 5.65;  Re = 0.445;  L = 2.90e-3;
Ka   = 2.39e-2; J = 0.481; B1 = 980;
K1   = 145000; mc = 80; Blin = 300;
n1   = 1430; n2 = 190; r = 0.1;

A = [ -(Re+R_)/L,       -n1/(Ka*L*n2),  0,     0;
       n2/(n1*J*Ka),    -B1/J,          0,     0;
       0,                K1/r,          0,    -K1;
       0,                0,             1/mc, -Blin/mc ];
B = [ 1/L; 0; 0; 0 ];
C = [ 0, 0, 0, 1 ];

sigma    = zeta_opt * wn_opt;
wd       = wn_opt * sqrt(1 - zeta_opt^2);
poles_opt = [-sigma+1j*wd; -sigma-1j*wd; -6*wn_opt; -7.2*wn_opt];

K_opt  = place(A, B, poles_opt);
Kr_opt = -1 / (C * inv(A - B*K_opt) * B);

fprintf('\n========================================\n');
fprintf('  GAINS\n');
fprintf('========================================\n');
fprintf('  K  = [%.6f  %.6f  %.6f  %.6f]\n', K_opt);
fprintf('  Kr = %.6f\n', Kr_opt);
