
close all; clc; clear all;


x0 = [0.5; 0.6];
lb = [0.5; 0.5];
ub = [15;  0.99];

[t_av, y_av, yref_av, OS_av, Tr_av, Mg_av] = systeme_ascenseur(x0(1), x0(2));
fprintf('--- AVANT ---\n');
fprintf('wn=%.2f  zeta=%.2f  Tr=%.2fs  OS=%.2f%%  Mg=%.2fdB  Mp=%.2fdeg\n', ...
    x0(1), x0(2), Tr_av, OS_av, Mg_av)

[c0, ~] = nlcon_ascenseur_ga(x0);
fprintf('Contraintes initiales (c<=0 = OK) : [%s]\n\n', ...
    strjoin(arrayfun(@(v) sprintf('%.2f',v), c0', 'UniformOutput', false), ', '));


% Configuration de l'Algorithme Genetique

options = optimoptions('ga', ...
    'Display', 'iter', ...
    'PopulationSize', 40, ...     % Taille de population
    'MaxGenerations', 50, ...     
    'PlotFcn', {@gaplotbestf, @gaplotdistance}, ... 
    'UseParallel', false);        

fprintf('--- Lancement Optimisation GA ---\n');

[x_opt_ga, fval_ga, exitflag_ga] = ga(...
    @funobj_ascenseur_ga, 2, [], [], [], [], ...
    lb, ub, @nlcon_ascenseur_ga, options);

fprintf('\nRésultat GA : wn = %.2f, zeta = %.2f\n', x_opt_ga(1), x_opt_ga(2));

% APRES
[t_op, y_op, ~, OS_op, Tr_op, Mg_op] = systeme_ascenseur(x_opt_ga(1), x_opt_ga(2));
fprintf('\n--- APRES ---\n');
fprintf('wn=%.2f  zeta=%.2f  Tr=%.2fs  OS=%.2f%%  Mg=%.2fdB  Mp=%.2fdeg\n', ...
    x_opt_ga(1), x_opt_ga(2), Tr_op, OS_op, Mg_op);

%Tableau comparatif
fprintf('\n%-25s %10s %10s %15s\n', 'Critere', 'Avant', 'Apres', 'Specification');
fprintf('%s\n', repmat('-', 1, 63));
fprintf('%-25s %10.3f %10.3f %15s\n', 'wn (rad/s)',           x0(1),  x_opt_ga(1), '—');
fprintf('%-25s %10.3f %10.3f %15s\n', 'zeta',                 x0(2),  x_opt_ga(2), '[0.5, 0.99]');
fprintf('%-25s %10.3f %10.3f %15s\n', 'Temps de reponse (s)', Tr_av,  Tr_op,    '<= 3 s');
fprintf('%-25s %10.3f %10.3f %15s\n', 'Depassement (%)',      OS_av,  OS_op,    '<= 5 %');
fprintf('%-25s %10.3f %10.3f %15s\n', 'Marge de gain (dB)',   Mg_av,  Mg_op,    '>= 6 dB');


% Verification finale des contraintes
fprintf('\n--- Verification contraintes apres optimisation ---\n');
[c_fin, ~] = nlcon_ascenseur_ga(x_opt_ga);
labels = {'OS<=5%', 'Tr<=3s', 'Mg>=6dB', 'Mp>=45deg', 'Mg>=0', 'Mp>=0'};
for i = 1:length(c_fin)
    if c_fin(i) <= 0
        statut = 'OK';
    else
        statut = '** VIOLEE **';
    end
    fprintf('  c(%d) = %+8.3f   [%s]  %s\n', i, c_fin(i), labels{i}, statut);
end

%Gains du correcteur optimal
R_=5.65; Re=0.445; L=2.90e-3; Ka=2.39e-2; J=0.481;
B1=980; K1=145000; mc=80; Blin=300; n1=1430; n2=190; r=0.1;
A=[-(Re+R_)/L,-n1/(Ka*L*n2),0,0; n2/(n1*J*Ka),-B1/J,0,0;
   0,K1/r,0,-K1; 0,0,1/mc,-Blin/mc];
B=[1/L;0;0;0]; C=[0,0,0,1];
sig=x_opt_ga(2)*x_opt_ga(1); wd=x_opt_ga(1)*sqrt(1-x_opt_ga(2)^2);
K_opt  = place(A, B, [-sig+1j*wd; -sig-1j*wd; -6*x_opt_ga(1); -7.2*x_opt_ga(1)]);
Kr_opt = -1/(C*((A-B*K_opt)\B));
fprintf('\nK  = [%.6f  %.6f  %.6f  %.6f]\n', K_opt);
fprintf('Kr = %.6f\n', Kr_opt);

%Trace comparatif
figure('Position', [100 100 700 480]);
plot(t_av, yref_av, 'r--', 'LineWidth', 1.5); hold on;
plot(t_av, y_av,    'b-',  'LineWidth', 2);
plot(t_op, y_op,    'g-',  'LineWidth', 2);
xlabel('Temps (s)');
ylabel('Vitesse cabine (m/s)');
title('Reponse indicielle — Avant vs Apres optimisation robuste');
legend('Consigne', ...
    sprintf('Avant  (wn=%.1f, Mg=%.1fdB)', x0(1),    Mg_av), ...
    sprintf('Apres  (wn=%.2f, Mg=%.1fdB)', x_opt_ga(1), Mg_op), ...
    'Location', 'southeast');
grid on;
