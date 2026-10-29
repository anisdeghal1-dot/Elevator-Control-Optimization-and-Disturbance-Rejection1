
close all; clc; clear all;

mc_nominal = 80;  % kg (masse nominale)

% Variations de masse pour l'analyse de robustesse
delta_pct  = [-20, -15, -10, -5, 0, 5, 10, 15, 20];  % en %
delta_frac  = delta_pct / 100;
mc_scenarios = mc_nominal * (1 + delta_frac);

% Bornes de recherche
x0 = [0.5; 0.6];  % point initial
lb = [0.5; 0.5];
ub = [15;  0.99];


fprintf('========================================\n');
fprintf('  POINT INITIAL (wn=%.1f, zeta=%.1f)\n', x0(1), x0(2));
fprintf('========================================\n');
fprintf('%-6s  %8s  %8s  %8s\n', 'Δmc(%)', 'Tr(s)', 'OS(%)', 'Mg(dB)');
fprintf('%s\n', repmat('-', 1, 38));
for k = 1:length(delta_pct)
    mc_k = mc_scenarios(k);
    [~,~,~,OS_k,Tr_k,Mg_k] = systeme_ascenseur_robuste(x0(1), x0(2), mc_k);
    flag = '';
    if OS_k > 5 || Tr_k > 3 || Mg_k < 6
        flag = ' �? VIOLATION';
    end
    fprintf('  %+4d%%  %8.3f  %8.3f  %8.3f%s\n', delta_pct(k), Tr_k, OS_k, Mg_k, flag);
end


% OPTIMISATION ROBUSTE (GA )
fprintf('\n========================================\n');
fprintf('  LANCEMENT OPTIMISATION ROBUSTE (GA)\n');
fprintf('  Critère : minimax sur 9 scénarios\n');
fprintf('========================================\n');

options_rob = optimoptions('ga', ...
    'Display',         'iter', ...
    'PopulationSize',  50, ...       
    'MaxGenerations',  80, ...       
    'EliteCount',      5, ...       
    'CrossoverFraction', 0.8, ...
    'PlotFcn',         {@gaplotbestf, @gaplotdistance}, ...
    'UseParallel',     false);

[x_rob, fval_rob, exitflag_rob] = ga(...
    @funobj_ascenseur_robuste, 2, [], [], [], [], ...
    lb, ub, @nlcon_ascenseur_robuste, options_rob);

fprintf('\nRésultat GA robuste : wn = %.4f rad/s, zeta = %.4f\n', x_rob(1), x_rob(2));


% ANALYSE DE ROBUSTESSE — CORRECTEUR ROBUSTE

fprintf('\n========================================\n');
fprintf('  ANALYSE ROBUSTESSE — CORRECTEUR ROBUSTE\n');
fprintf('  wn=%.3f rad/s, zeta=%.3f\n', x_rob(1), x_rob(2));
fprintf('========================================\n');
fprintf('%-6s  %8s  %8s  %8s\n', 'Δmc(%)', 'Tr(s)', 'OS(%)', 'Mg(dB)');
fprintf('%s\n', repmat('-', 1, 38));

Tr_rob_all = zeros(1, length(delta_pct));
OS_rob_all = zeros(1, length(delta_pct));
Mg_rob_all = zeros(1, length(delta_pct));

for k = 1:length(delta_pct)
    mc_k = mc_scenarios(k);
    [~,~,~,OS_k,Tr_k,Mg_k] = systeme_ascenseur_robuste(x_rob(1), x_rob(2), mc_k);
    Tr_rob_all(k) = Tr_k;
    OS_rob_all(k) = OS_k;
    Mg_rob_all(k) = Mg_k;
    flag = '';
    if OS_k > 5 || Tr_k > 3 || Mg_k < 6
        flag = ' �? VIOLATION';
    end
    fprintf('  %+4d%%  %8.3f  %8.3f  %8.3f%s\n', delta_pct(k), Tr_k, OS_k, Mg_k, flag);
end



fprintf('\n%s\n', repmat('=', 1, 70));
fprintf('  TABLEAU COMPARATIF SYNTHÉTIQUE\n');
fprintf('%s\n', repmat('=', 1, 70));
fprintf('%-25s %15s %15s %10s\n', 'Critère', 'Non-robuste', 'Robuste', 'Spec.');
fprintf('%s\n', repmat('-', 1, 70));
fprintf('%-25s %15.4f %15.4f %10s\n', 'wn (rad/s)',        x_rob(1), '—');
fprintf('%-25s %15.4f %15.4f %10s\n', 'zeta',              x_rob(2), '[0.5,0.99]');
fprintf('%s\n', repmat('-', 1, 70));
fprintf('  [Cas nominal mc=80kg]\n');
idx_nom = find(delta_pct == 0);
fprintf('%-25s %15.3f %15.3f %10s\n', 'Tr nominal (s)', Tr_rob_all(idx_nom), '<= 3 s');
fprintf('%-25s %15.3f %15.3f %10s\n', 'OS nominal (%%)', OS_rob_all(idx_nom), '<= 5 %%');
fprintf('%-25s %15.3f %15.3f %10s\n', 'Mg nominal (dB)', Mg_rob_all(idx_nom), '>= 6 dB');
fprintf('%s\n', repmat('-', 1, 70));
fprintf('  [Pire cas sur ±20%%]\n');
fprintf('%-25s %15.3f %15.3f %10s\n', 'Tr pire cas (s)',   max(Tr_rob_all), '<= 3 s');
fprintf('%-25s %15.3f %15.3f %10s\n', 'OS pire cas (%%)',  max(OS_rob_all), '<= 5 %%');
fprintf('%-25s %15.3f %15.3f %10s\n', 'Mg pire cas (dB)',  min(Mg_rob_all), '>= 6 dB');
fprintf('%s\n', repmat('=', 1, 70));

%% =========================================================
%  SECTION 5 : GAINS DU CORRECTEUR ROBUSTE
%  =========================================================
fprintf('\n--- Gains du correcteur robuste (mc nominal) ---\n');
R_=5.65; Re=0.445; L=2.90e-3; Ka=2.39e-2; J=0.481;
B1=980; K1=145000; Blin=300; n1=1430; n2=190; r=0.1;
A_nom=[-(Re+R_)/L,-n1/(Ka*L*n2),0,0; n2/(n1*J*Ka),-B1/J,0,0;
   0,K1/r,0,-K1; 0,0,1/mc_nominal,-Blin/mc_nominal];
B_=[1/L;0;0;0]; C_=[0,0,0,1];
sig=x_rob(2)*x_rob(1); wd_=x_rob(1)*sqrt(1-x_rob(2)^2);
K_rob  = place(A_nom, B_, [-sig+1j*wd_; -sig-1j*wd_; -6*x_rob(1); -7.2*x_rob(1)]);
Kr_rob = -1/(C_*((A_nom-B_*K_rob)\B_));
fprintf('K  = [%.6f  %.6f  %.6f  %.6f]\n', K_rob);
fprintf('Kr = %.6f\n', Kr_rob);


% Figure 1 : Réponses indicielles pour toutes les variations de masse
figure('Position', [50 50 900 600], 'Name', 'Robustesse — Réponses indicielles');
colors_rob = autumn(length(delta_pct));
figure
title('Correcteur ROBUSTE (minimax)'); hold on; grid on;
for k = 1:length(delta_pct)
    mc_k = mc_scenarios(k);
    [t_k, y_k, yref_k] = systeme_ascenseur_robuste(x_rob(1), x_rob(2), mc_k);
    if delta_pct(k) == 0
        plot(t_k, y_k, 'k-', 'LineWidth', 2.5, 'DisplayName', 'Nominal');
    else
        plot(t_k, y_k, '-', 'Color', colors_rob(k,:), 'LineWidth', 1.2, ...
            'DisplayName', sprintf('\\Deltamc=%+d%%', delta_pct(k)));
    end
end
plot(t_k, yref_k, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Consigne');
xlabel('Temps (s)'); ylabel('Vitesse cabine (m/s)');
legend('show', 'Location', 'southeast', 'FontSize', 8);
yline(1.05, 'r:', 'LineWidth', 0.8);
yline(0.95, 'r:', 'LineWidth', 0.8);

sgtitle('Robustesse aux perturbations endogènes (variation de masse mc)', 'FontSize', 13);

%Figure 2 : Indicateurs de performance vs variation de masse
figure('Position', [100 100 1000 350], 'Name', 'Indicateurs de robustesse');


subplot(1,2,1);
plot(delta_pct, OS_rob_all, 'g-s', 'LineWidth', 1.5, 'MarkerSize', 6, 'DisplayName', 'Robuste');
yline(5.0, 'r--', 'Spec: 5%', 'LineWidth', 1.5);
xlabel('Variation de masse (%)'); ylabel('Dépassement (%)');
title('Dépassement OS'); legend('show'); grid on;

subplot(1,2,2);
plot(delta_pct, Mg_rob_all, 'g-s', 'LineWidth', 1.5, 'MarkerSize', 6, 'DisplayName', 'Robuste');
yline(6.0, 'r--', 'Spec: 6dB', 'LineWidth', 1.5);
xlabel('Variation de masse (%)'); ylabel('Marge de gain (dB)');
title('Marge de gain M_g'); legend('show'); grid on;

sgtitle('Indicateurs de performance vs variation de masse — Comparaison robuste/non-robuste', 'FontSize', 12);
