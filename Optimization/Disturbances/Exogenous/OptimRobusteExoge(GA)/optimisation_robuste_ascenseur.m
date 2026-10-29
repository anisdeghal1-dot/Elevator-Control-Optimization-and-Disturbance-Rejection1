%OPTIMISATION_ROBUSTE_ASCENSEUR
%  Optimisation du correcteur par retour d'etat avec perturbation exogene.
%  Perturbation : petit echelon +0.5V en entree a t=5s.

close all; clc; clear all;


x0 = [10; 0.7];
lb = [5;  0.5];
ub = [50; 0.99];


% POINTS INITIAUX

fprintf('========================================\n');
fprintf('  POINT INITIAL (wn=%.1f, zeta=%.2f)\n', x0(1), x0(2));
fprintf('========================================\n');
[~,~,~,OS0,Tr0,Mg0] = systeme_ascenseur_robuste(x0(1), x0(2));
flag0 = '';
if OS0>5 || Tr0>3 || Mg0<6, flag0=' <- VIOLATION'; end
fprintf('Tr=%.3f s | OS=%.3f%% | Mg=%.3f dB%s\n', Tr0, OS0, Mg0, flag0);


% OPTIMISATION GA

fprintf('\n========================================\n');
fprintf('  OPTIMISATION GA\n');
fprintf('  Perturbation exogene : +0.5V en entree a t=5s\n');
fprintf('========================================\n');

options_ga = optimoptions('ga', ...
    'Display',           'iter', ...
    'PopulationSize',    60, ...
    'MaxGenerations',    100, ...
    'EliteCount',        6, ...
    'CrossoverFraction', 0.8, ...
    'FunctionTolerance', 1e-6, ...
    'PlotFcn',           {@gaplotbestf, @gaplotdistance}, ...
    'UseParallel',       false);

[x_opt, fval_opt, exitflag_opt] = ga(...
    @funobj_ascenseur_robuste, 2, [], [], [], [], ...
    lb, ub, @nlcon_ascenseur_robuste, options_ga);

fprintf('\nResultat GA : wn=%.4f rad/s  zeta=%.4f\n', x_opt(1), x_opt(2));


%Simulation des  RESULTATS

[t_opt, y_opt, yref_opt, OS_opt, Tr_opt, Mg_opt] = systeme_ascenseur_robuste(x_opt(1), x_opt(2));
[t_ini, y_ini, ~,         OS_ini, Tr_ini, Mg_ini] = systeme_ascenseur_robuste(x0(1),    x0(2));

%Tableau comparatif
fprintf('\n%s\n', repmat('=',1,65));
fprintf('%-28s %15s %15s\n','Critere','Initial','Optimise GA');
fprintf('%s\n', repmat('-',1,65));
fprintf('%-28s %15.4f %15.4f\n','wn (rad/s)',  x0(1),   x_opt(1));
fprintf('%-28s %15.4f %15.4f\n','zeta',        x0(2),   x_opt(2));
fprintf('%s\n', repmat('-',1,65));
fprintf('%-28s %15.3f %15.3f\n','Tr (s)',      Tr_ini,  Tr_opt);
fprintf('%-28s %15.3f %15.3f\n','OS (%%)',     OS_ini,  OS_opt);
fprintf('%-28s %15.3f %15.3f\n','Mg (dB)',     Mg_ini,  Mg_opt);
fprintf('%s\n', repmat('=',1,65));
fprintf('Specifications : Tr<=3s | OS<=5%% | Mg>=6dB\n');

%Gains du correcteur optimise
fprintf('\n--- Gains du correcteur optimise ---\n');
R_=5.65; Re=0.445; L=2.90e-3; Ka=2.39e-2; J=0.481;
B1=980; K1=145000; Blin=300; n1=1430; n2=190; r=0.1; mc=80;
A_=[-(Re+R_)/L,-n1/(Ka*L*n2),0,0; n2/(n1*J*Ka),-B1/J,0,0;
    0,K1/r,0,-K1; 0,0,1/mc,-Blin/mc];
B__=[1/L;0;0;0]; C__=[0,0,0,1];
sig=x_opt(2)*x_opt(1); wd_=x_opt(1)*sqrt(1-x_opt(2)^2);
K_opt  = place(A_,B__,[-sig+1j*wd_;-sig-1j*wd_;-6*x_opt(1);-7.2*x_opt(1)]);
Kr_opt = -1/(C__*((A_-B__*K_opt)\B__));
fprintf('K  = [%.4f  %.4f  %.4f  %.4f]\n', K_opt);
fprintf('Kr = %.4f\n', Kr_opt);

%% =========================================================
%  SECTION 3 : FIGURES
%  =========================================================

% Figure 1Reponse complete
figure('Position',[50 50 850 480],'Name','Comparaison Initial vs Optimise');
hold on; grid on;
plot(t_ini, y_ini,  'b--','LineWidth',1.8,...
    'DisplayName',sprintf('Initial (wn=%.1f, \\zeta=%.2f)',x0(1),x0(2)));
plot(t_opt, y_opt,  'g-', 'LineWidth',2.2,...
    'DisplayName',sprintf('Optimise GA (wn=%.2f, \\zeta=%.3f)',x_opt(1),x_opt(2)));
plot(t_opt,yref_opt,'r--','LineWidth',1.2,'DisplayName','Consigne');
yline(1.05,'r:','LineWidth',0.8,'HandleVisibility','off');
yline(0.95,'r:','LineWidth',0.8,'HandleVisibility','off');
xline(5.0,'m:','LineWidth',1.2,'DisplayName','Perturbation +0.5V a t=5s');
xlabel('Temps (s)'); ylabel('Vitesse cabine (m/s)');
title('Reponse indicielle avec perturbation exogene en entree (+0.5V a t=5s)');
legend('show','Location','southeast','FontSize',9);

%Zoom sur le rejet de perturbation
figure('Position',[100 100 850 420],'Name','Zoom Rejet de Perturbation');
hold on; grid on;
idx_i = t_ini >= 4;
idx_o = t_opt >= 4;
plot(t_ini(idx_i),y_ini(idx_i),'b--','LineWidth',1.8,...
    'DisplayName',sprintf('Initial (wn=%.1f, \\zeta=%.2f)',x0(1),x0(2)));
plot(t_opt(idx_o),y_opt(idx_o),'g-','LineWidth',2.2,...
    'DisplayName',sprintf('Optimise GA (wn=%.2f, \\zeta=%.3f)',x_opt(1),x_opt(2)));
yline(1.0,'r--','LineWidth',1.2,'DisplayName','Consigne');
yline(1.05,'r:','LineWidth',0.8,'HandleVisibility','off');
yline(0.95,'r:','LineWidth',0.8,'HandleVisibility','off');
xline(5.0,'m:','LineWidth',1.2,'DisplayName','Perturbation +0.5V');
xlabel('Temps (s)'); ylabel('Vitesse cabine (m/s)');
title('Zoom rejet de perturbation exogene (t=4s a t=10s)');
legend('show','Location','northeast','FontSize',9);
