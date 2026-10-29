

close all; clc; clear all;

mc_nominal = 80;                          % kg  masse nominale
delta_pct  = [-20,-15,-10,-5,0,5,10,15,20]; % scenarios (%)
mc_scen    = mc_nominal * (1 + delta_pct/100);

x0 = [0.5; 0.6]; 
lb = [0.5; 0.5];
ub = [15;  0.99];

x_nr = [4.569; 0.687];   % wn, zeta issus de optimisation_ascenseur.m

fprintf('=========================================\n');
fprintf('  CORRECTEUR NON-ROBUSTE (patternsearch)\n');
fprintf('  wn=%.3f rad/s  zeta=%.3f\n', x_nr(1), x_nr(2));
fprintf('=========================================\n');
fprintf('%-6s  %8s  %8s  %8s\n', 'Δmc(%)', 'Tr(s)', 'OS(%)', 'Mg(dB)');
fprintf('%s\n', repmat('-',1,38));

Tr_nr = zeros(1,9); OS_nr = zeros(1,9); Mg_nr = zeros(1,9);
for k = 1:9
    [~,~,~,OS_k,Tr_k,Mg_k] = systeme_ascenseur_robuste(x_nr(1),x_nr(2),mc_scen(k));
    Tr_nr(k)=Tr_k; OS_nr(k)=OS_k; Mg_nr(k)=Mg_k;
    flag=''; if OS_k>5||Tr_k>1||Mg_k<6, flag=' �? VIOLATION'; end
    fprintf('  %+4d%%  %8.3f  %8.3f  %8.3f%s\n',delta_pct(k),Tr_k,OS_k,Mg_k,flag);
end

%  OPTIMISATION ROBUSTE� PATTERNSEARCH

fprintf('\n=========================================\n');
fprintf('  LANCEMENT PATTERNSEARCH ROBUSTE\n');
fprintf('  Critère minimax sur 9 scénarios de masse\n');
fprintf('=========================================\n');

options_rob = optimoptions('patternsearch', ...
    'Algorithm',              'classic',            ...
    'Display',                'iter',               ...
    'MaxFunctionEvaluations', 5000,                 ...  
    'MeshTolerance',          1e-6,                 ...
    'StepTolerance',          1e-6,                 ...
    'SearchFcn',              @MADSPositiveBasis2N, ...
    'PlotFcn', {@psplotbestf, @psplotmeshsize, ...
                @psplotfuncount, @psplotbestx});

[x_rob, fval_rob, exitflag_rob] = patternsearch( ...
    @funobj_ascenseur_robuste, x0, [], [], [], [], ...
    lb, ub, @nlcon_ascenseur_robuste, options_rob);

fprintf('\nCode de sortie : %d  |  fval = %.6f\n', exitflag_rob, fval_rob);
fprintf('Résultat robuste : wn = %.4f rad/s, zeta = %.4f\n', x_rob(1), x_rob(2));


fprintf('\n=========================================\n');
fprintf('  CORRECTEUR ROBUSTE (patternsearch minimax)\n');
fprintf('  wn=%.3f rad/s  zeta=%.3f\n', x_rob(1), x_rob(2));
fprintf('=========================================\n');
fprintf('%-6s  %8s  %8s  %8s\n', 'Δmc(%)', 'Tr(s)', 'OS(%)', 'Mg(dB)');
fprintf('%s\n', repmat('-',1,38));

Tr_rob = zeros(1,9); OS_rob = zeros(1,9); Mg_rob = zeros(1,9);
t_rob_all = cell(1,9); y_rob_all = cell(1,9);
t_nr_all  = cell(1,9); y_nr_all  = cell(1,9);

for k = 1:9
    [t_r,y_r,~,OS_k,Tr_k,Mg_k] = systeme_ascenseur_robuste(x_rob(1),x_rob(2),mc_scen(k));
    Tr_rob(k)=Tr_k; OS_rob(k)=OS_k; Mg_rob(k)=Mg_k;
    t_rob_all{k}=t_r; y_rob_all{k}=y_r;
    flag=''; if OS_k>5||Tr_k>1||Mg_k<6, flag=' �? VIOLATION'; end
    fprintf('  %+4d%%  %8.3f  %8.3f  %8.3f%s\n',delta_pct(k),Tr_k,OS_k,Mg_k,flag);

    [t_n,y_n] = deal([]); 
    [t_n,y_n,~] = systeme_ascenseur_robuste(x_nr(1),x_nr(2),mc_scen(k));
    t_nr_all{k}=t_n; y_nr_all{k}=y_n;
end


%  TABLEAU COMPARATIF 

idx0 = find(delta_pct==0);
fprintf('\n%s\n', repmat('=',1,68));
fprintf('  TABLEAU COMPARATIF SYNTHÉTIQUE\n');
fprintf('%s\n', repmat('=',1,68));
fprintf('%-28s %15s %15s\n','Critère','Non-robuste','Robuste');
fprintf('%s\n', repmat('-',1,68));
fprintf('%-28s %15.4f %15.4f\n','wn (rad/s)',      x_nr(1),  x_rob(1));
fprintf('%-28s %15.4f %15.4f\n','zeta',            x_nr(2),  x_rob(2));
fprintf('%s\n', repmat('-',1,68));
fprintf('  [Cas nominal — mc = 80 kg]\n');
fprintf('%-28s %15.3f %15.3f\n','Tr nominal (s)',  Tr_nr(idx0), Tr_rob(idx0));
fprintf('%-28s %15.3f %15.3f\n','OS nominal (%%)',  OS_nr(idx0), OS_rob(idx0));
fprintf('%-28s %15.3f %15.3f\n','Mg nominal (dB)', Mg_nr(idx0), Mg_rob(idx0));
fprintf('%s\n', repmat('-',1,68));
fprintf('  [Pire cas — sur ±20%%]\n');
fprintf('%-28s %15.3f %15.3f\n','Tr pire cas (s)',  max(Tr_nr), max(Tr_rob));
fprintf('%-28s %15.3f %15.3f\n','OS pire cas (%%)',  max(OS_nr), max(OS_rob));
fprintf('%-28s %15.3f %15.3f\n','Mg pire cas (dB)', min(Mg_nr), min(Mg_rob));
fprintf('%s\n', repmat('=',1,68));


%  GAINS DU CORRECTEUR ROBUSTE

R_=5.65; Re=0.445; L=2.90e-3; Ka=2.39e-2; J=0.481;
B1=980; K1=145000; Blin=300; n1=1430; n2=190; r=0.1; mc=80;
A=[-(Re+R_)/L,-n1/(Ka*L*n2),0,0; n2/(n1*J*Ka),-B1/J,0,0;
   0,K1/r,0,-K1; 0,0,1/mc,-Blin/mc];
B_=[1/L;0;0;0]; C_=[0,0,0,1];
sig=x_rob(2)*x_rob(1); wd_=x_rob(1)*sqrt(1-x_rob(2)^2);
K_rob  = place(A,B_,[-sig+1j*wd_;-sig-1j*wd_;-6*x_rob(1);-7.2*x_rob(1)]);
Kr_rob = -1/(C_*((A-B_*K_rob)\B_));
fprintf('\nK  = [%.6f  %.6f  %.6f  %.6f]\n', K_rob);
fprintf('Kr = %.6f\n', Kr_rob);


%  FIGURE 1 : Reponses indicielles toutes variations de masse

colors = cool(9);

figure('Position',[50 50 1100 500],'Name','Robustesse — Réponses indicielles');

subplot(1,2,1);
title('Correcteur NON-ROBUSTE (patternsearch nominal)','FontSize',11);
hold on; grid on;
for k = 1:9
    lw = 1.2; ls = '-';
    if delta_pct(k)==0, lw=2.5; end
    h=plot(t_nr_all{k}, y_nr_all{k}, ls, 'Color', colors(k,:), 'LineWidth', lw);
    if delta_pct(k)==0, set(h,'DisplayName','Nominal (80 kg)','Color','k');
    else, set(h,'DisplayName',sprintf('\\Deltamc=%+d%%',delta_pct(k))); end
end
yline(1,'r--','LineWidth',1.5,'DisplayName','Consigne');
yline(1.05,'r:','LineWidth',0.8,'HandleVisibility','off');
yline(0.95,'r:','LineWidth',0.8,'HandleVisibility','off');
xlabel('Temps (s)'); ylabel('Vitesse cabine (m/s)');
legend('show','Location','southeast','FontSize',7);
ylim([-0.05 1.35]);

subplot(1,2,2);
title('Correcteur ROBUSTE — patternsearch minimax','FontSize',11);
hold on; grid on;
for k = 1:9
    lw = 1.2;
    if delta_pct(k)==0, lw=2.5; end
    h=plot(t_rob_all{k}, y_rob_all{k}, '-', 'Color', colors(k,:), 'LineWidth', lw);
    if delta_pct(k)==0, set(h,'DisplayName','Nominal (80 kg)','Color','k');
    else, set(h,'DisplayName',sprintf('\\Deltamc=%+d%%',delta_pct(k))); end
end
yline(1,'r--','LineWidth',1.5,'DisplayName','Consigne');
yline(1.05,'r:','LineWidth',0.8,'HandleVisibility','off');
yline(0.95,'r:','LineWidth',0.8,'HandleVisibility','off');
xlabel('Temps (s)'); ylabel('Vitesse cabine (m/s)');
legend('show','Location','southeast','FontSize',7);
ylim([-0.05 1.35]);

sgtitle('Robustesse aux perturbations endogènes — Variation de masse mc',...
    'FontSize',13,'FontWeight','bold');


%  FIGURE 2 : Performances vs variation de masse

figure('Position',[100 600 1050 320],'Name','Indicateurs de robustesse');

subplot(1,3,1);
plot(delta_pct, Tr_nr, 'b-o','LineWidth',1.8,'MarkerSize',6,'DisplayName','Non-robuste'); hold on;
plot(delta_pct, Tr_rob,'g-s','LineWidth',1.8,'MarkerSize',6,'DisplayName','Robuste');
yline(1.0,'r--','Spec: 1s','LineWidth',1.5,'LabelHorizontalAlignment','left');
xlabel('Variation de masse (%)'); ylabel('Tr (s)');
title('Temps de réponse'); legend('show','Location','best'); grid on;

subplot(1,3,2);
plot(delta_pct, OS_nr, 'b-o','LineWidth',1.8,'MarkerSize',6,'DisplayName','Non-robuste'); hold on;
plot(delta_pct, OS_rob,'g-s','LineWidth',1.8,'MarkerSize',6,'DisplayName','Robuste');
yline(5.0,'r--','Spec: 5%','LineWidth',1.5,'LabelHorizontalAlignment','left');
xlabel('Variation de masse (%)'); ylabel('OS (%)');
title('Dépassement'); legend('show','Location','best'); grid on;

subplot(1,3,3);
plot(delta_pct, Mg_nr, 'b-o','LineWidth',1.8,'MarkerSize',6,'DisplayName','Non-robuste'); hold on;
plot(delta_pct, Mg_rob,'g-s','LineWidth',1.8,'MarkerSize',6,'DisplayName','Robuste');
yline(6.0,'r--','Spec: 6dB','LineWidth',1.5,'LabelHorizontalAlignment','left');
xlabel('Variation de masse (%)'); ylabel('Mg (dB)');
title('Marge de gain'); legend('show','Location','best'); grid on;

sgtitle('Indicateurs de performance vs variation de masse mc',...
    'FontSize',12,'FontWeight','bold');
