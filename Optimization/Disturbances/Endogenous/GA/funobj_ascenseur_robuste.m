function fob = funobj_ascenseur_robuste(x)

    wn   = x(1);
    zeta = x(2);

    %Masse nominale
    mc_nominal = 80; % kg

   
    delta_mc = [-0.20, -0.15, -0.10, -0.05, 0, +0.05, +0.10, +0.15, +0.20];
    n_scenarios = length(delta_mc);

    scores = zeros(1, n_scenarios);

    for k = 1:n_scenarios
        mc_k = mc_nominal * (1 + delta_mc(k));

        % Simulation avec la masse perturbee
        [t, y, yref, OS, Tr, Mg] = systeme_ascenseur_robuste(wn, zeta, mc_k);

        % Erreur de tracking
        erreur = trapz(t, abs(yref - y));

        
        pen_Mg = max(0, 6  - Mg)^2;   % Marge de gain >= 6 dB
        pen_OS = max(0, OS - 5)^2;    % Dépassement <= 5 %
        pen_Tr = max(0, Tr - 3.0)^2;  % Temps de réponse <= 3 s

        %Effort de commande
        pen_u = 0.01 * wn^2;

        
        scores(k) = erreur + 0.2*pen_Tr + 0.1*pen_OS + 0.1*pen_Mg + pen_u;
    end

    
  
    fob = 0.6 * max(scores) + 0.4 * mean(scores);
end
