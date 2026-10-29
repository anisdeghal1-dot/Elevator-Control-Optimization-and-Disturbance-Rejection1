function [c, ceq] = nlcon_ascenseur_robuste(x)

    wn   = x(1);
    zeta = x(2);

    % Masse nominale et scenarios de variation
    mc_nominal = 80; % kg
    delta_mc   = [-0.20, -0.15, -0.10, -0.05, 0, +0.05, +0.10, +0.15, +0.20];
    n_scenarios = length(delta_mc);

  
    R=5.65; Re=0.445; L=2.90e-3; Ka=2.39e-2; J=0.481;
    B1=980; K1=145000; Blin=300; n1=1430; n2=190; r=0.1;

    % Contraintes : pire cas parmi tous les scenarios
    c_OS  = zeros(1, n_scenarios);  
    c_Tr  = zeros(1, n_scenarios); 
    c_Mg  = zeros(1, n_scenarios);  
    c_stab= zeros(1, n_scenarios); 

    for k = 1:n_scenarios
        mc_k = mc_nominal * (1 + delta_mc(k));

        % Simulation avec la masse perturbe
        [~, ~, ~, OS, Tr, Mg] = systeme_ascenseur_robuste(wn, zeta, mc_k);

        % Contraintes de performance
        c_OS(k)   = OS - 5;           % Dépassement <= 5 %
        c_Tr(k)   = Tr - 3.0;         % Temps de réponse <= 3 s
        c_Mg(k)   = 6 - Mg;           % Marge de gain >= 6 dB

        % Calcul de la stabilite avec la masse perturbe
        A_k = [-(Re+R)/L,     -n1/(Ka*L*n2), 0,     0;
                n2/(n1*J*Ka), -B1/J,          0,     0;
                0,             K1/r,           0,    -K1;
                0,             0,              1/mc_k, -Blin/mc_k];
        B_k = [1/L; 0; 0; 0];

        sigma = zeta * wn;
        wd    = wn * sqrt(max(0, 1 - zeta^2));

        try
            K_k    = place(A_k, B_k, [-sigma+1j*wd; -sigma-1j*wd; -6*wn; -7.2*wn]);
            lambda = eig(A_k - B_k*K_k);
            c_stab(k) = max(real(lambda)) + 0.1;  % Marge de stabilite 0.1
        catch
            c_stab(k) = 10;  % si placement impossible
        end
    end


    c = zeros(4, 1);
    c(1) = max(c_OS);    % Le pire depassement doit rester <= 5 %
    c(2) = max(c_Tr);    % Le pire Tr doit rester <= 3 s
    c(3) = max(c_Mg);    % La pire marge de gain doit rester >= 6 dB
    c(4) = max(c_stab);  % Le systeme doit rester stable dans tous les cas

    ceq = [];
end
