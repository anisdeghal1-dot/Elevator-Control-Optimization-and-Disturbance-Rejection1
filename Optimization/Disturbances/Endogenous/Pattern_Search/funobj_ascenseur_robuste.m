function fob = funobj_ascenseur_robuste(x)


    wn   = x(1);
    zeta = x(2);

    mc_nominal = 80; % kg
    delta_mc   = [-0.20, -0.15, -0.10, -0.05, 0, +0.05, +0.10, +0.15, +0.20];
    n_sc       = length(delta_mc);
    scores     = zeros(1, n_sc);

    for k = 1:n_sc
        mc_k = mc_nominal * (1 + delta_mc(k));
        [t, y, yref, OS, Tr, Mg] = systeme_ascenseur_robuste(wn, zeta, mc_k);

        erreur  = trapz(t, abs(yref - y));
        pen_Mg  = max(0, 10 - Mg)^2;
        pen_OS  = max(0, OS - 5)^2;
        pen_Tr  = max(0, Tr - 1.0)^2;
        pen_u   = 0.01 * wn^2;

        scores(k) = erreur + 0.2*pen_Tr + 0.1*pen_OS + 0.1*pen_Mg + pen_u;
    end

   
    fob = 0.6 * max(scores) + 0.4 * mean(scores);
end
