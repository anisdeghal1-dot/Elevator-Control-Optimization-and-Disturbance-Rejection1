function [c, ceq] = nlcon_ascenseur_robuste(x)

    wn   = x(1);
    zeta = x(2);

    mc_nominal = 80;
    delta_mc   = [-0.20, -0.15, -0.10, -0.05, 0, +0.05, +0.10, +0.15, +0.20];
    n_sc       = length(delta_mc);

    
    R=5.65; Re=0.445; L=2.90e-3; Ka=2.39e-2; J=0.481;
    B1=980; K1=145000; Blin=300; n1=1430; n2=190; r=0.1;
    mc_nom = 80;

    A_nom = [-(Re+R)/L,     -n1/(Ka*L*n2), 0,        0;
              n2/(n1*J*Ka), -B1/J,          0,        0;
              0,             K1/r,          0,       -K1;
              0,             0,             1/mc_nom, -Blin/mc_nom];
    B_ = [1/L; 0; 0; 0];

    sigma = zeta * wn;
    wd    = wn * sqrt(max(0, 1 - zeta^2));

    
    v_OS   = zeros(1, n_sc);
    v_Tr   = zeros(1, n_sc);
    v_Mg   = zeros(1, n_sc);
    v_stab = zeros(1, n_sc);

    for k = 1:n_sc
        mc_k = mc_nominal * (1 + delta_mc(k));

        
        [~,~,~,OS,Tr,Mg] = systeme_ascenseur_robuste(wn, zeta, mc_k);
        v_OS(k) = OS - 5;      % <= 0 requis
        v_Tr(k) = Tr - 1.0;   % <= 0 requis  (meme spec que nlcon_ascenseur.m)
        v_Mg(k) = 6 - Mg;     % <= 0 requis

        
        A_pert = [-(Re+R)/L,     -n1/(Ka*L*n2), 0,       0;
                   n2/(n1*J*Ka), -B1/J,          0,       0;
                   0,             K1/r,          0,      -K1;
                   0,             0,             1/mc_k, -Blin/mc_k];
        try
            K_k    = place(A_nom, B_, [-sigma+1j*wd; -sigma-1j*wd; -6*wn; -7.2*wn]);
            lambda = eig(A_pert - B_*K_k);
            v_stab(k) = max(real(lambda)) + 0.5;  % même marge que nlcon_ascenseur.m
        catch
            v_stab(k) = 10;
        end
    end

    
    c = zeros(4, 1);
    c(1) = max(v_OS);    % pire dépassement
    c(2) = max(v_Tr);    % pire temps de réponse
    c(3) = max(v_Mg);    % pire marge de gain
    c(4) = max(v_stab);  % pire stabilité

    ceq = [];
end
