function [c, ceq] = nlcon_ascenseur_ga(x)
    wn   = x(1);
    zeta = x(2);

    % Recuperation des perfs du systeme
    [~, ~, ~, OS, Tr, Mg] = systeme_ascenseur(wn, zeta);

   
    R=5.65; Re=0.445; L=2.90e-3; Ka=2.39e-2; J=0.481;
    B1=980; K1=145000; mc=80; Blin=300; n1=1430; n2=190; r=0.1;

    A=[-(Re+R)/L, -n1/(Ka*L*n2), 0, 0;
       n2/(n1*J*Ka), -B1/J, 0, 0;
       0, K1/r, 0, -K1;
       0, 0, 1/mc, -Blin/mc];
    B=[1/L; 0; 0; 0];

    
    sigma = zeta * wn;
    wd    = wn * sqrt(max(0, 1 - zeta^2));
    K     = place(A, B, [-sigma+1j*wd; -sigma-1j*wd; -6*wn; -7.2*wn]);
    lambda = eig(A - B*K);

    %Cintraintes
    c = zeros(4,1);
    c(1) = OS - 5;                   % Depassement max 5%
    c(2) = Tr - 3.0;                 % Temps de reponse max 3s[cite: 3]
    c(3) = 6 - Mg;                  % Marge de gain min 6dB
    c(4) = max(real(lambda)) + 0.1;  % Stabilite (marge de 0.1)

    ceq = [];
end