function [c, ceq] = nlcon_ascenseur_robuste(x)
% NLCON_ASCENSEUR_ROBUSTE  Contraintes non-lineaires.


    wn   = x(1);
    zeta = x(2);

    R=5.65; Re=0.445; L=2.90e-3; Ka=2.39e-2; J=0.481;
    B1=980; K1=145000; Blin=300; n1=1430; n2=190; r=0.1; mc=80;

    [~, ~, ~, OS, Tr, Mg] = systeme_ascenseur_robuste(wn, zeta);

    c_OS = OS - 5;
    c_Tr = Tr - 3.0;
    c_Mg = 6  - Mg;

    A = [-(Re+R)/L,     -n1/(Ka*L*n2), 0,    0;
          n2/(n1*J*Ka), -B1/J,          0,    0;
          0,             K1/r,          0,   -K1;
          0,             0,             1/mc, -Blin/mc];
    B_ = [1/L; 0; 0; 0];

    sigma = zeta*wn;
    wd    = wn*sqrt(max(0,1-zeta^2));

    try
        K_     = place(A, B_, [-sigma+1j*wd; -sigma-1j*wd; -6*wn; -7.2*wn]);
        lambda = eig(A - B_*K_);
        c_stab = max(real(lambda)) + 0.1;
    catch
        c_stab = 10;
    end

    c   = [c_OS; c_Tr; c_Mg; c_stab];
    ceq = [];
end
