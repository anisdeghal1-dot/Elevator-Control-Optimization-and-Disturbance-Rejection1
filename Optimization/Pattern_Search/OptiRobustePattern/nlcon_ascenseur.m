function [c, ceq] = nlcon_ascenseur(x)

wn   = x(1);
zeta = x(2);

[~, ~, ~, OS, Tr, Mg] = systeme_ascenseur(wn, zeta);


R=5.65; Re=0.445; L=2.90e-3; Ka=2.39e-2; J=0.481;
B1=980; K1=145000; mc=80; Blin=300; n1=1430; n2=190; r=0.1;

A=[-(Re+R)/L,-n1/(Ka*L*n2),0,0;
   n2/(n1*J*Ka),-B1/J,0,0;
   0,K1/r,0,-K1;
   0,0,1/mc,-Blin/mc];
B=[1/L;0;0;0];

sigma = zeta * wn;
wd    = wn * sqrt(max(0,1-zeta^2));
K  = place(A, B, [-sigma+1j*wd; -sigma-1j*wd; -6*wn; -7.2*wn]);

Acl = A - B*K;
lambda = eig(Acl);

%Contraintes
c = zeros(4,1);

c(1) = OS - 5;          % OS <= 5%
c(2) = Tr - 1.0;        % Tr <= 1s
c(3) = 6 - Mg;         % Mg >= 10 dB (robuste réaliste)
c(4) = max(real(lambda)) + 0.5;  % stabilité forte

ceq = [];
end