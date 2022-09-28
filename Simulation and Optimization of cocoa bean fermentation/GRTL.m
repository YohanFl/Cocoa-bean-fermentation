%Growth rate multiplier Thornton & Lessem.

function [gr] = GRTL(T,Tmin,Tmax,Topt,K1,K4)

if nargin < 5

  K1 = 0.1;
  K4 = 0.1;

endif

 K2 = 0.98;
 K3 = 0.98;

 g1 = 1/(Topt-Tmin)*log(K2*(1-K1)/K1/(1-K2));
 g2 = 1/(Tmax-Topt)*log(K3*(1-K4)/K4/(1-K3));

 Ka = K1*exp(g1*(T-Tmin))/(1+K1*(exp(g1*(T-Tmin))-1));
 Kb = K4*exp(g2*(Tmax-T))/(1+K4*(exp(g2*(Tmax-T))-1));

 gr = Ka*Kb;

endfunction
