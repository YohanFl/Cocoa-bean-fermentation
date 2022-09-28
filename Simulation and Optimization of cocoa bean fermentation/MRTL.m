%Mortality rate multiplier Thornton & Lessem.

function [mr] = MRTL(T,Tmin,Tmax,Topt,K1,K4)


if nargin < 5

 K1 = 0.1;
 K4 = 0.1;

endif

 K2 = 0.98;
 K3 = 0.98;

 g2 = 1/(Tmax-Topt)*log(K3*(1-K4)/K4/(1-K3));

 Kb = K4*exp(g2*(Tmax-T))/(1+K4*(exp(g2*(Tmax-T))-1));

 mr = 1 - Kb;

endfunction
