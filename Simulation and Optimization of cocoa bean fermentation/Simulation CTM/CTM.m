%Function CTM gets the temperature T, the minimal growth temperature Tmin, the
%maximal growth temperature Tmax and the optimal growth temperature Topt as
%arguments and returns the growth rate multiplier tau.

function [tau] = CTM(T,Tmin,Tmax,Topt)

%Bernard & Reymond condition

if (Topt<=(Tmin+Tmax)/2)
  error('CMT fail. Topt too small')
endif


if (T>Tmin)&&(T<Tmax)
tau = (T-Tmax)*(T-Tmin)^2/((Topt-Tmin)*((Topt-Tmin)*(T-Topt)-(Topt-Tmax)*(Topt+Tmin-2*T)));
else tau = 0;
endif

endfunction
