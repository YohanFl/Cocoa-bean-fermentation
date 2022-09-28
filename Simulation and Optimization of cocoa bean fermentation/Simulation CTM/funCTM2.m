%Function funCMT2 (Cardinal Temperature Model). Only the growth rates are here
% following the CMT model. The mortality rates follow the Arrhenius model. The
%function gets x and theta as argument and returns dx/dt.

function [dxdt]=funCTM2(x,theta)

dxdt = zeros(9,1);

v1 = theta(12)*CTM(x(9),theta(35),theta(32),theta(38))*x(1)*x(6)/(theta(17)+x(1));
v2 = theta(13)*CTM(x(9),theta(35),theta(32),theta(38))*x(2)*x(6)/(theta(18)+x(2));
v3 = theta(14)*CTM(x(9),theta(36),theta(33),theta(39))*x(1)*x(7)/(theta(19)+x(1));
v4 = theta(15)*CTM(x(9),theta(37),theta(34),theta(40))*x(3)*x(8)/(theta(20)+x(3));
v5 = theta(16)*CTM(x(9),theta(37),theta(34),theta(40))*x(4)*x(8)/(theta(21)*x(8)+x(4));
v6 = theta(22)*exp(-theta(25)/x(9))*x(6)*x(3);
v7 = theta(23)*exp(-theta(25)/x(9))*x(7)*x(4);
v8 = theta(24)*exp(-theta(25)/x(9))*x(8)*x(5)^2;


dxdt(1) = - theta(1)*v1 - theta(2)*v3;
dxdt(2) = - theta(3)*v2;
dxdt(3) = theta(4)*v1 + theta(5)*v2 - theta(6)*v4;
dxdt(4) = theta(7)*v3 - theta(8)*v5;
dxdt(5) = theta(9)*v3 + theta(10)*v4 +theta(11)*v5;
dxdt(6) = v1 + v2 - v6;
dxdt(7) = v3 - v7;
dxdt(8) = v4 + v5 - v8;
dxdt(9) = theta(28)*(theta(1)*v1+theta(2)*v3) + theta(29)*theta(3)*v2 + theta(30)*theta(6)*v4 + theta(31)*theta(8)*v5 - theta(26)*(x(9) - theta(27));

endfunction
