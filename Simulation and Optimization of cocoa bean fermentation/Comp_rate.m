%Script to compare growth and mortality rates

x=0:.1:100;
y1 = [];
y2 = [];
AAB1 = [];
LAB = [];
AAB2 = [];
kY = [];
kY2 = [];
kLAB = [];
kAAB = [];
k1 = [];
k2 = [];
k3 =[];
k4 =[];
k5 =[];
mu1 = [];
mu2 = [];
mu3 = [];
mu4 = [];
mu5 = [];
m1 = [];
m2 = [];
m3 = [];

% Growth rates CTM
for i=1:size(x,2)
y1 = [y1,0.46101*CTM(x(i),3.15,41.79,29.32)];
endfor

for i=1:size(x,2)
y2 = [y2,0.62543*CTM(x(i),3.15,41.79,29.32)];
endfor

for i=1:size(x,2)
LAB = [LAB,0.60705*CTM(x(i),15,46.5,37.5)];
endfor

for i=1:size(x,2)
AAB1 = [AAB1,0.57379*CTM(x(i),8,35,30.9)];
endfor

for i=1:size(x,2)
AAB2 = [AAB2,0.016147*CTM(x(i),8,35,30.9)];
endfor

%Kill rate CTM

for i=1:size(x,2)
m1 = [m1,0.046517*exp(-12.692/x(i))];
endfor

for i=1:size(x,2)
m2 = [m2,0.0076645*exp(-12.692/x(i))];
endfor

for i=1:size(x,2)
m3 = [m3,0.009695*exp(-12.692/x(i))];
endfor



%kill rates TL

for i=1:size(x,2)
kY = [kY,0.046517*MRTL(x(i),3.15,41.79,29.32)];
endfor

for i=1:size(x,2)
kY2 = [kY2,0.046517*MRTL(x(i),3.15,41.79,29.32)];
endfor

for i=1:size(x,2)
kAAB = [kAAB,0.009695*MRTL(x(i),8,35,30.9)];
endfor

for i=1:size(x,2)
kLAB = [kLAB,0.0076645*MRTL(x(i),15,46.5,37.5)];
endfor

%Growth rates Thornton & Lessem

for i=1:size(x,2)
mu1 = [mu1,(0.46101+0.46101*4/100)*GRTL(x(i),3.5,41.79,29.32,0.005,0.005)];
endfor

for i=1:size(x,2)
mu2 = [mu2,(0.62543+0.62543*4/100)*GRTL(x(i),3.5,41.79,28.32,0.005,0.005)];
endfor

for i=1:size(x,2)
mu3 = [mu3,(0.60705+0.60705*4/100)*GRTL(x(i),15,46.5,37.5,0.005,0.005)];
endfor

for i=1:size(x,2)
mu4 = [mu4,(0.57379+0.57379*4/100)*GRTL(x(i),8,43,35,0.005,0.005)];
endfor

for i=1:size(x,2)
mu5 = [mu5,(0.016147+0.016147*4/100)*GRTL(x(i),8,43,35,0.005,0.005)];
endfor


%Growth rates Arrhenius

for i=1:size(x,2)
k1 = [k1,0.46101*exp(-17.723/x(i))];
endfor
for i=1:size(x,2)
k2 = [k2,0.62543*exp(-17.723/x(i))];
endfor

for i=1:size(x,2)
k3 = [k3,0.60705*exp(-17.723/x(i))];
endfor

for i=1:size(x,2)
k4 = [k4,0.57379*exp(-17.723/x(i))];
endfor

for i=1:size(x,2)
k5 = [k5,0.016147*exp(-17.723/x(i))];
endfor

close all


figure
title('growth rate comparaison')
plot(x,y1,'b',x,y2,'r',x,LAB,'c',x,AAB1,'k',x,AAB2,'g',x,mu1,'--b',x,mu2,'--r',x,mu3,'--c',x,mu4,'--k',x,mu5,'--g')
legend('Y1_{CTM}','Y2_{CTM}','LAB_{CTM}','AAB1_{CTM}','AAB2_{CTM}','Y1_{TL}','Y2_{TL}','LAB_{TL}','AAB1_{TL}','AAB2_{TL}')

figure
title('mortality rate comparaison')
plot(x,m1,'b',x,m2,'r',x,m3,'k',x,kY,'--b',x,kLAB,'--r',x,kAAB,'--k')
legend('Y_{CTM}','LAB_{CTM}','AAB_{CTM}','Y_{TL}','LAB_{TL}','AAB_{TL}')


