clear all
close all

a=0;
b=2;

x=a:1e-3:b;
u=0.630950235*x-1.11041618*x.^2+0.83497053*x.^3;
d2u = -1.11041618*2 + 0.83497053*3*2*x;
R = d2u + u - 3*x.^2;

y=6*cos(x)+3*(x.^2-2);

figure
plot(x,y,x,u)
legend("exact solution","approximate solution")
grid minor

figure
plot(x,R)
legend("Residue")
grid minor