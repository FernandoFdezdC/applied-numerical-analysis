clear all
close all

a=0;
b=2;

x=a:1e-3:b;
u=0.23245614*x-0.807017544*x.^2+0.782894737*x.^3;

y=6*cos(x)+3*(x.^2-2);

figure
plot(x,y,x,u)
legend("exact solution","approximate solution")
grid minor