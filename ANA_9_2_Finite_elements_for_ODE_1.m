clear all
close all

x = [0,0.4,0.7,0.9,1.1,1.3,1.6,2];
u = [0, -0.00240701581, 0.04334356, 0.137134002, 0.323223755, 0.641900644, 1.47588789, 3.5];

y=6*cos(x)+3*(x.^2-2);

figure
plot(x,y,x,u)
legend("exact solution","approximate solution")
grid minor