clear all
close all

x = 2:0.5:4;
u = [0.13343914, 0.122804456, 0.0995944283, 0.0758214214, 0.0563783214];

y = exp(-x).*(x-1);

figure
plot(x,y,x,u)
legend("exact solution","approximate solution")
grid minor