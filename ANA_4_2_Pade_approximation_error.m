clear all
close all

x = [-10:0.01:5];
f = exp(x);
R = (1+x/2+x.^2/10+x.^3/120)./(1-x/2+x.^2/10-x.^3/120); % Padé approximation


figure
plot(x,f,"LineWidth",0.8)
hold on
plot(x,R,"LineWidth",0.8)
legend("real function","approximated function")
axis ([-10,5,0,1])
grid minor

figure
plot(x,abs(f-R),"LineWidth",0.8)
legend("error of the approximation")
axis ([-10,5,0,1])
grid minor