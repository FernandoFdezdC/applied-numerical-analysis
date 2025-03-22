clear all
close all

x=-2:0.01:2;
y=cos(x).^10;
coeff=[1, -7.71771e-014, -4.61032, 3.65448e-013, 7.77931, -4.54332e-013, -5.66846, 1.91147e-013, 1.77595, -2.47997e-014, -0.194228];
y_l=coeff(1);
for i=2:11
    y_l=y_l+coeff(i)*x.^(i-1);
end

figure
plot(x,y,x,y_l)
legend("cos^{10}(x)","Lagrange interpolating polynomial of degree 10")
title("Lagrange interpolation polynomials");
xlabel("X");
ylabel("Y");
axis([-2 2 0 1.2])
grid minor
