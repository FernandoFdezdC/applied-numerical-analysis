clear all

t = [0:0.1:50];             % time vector

w = sqrt([1.25;4.41511;2.93412;0.150768]);          % vector of normal mode frequencies

V = [0.57735 -0.428525 -0.656539 0.228013
0.57735 0.656538 0.228013 0.428525
-9.69502e-009 -0.57735 0.57735 0.57735
-0.57735 0.228013 -0.428525 0.656539];      % eigenvector matrix

x0 = [0.1;0.2;0.3;0.4];         % Initial positions
v0 = [1;2;-1;-2];              % Initial velocities

u = inv(V)*x0;
v = inv(V)*v0;

phi = atan((w.*u)./v);          % initial phases vector

C = diag(u./sin(phi));          % matrix with proper scalar multiplication for each eigenvector
x = V*C*sin(w*t + phi);

figure
plot(t,x(1,:),t,x(2,:),t,x(3,:),t,x(4,:))
legend("x_1(t)","x_2(t)","x_3(t)","x_4(t)")
title("Position for each mass x_{i} as a function of time t");
xlabel("time (s)");
ylabel("x_{i} (m)");
grid minor
