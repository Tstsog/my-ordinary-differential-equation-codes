%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code solve 1D ordinary differential equation (ODE) using the Euler's method [1]. 
%
% 1D ODE:  dx/dt = t/x, x(0) = 1 and t = [0,5];
% Euler's method with finite difference scheme: x(i+1) = x(i) + dt * (t(i)/x(i))
%
% An exact solution: x(t) = (1 + t*t)^(1/2);
%
% Reference book: D. Bradie, A friendly introduction to numerical analysis 
%
% Written by Tsogbayar Tsednee (PhD), Institute of Physics and Technology, Mongolian Academy of Sciences 
%
% Contact: tsog215@gmail.com
% Date: December 25, 2025
%
%%%
function [] = euler_example_3
%
clc; clear euler_example_3
%
N = 10; % number of steps
tb = 0.; % initial value of x
te = 5.; % final value of x
dt = (te - tb)/N; % step size
%
t = zeros(1,N);
x = zeros(1,N);
%
t(1) = 0.;
x(1) = 1;
%
for i = 1:N
    t(i+1) = tb + dt * i;
    x(i+1) = x(i) + dt * (t(i)/x(i));
end
%
x_exact = (1. + t.*t).^(1/2);
[t', x', x_exact']
%
[   0         1.0000    1.0000
    0.5000    1.0000    1.1180
    1.0000    1.2500    1.4142
    1.5000    1.6500    1.8028
    2.0000    2.1045    2.2361
    2.5000    2.5797    2.6926
    3.0000    3.0643    3.1623
    3.5000    3.5538    3.6401
    4.0000    4.0462    4.1231
    4.5000    4.5405    4.6098
    5.0000    5.0360    5.0990];
%

%%%
figure(1)
hold on
plot(t, x, 'b-', LineWidth=1.5)
plot(t, x_exact, 'ro', LineWidth=1.5)
xlabel('$t$','Interpreter','latex') % ,'fontsize',16
ylabel('$x(t)$','Interpreter','latex') % , 'Rotation',0 ,'Rotation',1
%axis([0. 3. 0.000 2.00])
set(gca,'FontSize',16)
hold off
box on

%%%
return
end
