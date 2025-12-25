%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code solve 1D ordinary differential equation (ODE) using the Euler's method [1]. 
%
% 1D ODE:  dy/dx = x + y, y(0) = 2 and x = [0,1];
% Euler's method with finite difference scheme: y(i+1) = y(i) + dx * (x(i) + y(i))
%
% An exact solution: y(x) = 3*exp(x) - x - 1;
%
% Reference book: D. Bradie, A friendly introduction to numerical analysis 
%
% Written by Tsogbayar Tsednee (PhD), Institute of Physics and Technology, Mongolian Academy of Sciences 
%
% Contact: tsog215@gmail.com
% Date: December 25, 2025
%
%%%
function [] = euler_example_1
%
clc; clear euler_example_1
%
N = 5; % number of steps
xb = 0.; % initial value of x
xe = 1.; % final value of x
dx = (xe - xb)/N; % step size
y(1) = 2;
%
for i = 1:N
    x(i+1) = xb + dx*i;
    y(i+1) = y(i) + dx * (x(i) + y(i));
end
%
y_exact = (3*exp(x) - x - 1);
%
[x', y', y_exact']
%
[        0    2.0000    2.0000
    0.2000    2.4000    2.4642
    0.4000    2.9200    3.0755
    0.6000    3.5840    3.8664
    0.8000    4.4208    4.8766
    1.0000    5.4650    6.1548
    ];


%%%
%%%
figure(1)
hold on
plot(x, y, 'b-', 'LineWidth', 1.5)
plot(x, y_exact, 'ro',  'LineWidth', 1.5)
xlabel('$x$','Interpreter','latex') % ,'fontsize',16
ylabel('$y(x)$','Interpreter','latex') % , 'Rotation',0 ,'Rotation',1
%axis([0. 3. 0.000 2.00])
set(gca,'FontSize',16)
hold off
box on

%%%
return
end
