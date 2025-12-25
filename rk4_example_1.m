%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code solve 1D ordinary differential equation (ODE) using the 4th-order Runge-Kutta method [1]. 
%
% 1D ODE:  dy/dx = x + y, y(0) = 2 and x = [0,1];
% Runge-Kutta method: y(i+1) = y(i) + (1/6) * (k1 + 2 * k2 + 2 * k3 + k4)
%                         k1 =  dx * f(x(i),          y(i))  
%                         k2 =  dx * f(x(i) + 0.5*dx, y(i) + 0.5*k1)  
%                         k3 =  dx * f(x(i) + 0.5*dx, y(i) + 0.5*k2)
%                         k4 =  dx * f(x(i) +     dx, y(i) +     k3)  
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
function [] = rk4_example_1
%
clc; clear rk4_example_1
%
N = 5; % number of steps
xb = 0.; % initial value of x
xe = 1.; % final value of x
dx = (xe - xb)/N; % step size
y(1) = 2;
%
for i = 1:N
    x(i+1) = xb + dx*i;
    %
    y_k1 = dx * rhs_function(x(i),          y(i));
    %
    y_k2 = dx * rhs_function(x(i) + 0.5*dx, y(i) + 0.5*y_k1);
    %
    y_k3 = dx * rhs_function(x(i) + 0.5*dx, y(i) + 0.5*y_k2);
    %
    y_k4 = dx * rhs_function(x(i) +     dx, y(i) +     y_k1);
    %
    y(i+1) = y(i) + (y_k1/6.) + (y_k2/3.) + (y_k3/3.) + (y_k4/6.);
end
%
y_exact = (3*exp(x) - x - 1);
%
[x', y', y_exact']
%
[        0    2.0000    2.0000
    0.2000    2.4620    2.4642
    0.4000    3.0701    3.0755
    0.6000    3.8565    3.8664
    0.8000    4.8605    4.8766
    1.0000    6.1303    6.1548];

%%%
%%%
figure(1)
hold on
plot(x, y, 'b-', LineWidth=1.5)
plot(x, y_exact, 'ro',  LineWidth=1.5)
xlabel('$x$','Interpreter','latex') % ,'fontsize',16
ylabel('$y(x)$','Interpreter','latex') % , 'Rotation',0 ,'Rotation',1
%axis([0. 3. 0.000 2.00])
set(gca,'FontSize',16)
hold off
box on

%%%
return
end
%

function dydx = rhs_function(x,y)
%
dydx = x + y;
%
return
end