%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code solve 1D ordinary differential equation (ODE) using the 4th-order Runge-Kutta method [1]. 
%
% 1D ODE:  dx/dt = t/x, x(0) = 1 and t = [0,5];
% Runge-Kutta method: y(i+1) = y(i) + (1/6) * (k1 + 2 * k2 + 2 * k3 + k4)
%                         k1 =  dt * f(t(i),          y(i))  
%                         k2 =  dt * f(t(i) + 0.5*tx, y(i) + 0.5*k1)  
%                         k3 =  dt * f(t(i) + 0.5*tx, y(i) + 0.5*k2)
%                         k4 =  dt * f(t(i) +     tx, y(i) +     k3)  
%
% An exact solution:  x(t) = (1 + t*t)^(1/2);
%
% Reference book: D. Bradie, A friendly introduction to numerical analysis 
%
% Written by Tsogbayar Tsednee (PhD), Institute of Physics and Technology, Mongolian Academy of Sciences 
%
% Contact: tsog215@gmail.com
% Date: December 25, 2025
%
%%%
function [] = rk4_example_3
%
clc; clear rk4_example_3
%
N = 5; % number of steps
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
%
for i = 1:N
    t(i+1) = tb + dt * i;
    x_k1 = dt * rhs_function(x(i),            t(i));
    x_k2 = dt * rhs_function(x(i) + 0.5*x_k1, t(i) + 0.5*dt);
    x_k3 = dt * rhs_function(x(i) + 0.5*x_k2, t(i) + 0.5*dt);
    x_k4 = dt * rhs_function(x(i) +     x_k3, t(i) +     dt);    
    %
    x(i+1) = x(i) + (x_k1/6.) + (x_k2/3.) + (x_k3/3.) + (x_k4/6.);
end
%
x_exact = (1. + t.*t).^(1/2);
[t', x', x_exact']
%
[        0    1.0000    1.0000
    1.0000    1.4190    1.4142
    2.0000    2.2394    2.2361
    3.0000    3.1646    3.1623
    4.0000    4.1249    4.1231
    5.0000    5.1005    5.0990];
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
%
function dxdt = rhs_function(x,t)
%
dxdt = t/x;
%
return
end
