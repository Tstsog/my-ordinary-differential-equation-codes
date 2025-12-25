%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code solve 1D ordinary differential equation (ODE) using the Euler's method [1]. 
%
% 1D ODE:  dx/dt = 1 + x/t, x(1) = 1 and t = [1,6];
% Euler's method with finite difference scheme: x(i+1) = x(i) + dt * (1 + x(i)/t(i))
%
% An exact solution: x(t) = t * (1 + ln(t));
%
% Reference book: D. Bradie, A friendly introduction to numerical analysis 
%
% Written by Tsogbayar Tsednee (PhD), Institute of Physics and Technology, Mongolian Academy of Sciences 
%
% Contact: tsog215@gmail.com
% Date: December 25, 2025
%
%%%
function [] = euler_example_2
%
clc; clear euler_example_2
%
N = 10; % number of steps
tb = 1.; % initial value of x
te = 6.; % final value of x
dt = (te - tb)/N; % step size
%
t = zeros(1,N);
x = zeros(1,N);
%
t(1) = 1.;
x(1) = 1;
%
for i = 1:N
    t(i+1) = tb + dt * i;
    x(i+1) = x(i) + dt * (1 + x(i)/t(i));
end
%
x_exact = t.* (1. + log(t));
[t', x', x_exact']
%
[   1.0000    1.0000    1.0000
    1.5000    2.0000    2.1082
    2.0000    3.1667    3.3863
    2.5000    4.4583    4.7907
    3.0000    5.8500    6.2958
    3.5000    7.3250    7.8847
    4.0000    8.8714    9.5452
    4.5000   10.4804   11.2683
    5.0000   12.1448   13.0472
    5.5000   13.8593   14.8761
    6.0000   15.6193   16.7506];
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
