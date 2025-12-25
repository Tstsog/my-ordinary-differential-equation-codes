%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code solve 1D coupled ordinary differential equation (ODE) using the 4th-order Runge-Kutta method [1]. 
%
% The Lorentz model: dx/dt = sigma * (y-x)
%                    dy/dt = r * x - y - x * z
%                    dz/dt = x * y - b * z
%
% The Runge-Kutta scheme for an equation:
%                        dy/dt = f(t,y) 
%                       y(i+1) = y(i) + (1/6) * (k1 + 2 * k2 + 2 * k3 + k4)
%                         k1 =  dt * f(t(i),          y(i))  
%                         k2 =  dt * f(t(i) + 0.5*tx, y(i) + 0.5*k1)  
%                         k3 =  dt * f(t(i) + 0.5*tx, y(i) + 0.5*k2)
%                         k4 =  dt * f(t(i) +     tx, y(i) +     k3)
%
% Reference book: D. Bradie, A friendly introduction to numerical analysis 
%
% Written by Tsogbayar Tsednee (PhD), Institute of Physics and Technology, Mongolian Academy of Sciences 
%
% Contact: tsog215@gmail.com
% Date: December 25, 2025
%
%%%
function [] = Lorentz_model_rk4
%
clc; clear Lorentz_model_rk4
%
N = 1000; % number of steps
tb = 0.; % initial value of x
te = 10.; % final value of x
dt = (te - tb)/N; % step size
%
t = zeros(1,N);
x = zeros(1,N); y = zeros(1,N); z = zeros(1,N);
%
t(1) = 0.;
x(1) = 1; y(1) = 1.; z(1) = 20.;
%
sigma = 10; r = 28; b = 8/3;
%
%
for i = 1:N
    t(i+1) = tb + dt * i;
    %
    [dx_dt_1, dy_dt_1, dz_dt_1] = rhs_function(x(i), y(i), z(i), sigma, r, b);
    x_k1 = dt * dx_dt_1;
    y_k1 = dt * dy_dt_1; 
    z_k1 = dt * dz_dt_1;
    %
    [dx_dt_2, dy_dt_2, dz_dt_2] = rhs_function(x(i) + 0.5*x_k1, y(i) + 0.5*y_k1, z(i) + 0.5*z_k1, sigma, r, b);
    x_k2 = dt * dx_dt_2;
    y_k2 = dt * dy_dt_2; 
    z_k2 = dt * dz_dt_2;
    %
    [dx_dt_3, dy_dt_3, dz_dt_3] = rhs_function(x(i) + 0.5*x_k2, y(i) + 0.5*y_k2, z(i) + 0.5*z_k2, sigma, r, b);
    x_k3 = dt * dx_dt_3;
    y_k3 = dt * dy_dt_3; 
    z_k3 = dt * dz_dt_3;
    %
    [dx_dt_4, dy_dt_4, dz_dt_4] = rhs_function(x(i) +     x_k3, y(i) +     y_k3, z(i) +     z_k3, sigma, r, b);
    x_k4 = dt * dx_dt_4;
    y_k4 = dt * dy_dt_4; 
    z_k4 = dt * dz_dt_4;
    %
    x(i+1) = x(i) + (x_k1/6) + (x_k2/3) + (x_k3/3) + (x_k4/6);
    y(i+1) = y(i) + (y_k1/6) + (y_k2/3) + (y_k3/3) + (y_k4/6);
    z(i+1) = z(i) + (z_k1/6) + (z_k2/3) + (z_k3/3) + (z_k4/6);
    %
end
figure(1)
hold on
plot(t, x, 'b-', LineWidth=1.5)
plot(t, y, 'g-', LineWidth=1.5)
plot(t, z, 'k-', LineWidth=1.5)
xlabel('$t$','Interpreter','latex') % ,'fontsize',16
ylabel('$x(t), y(t), z(t)$','Interpreter','latex') % , 'Rotation',0 ,'Rotation',1
%axis([0. 3. 0.000 2.00])
set(gca,'FontSize',16)
hold off
box on
%
figure(2)
plot3(x,y,z, 'b-', LineWidth=1.5)
xlabel('$x$','Interpreter','latex')
ylabel('$y$','Interpreter','latex')
zlabel('$z$','Interpreter','latex')
set(gca,'FontSize',16)
box on

%%%
return
end

%%%
function [dx_dt, dy_dt, dz_dt] = rhs_function(x, y, z, sigma, r, b)
%
dx_dt = sigma.*(y - x);
dy_dt = (r.*x - y - x.*z);
dz_dt = (x.*y - b.*z);
%
return
end