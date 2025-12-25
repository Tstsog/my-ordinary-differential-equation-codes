%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code solve 1D ordinary differential equation (ODE) using the 4th-order Runge-Kutta method [1]. 
%
% The Newtonian 2nd law: f = m*a = m * dˆ2x/dtˆ2 = − k * x
%                    dx/dt = v          
%                    dv/dt = a = −k * x/m 
%
% The Runge-Kutta scheme for an equation:
%                        dy/dt = f(t,y) 
%                       y(i+1) = y(i) + (1/6) * (k1 + 2 * k2 + 2 * k3 + k4)
%                           k1 =  dt * f(t(i),          y(i))  
%                           k2 =  dt * f(t(i) + 0.5*dt, y(i) + 0.5*k1)  
%                           k3 =  dt * f(t(i) + 0.5*dt, y(i) + 0.5*k2)
%                           k4 =  dt * f(t(i) +     dt, y(i) +     k3)  
%
% Reference book: D. Bradie, A friendly introduction to numerical analysis 
%
% Written by Tsogbayar Tsednee (PhD), Institute of Physics and Technology, Mongolian Academy of Sciences 
%
% Contact: tsog215@gmail.com
% Date: December 25, 2025
%
%%%
function [] = Newtonian_2nd_law_rk4
%
clc; clear Newtonian_2nd_law_rk4
%
Nt = 512; % number of steps
tb = 0.; % initial value of x
te = 50.; % final value of x
dt = (te - tb)/Nt; % step size
%
t = zeros(1,Nt);
x = zeros(1,Nt); v = zeros(1,Nt); 
%
x(1) = 0.0; % initial position
v(1) = 0.5; % initial velosity
%
k = 0.5; % force constant
m = 1.5; % mass of object
%
for i = 1:Nt
    t(i+1) = tb + dt * i;
    %
    [dx_dt_1, dv_dt_1] = rhs_function(x(i), v(i), k, m);
    x_k1 = dt * dx_dt_1;
    v_k1 = dt * dv_dt_1;
    %
    [dx_dt_2, dv_dt_2] = rhs_function(x(i) + 0.5*x_k1, v(i) + 0.5*v_k1, k, m);
    x_k2 = dt * dx_dt_2;
    v_k2 = dt * dv_dt_2;
    %
    [dx_dt_3, dv_dt_3] = rhs_function(x(i) + 0.5*x_k2, v(i) + 0.5*v_k2, k, m);    
    x_k3 = dt * dx_dt_3;
    v_k3 = dt * dv_dt_3;
    %
    [dx_dt_4, dv_dt_4] = rhs_function(x(i) +     x_k3, v(i) +     v_k3, k, m);    
    x_k4 = dt * dx_dt_4;
    v_k4 = dt * dv_dt_4;
    %
    x(i+1) = x(i) + (x_k1/6) + (x_k2/3) + (x_k3/3) + (x_k4/6);
    v(i+1) = v(i) + (v_k1/6) + (v_k2/3) + (v_k3/3) + (v_k4/6);
    %
end
%
figure(1)
hold on
plot(t, x, 'b-', LineWidth=1.5)
plot(t, v, 'g-', LineWidth=1.5)
xlabel('$t$','Interpreter','latex') % ,'fontsize',16
ylabel('$x(t), v(t)$','Interpreter','latex') % , 'Rotation',0 ,'Rotation',1
%axis([0. 3. 0.000 2.00])
set(gca,'FontSize',16)
hold of
box on
%
figure(2)
plot3(t,x,v, 'b-', LineWidth=1.5)
xlabel('$t$','Interpreter','latex') % ,'fontsize',16
ylabel('$x$','Interpreter','latex') % ,'fontsize',16
zlabel('$v$','Interpreter','latex') % ,'fontsize',16
set(gca,'FontSize',16)
box on

%%%
return
end
%
function [dx_dt, dv_dt] = rhs_function(x,v,k,m)
%
dx_dt = v;
dv_dt = -(k/m).*x;
%
return
end