%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code solve the Newtonian 2dn law (f = m * a) using the Euler's method [1]. 
%
% The Newtonian 2nd law: f = m*a = m * dˆ2x/dtˆ2 = − k * x
%
% Euler's method with finite difference scheme: 
% dx/dt = v            => x(i+1) = x(i) + dt * v(i)
% dv/dt = a = −k * x/m => v(i+1) = v(i) + dt *(−k/m) * x(i)
%
% Reference book: D. Bradie, A friendly introduction to numerical analysis 
%
% Written by Tsogbayar Tsednee (PhD), Institute of Physics and Technology, Mongolian Academy of Sciences 
%
% Contact: tsog215@gmail.com
% Date: December 25, 2025
%
%%%
function [] = Newtonian_2nd_law_Euler
%
clc; clear Newtonian_2nd_law_Euler
%
Nt = 512; % number of steps
tb = 0.; % initial value of x
te = 50.; % final value of x
dt = (te - tb)/Nt; % step size
%
t = zeros(1,Nt);
x = zeros(1,Nt); v = zeros(1,Nt); 
%
x(1) = 0.; % initial position
v(1) = 0.5; % initial velosity
%
k = 0.5; % force constant
m = 1.5; % mass of object
%
for i = 1:Nt
    t(i+1) = tb + dt * i;
    x(i+1) = x(i) + dt * v(i);
    v(i+1) = v(i) + dt * (-k/m) * x(i);
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
hold off
box on

%%%
return
end
