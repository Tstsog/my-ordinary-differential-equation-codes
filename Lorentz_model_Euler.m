%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code solve 1D coupled ordinary differential equation (ODE) using the Euler's method [1]. 
%
% The Lorentz model: dx/dt = sigma * (y-x)
%                    dy/dt = r * x - y - x * z
%                    dz/dt = x * y - b * z
% Euler's method with finite difference scheme: 
%                    x(i+1) = x(i) + dt * sigma * (y(i) - x(i))
%                    y(i+1) = y(i) + dt * (r * x(i) - y(i) - x(i) * z(i))
%                    z(i+1) = z(i) + dt * (x(i) * y(i) - b * z(i)))
%
% Reference book: D. Bradie, A friendly introduction to numerical analysis 
%
% Written by Tsogbayar Tsednee (PhD), Institute of Physics and Technology, Mongolian Academy of Sciences 
%
% Contact: tsog215@gmail.com
% Date: December 25, 2025
%
%%%
function [] = Lorentz_model_Euler
%
clc; clear Lorentz_model_Euler
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
for i = 1:N
    t(i+1) = tb + dt * i;
    x(i+1) = x(i) + dt * sigma* (y(i) - x(i)); % for x
    y(i+1) = y(i) + dt * (r * x(i) - y(i) - x(i) * z(i)); % for y    
    z(i+1) = z(i) + dt * (x(i) * y(i) - b * z(i)); % for y        
end
%

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
