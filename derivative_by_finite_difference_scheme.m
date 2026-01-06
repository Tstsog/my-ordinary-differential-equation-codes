%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code computes the 1st and 2nd derivative of funcion f(x) using finite difference method[1]. 
%
% Function: fx = cos(x)*exp(x). 
%
% Reference book: D. Bradie, A friendly introduction to numerical analysis 
%
% Written by Tsogbayar Tsednee (PhD), Institute of Physics and Technology, Mongolian Academy of Sciences 
%
% Contact: tsog215@gmail.com
% Date: January 6, 2026
%
%%%
%%% 2026.1.6
function [] = derivative_by_finite_difference_scheme
%
clc; clear derivative_by_finite_difference_scheme
%
x0 = 0.5; % a point at which the numerical derivatives are computed
%
dx = 0.5; % a step size 
%
fx_1st_derivative_by_forward = (fx(x0+dx) - fx(x0))/dx; 
fx_1st_derivative_by_centeral = (fx(x0+dx) - fx(x0-dx))/(2*dx);
fx_1st_derivative_by_backward = (fx(x0) - fx(x0-dx))/(dx); 
%
fx_1st_derivative_analytic = -sin(x0)*exp(x0) + cos(x0)*exp(x0);
%%%

[dx, fx_1st_derivative_by_forward, fx_1st_derivative_by_centeral, fx_1st_derivative_by_backward, fx_1st_derivative_analytic]
%
[
%[dx, function_1st_derivative_by_forward, function_1st_derivative_by_centeral, function_1st_derivative_by_backward, function_1st_derivative_analytic]    
0.5000    0.0436    0.4687    0.8938    0.6564
0.2500    0.4084    0.6098    0.8111    0.6564
0.1250    0.5455    0.6448    0.7441    0.6564
0.0625    0.6041    0.6535    0.7030    0.6564
0.0312    0.6310    0.6557    0.6804    0.6564
0.0156    0.6439    0.6563    0.6686    0.6564
0.0078    0.6502    0.6564    0.6626    0.6564
];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%% the 2nd-order derivative by finite difference 
%
fx_2st_derivative_by_forward = (fx(x0+2*dx) - 2*fx(x0+dx) + fx(x0))/dx^2; 
fx_2st_derivative_by_centeral = (fx(x0+dx) - 2*fx(x0) + fx(x0-dx))/dx^2;
fx_2st_derivative_by_backward = (fx(x0) - 2*fx(x0-dx) + fx(x0-2*dx))/dx^2; 
%
fx_2st_derivative_analytic = -2*sin(x0)*exp(x0);
%
[dx, fx_2st_derivative_by_forward , fx_2st_derivative_by_centeral, fx_2st_derivative_by_backward, fx_2st_derivative_analytic]
%

[
%[dx, fx_2st_derivative_by_forward , fx_2st_derivative_by_centeral, fx_2st_derivative_by_backward, fx_2st_derivative_analytic]    
0.5000   -4.6939   -1.7003   -0.0833   -1.5809
0.2500   -2.9182   -1.6110   -0.6612   -1.5809
0.1250   -2.1941   -1.5884   -1.0729   -1.5809
0.0625   -1.8739   -1.5828   -1.3142   -1.5809
0.0312   -1.7240   -1.5813   -1.4443   -1.5809
0.0156   -1.6516   -1.5810   -1.5118   -1.5809
0.0078   -1.6160   -1.5809   -1.5461   -1.5809
];

%%%
return
end
%%%

function [fx] = fx(x)
%
fx = cos(x)*exp(x);
%%%
return
end