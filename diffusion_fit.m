function [FtFo_fit] = diffusion_fit(x,xdata)
% this function fits a fluorescence recovery curve to the diffusion model,
% with the focal volume dimensions as wz = 5.811um and wr = 0.6455um. This
% fitting equation is designed to be used with lsqcurvefit()
% Inputs : 
%   x: a vector containing the four fitting parameters
%   xdata: time vector 
% Outputs :
%   FtFo_fit: a vector of points that fit the fluorescence recovery
%   data
n = 0:1:10; % summation index 
R = 81.0418; % square of wz/wr

t = xdata;
td = x(1); % recovery due to diffusion, tauD
beta = x(2); % bleach depth parameter

for i = 1:1:length(n)
   FtFo_fit(:,i) = (((-beta)^n(i))./(factorial(n(i)))).*(1./(1+n(i)+2*n(i)*t/td)).*...
        (1./(sqrt(1+n(i)+2*n(i)*t/(R*td))));
end 
FtFo_fit = sum(FtFo_fit,2); % summation over n
end