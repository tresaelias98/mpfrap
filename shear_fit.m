function [FtFo_fit] = shear_fit(x,xdata)
% this function fits a fluorescence recovery curve to the shear flow model,
% with the focal volume dimensions as wz = 5.811um and wr = 0.6455um. This
% fitting equation is designed to be used with lsqcurvefit()
% Inputs : 
%   x: a vector containing the four fitting parameters
%   xdata: time vector 
% Outputs :
%   FtFo_fit: a vector of points that fit the fluorescence recovery
%   data

n = 0:1:20; % summation index
R = 81.0418; % square of wz/wr

td = x(1); % recovery due to diffusion, tauD
tv = x(2); % recovery due to flow, tauV
ty = x(3); % recovery due to shear, tauGamma
beta = x(4); % bleach depth parameter

t = xdata;

for i = 1:1:length(n)
    
    
    A(:,i) = (((-beta)^n(i))/factorial(n(i)))*...
        (1./(1+2*n(i)*t/td)).*(1./(sqrt(1+2*n(i)*t/(R*td))));
    
    J1(:,i) = 1./((n(i)./(1+2*n(i)*t/td))+1);
    J2(:,i) = 1./(sqrt(1+(n(i)./(1+2*n(i)*t./(R*td))+...
        (n(i)*R*(t/ty).^2)./(1+n(i)+2*n(i)*t./td))));
    J(:,i) = J1(:,i) .* J2(:,i);
    
    S1(:,i) = (-4*n(i)*(t/tv).^2)./(1+n(i)+2*n(i)*t/td);
    S2(:,i) = (n(i)*R*(t/ty).^2)./(1+n(i)+2*n(i)*t/td);
    S3(:,i) = 1+(n(i)./(1+2*n(i)*t/(R*td)))+((n(i)*R*(t/ty).^2)./(1+n(i)+2*n(i)*t/td));
    S(:,i) = exp(S1(:,i).*(1-(S2(:,i)./S3(:,i))));  
    
    FtFo_fit(:,i) = A(:,i).*J(:,i).*S(:,i);
    

end
FtFo_fit = sum(FtFo_fit,2); % summation over n 
end


