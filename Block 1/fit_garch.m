function [theta_fin,sigmat] = fit_garch(theta,r)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here



T = size(r,1); % NOBS
sigmat = zeros(T+1,1);
e = zeros(T+1,1);

% etape 1 : tranformer les paramètres pour recuperer ceux du modèles

mu = theta(1,1)/10^3;
sigma2 = theta(4,1)/10^4;
alpha1 = theta(2,1);
beta1 = theta(3,1) - theta(2,1);
omega  = sigma2*(1-theta(3,1));

theta_fin = [mu omega alpha1 beta1]';


% initialisation de 
sigmat(1,1) = sigma2; %la première valeur de la variance conditionnelle est la variance inconditionnelle
e = [ 0 ; r-mu];%


for t=2:T+1
    sigmat(t,1)=omega+alpha1*e(t-1)^2+beta1*sigmat(t-1);
    %zt(t,1)=e(t,1)/sqrt(sigmat(t,1));
    %l(t,1) = -0.5*(log(2*pi)+ log(sigmat(t,1)) + zt(t,1)^2);
%LL = LL + l(t,1);

end
sigmat = sigmat(2:end);

end

