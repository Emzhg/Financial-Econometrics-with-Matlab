function [LL] = LLgarch11(theta,r)
%LLgarch11 evaluate the LL of the C+GARCH(1,1) model

%***********************************************
% input
%*********************************************
% theta : vector(4,1) of unknown parameters,
% r : vector(T,1) of observed returns

%************************************************
%output
%**********************************************
%LL : value of the log-likelihood for theta and r

T = size(r,1); % NOBS
l= zeros(T+1,1);% vraisemblance à chaque date
zt= zeros(T+1,1); %standardised residuals
sigmat = zeros(T+1,1);% conditional variance

% etape 1 : tranformer les paramètres pour recuperer ceux du modèles

mu = theta(1,1)/10^3;
sigma2 = theta(4,1)/10^4;
alpha1 = theta(2,1);
beta1 = theta(3,1) - theta(2,1);
omega  = sigma2*(1-theta(3,1));

% calcul des résidus = rendements corrigés de leur moyenne
e = r -mu; 

% initialisation de 
sigmat(1,1) = sigma2; %la première valeur de la variance conditionnelle est la variance inconditionnelle
e = [ 0 ; e];% la première valeur des résidus est égal à leur espérance soit 0


% Calcul de la log-vraisemblance à chaque date t
for t=2:T+1
    sigmat(t,1)=omega+alpha1*e(t-1)^2+beta1*sigmat(t-1);
    zt(t,1)=e(t,1)/sqrt(sigmat(t,1));
    l(t,1) = -0.5*(log(2*pi)+ log(sigmat(t,1)) + zt(t,1)^2);
end

% On va utiliser l'algorithme fminunc de minimisation sous contrainte
% Maximiser LL  <=> minimiser - LL

LL = - sum(l(2:end)); % -Log-vraisemblance : fonction objectif à minimiser.
