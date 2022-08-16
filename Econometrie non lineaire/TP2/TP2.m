clc;
close all;
clearvars;
clear;
%% TP2 : Xt suit un processus SETAR(1,1)
% DGP: Xt = voir la question 1
% Etape 1 : Identitification - supposer quer P1 = P2 = P
% Etape 2 : estimation
% Etape 3 : Test de spécification
% Etape 4: Test de linéarité

%% Exercice 1 : Simulation d un processus SETAR(1,1)


c1 = -0.5;
phi1 = 0.5;
sigma1 = 1;

c2 = 0.5;
phi2 = 0.7;
sigma2 = 1;

T = 500;
x =  zeros(T,1);
%indic_bas =  zeros(T,1);


for t = 2:1:500
    x(t) = (c1 + phi1*x(t-1)+sigma1*randn) * (x(t-1)<=0) + ...
        (c2 + phi2*x(t-1)+sigma2*randn)*(x(t-1)>0);
end

indic_bas(2:T) = x(1:T-1)<=0;

figure(1)
subplot(2,1,1), plot(x);
subplot(2,1,2);
    bar(indic_bas,1);

%indic_haut = zeros(T,1);
%indic_haut(2:T) = x(1:T-1)>0;

%% Etape 2 : Identification

% => On identifie un AR(1,1), Plin = 1;
% => On identifie SETAR(1,1) Hyp: p1 = P2 = p

figure(2)
subplot(2,1,1);autocorr(x)
subplot(2,1,2); parcorr(x)

%% Etape 3 : Estimation du processus SETAR
% seuil = 0; % Xt
% 
% y = x(2:T,1); 
% 
% cte = ones(T-1,1); % Constante c
% 
% x_1 = x(1:T-1,1) ; % Xt-1
% 
% indic = x_1 <= seuil;
% 
% X = [cte .*indic x_1.*indic cte.*(1-indic) x_1.*(1-indic)] ; % X
% bet = inv(X'*X)*X'*y; 
% bet2 = regress(y,X);
% e = y-X*bet;
% vare = e'*e/(T-1);

seuil = 0; % Xt

y = x(2:T,1); 

cte = ones(T-1,1); % Constante c

x_1 = x(1:T-1,1) ; % Xt-1

varmin = 10^10;
%seuil = 0; % Xt
seuil_tri=sort(x-1);
seuil_range = seuil_tri(round(0.15*(T-1)):round(0.85*(T-1)));

for seuil = seuil_range'
    indic = x_1 <= seuil;
    X = [cte .*indic x_1.*indic cte.*(1-indic) x_1.*(1-indic)] ; % X 
    bet = regress(y,X);
    e = y-X*bet;
    vare = e'*e/(T-1);
    if vare<varmin
        seuil_opt = seuil;
        e_opt = e;
        X_opt = X;
        varmin = vare;
    end
end

varbet = varmin * inv(X_opt'*X_opt);
sigbet= sqrt(diag(varbet));
tstat = bet_opt./ sigbet;
pval = 2*(1-normdf(abs(tstat)));
bet_dgp = [-0.5;0.5;0.5;0.7];
disp([bet_dgp bet_opt tstat pval]);

[~,pvalar] = lbqtest(e_opt, 'lags',[5 10 20]);
[~,pvalarch] = lbqtest(e_opt.^2, 'lags',[5 10 20]);
[~,pvalnorm] = jbtest(e_opt);











