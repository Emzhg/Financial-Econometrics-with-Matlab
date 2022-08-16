function [LMstat,cvalueLM,pvalLM,hLM,Fstat,cvalueF,pvalF,hF] = ARCH_test(r,q,alpha)
%archtest : test d'un effet arch

%*****************************************************
%  inputs
%*****************************************************
% r = vecteur(T,1) des rendements
% q = nombre de retards de l'équation ARCH
% alpha = risque de première espece


%*****************************************************
% outputs
%*****************************************************
% LMstat : statistique du test du multiplicateur de Lagrange
% cvalueLM : seuil critique pour le test du multiplicateur de Lagrange (loi
% chi-deux(q))
% pvalLM : probabilité critique du test du multiplicateur de Lagrange
% hLM : 1 si rejet H0, 0 sinon pour le test du multiplicateur de Lagrange
% Fstat : statistique de Fisher
% cvalueF : seuil critique pour le test de Fisher
% pvalF : probabilité critique test de Fisher
% hF : 1 si rejet H0, 0 sinon.

esq  = (r -mean(r)).^2;

% Regression des rendements au carré sur leurs q valeurs retardés
X = lagmatrix(esq,1:q);
Reg = X(q+1:end,:);
Y = esq(q+1:end,1);
T = size(Y,1);%nombre d'observations disponibles pour la regression
[~,~,res,~,stat] = regress(Y,[ones(size(Y)) Reg]);

LMstat = T*stat(1); % LM_stat = T*R2
scr_nc = res'*res; % somme des carrés des résidus du modèle non contraint
scr_c = (esq-mean(esq))'*(esq-mean(esq));% somme des carrés des résidus modèle contraint
Fstat = ((scr_c - scr_nc)/scr_c)*((T-q-1)/q);

cvalueLM = chi2inv(1-alpha,q);
cvalueF = finv(1-alpha,q,T-q-1);

pvalF = 1 -fcdf(Fstat,q,T-q-1);
pvalLM = 1-chi2cdf(LMstat,q);

hLM = pvalF<alpha;
hF = pvalLM<alpha;
end


