clear; close all; clc

%% Q1
A = [-1 4 0 ; 2 0 1 ; 2 6 2];
B = [2 1 5 ; 4 0 9 ; -1 3 6 ];
C = ones(2,2);

A*B; %Produit matriciel
A.*B; % Porduit élément par élément
A>1; %(zij) avec zij = 1 si aij >1 , 0 sinon
A(A>1)=5; %M(bij) % avec bij= 5 si oij > 1, aij sinon

Test1 = A*B;
Test2 = A.*B;
Test3 = A.^ 2;
Test4 = A>1;
%Test5 = A(A>1) 0

v1 = A(2,:);
v2 = B(:,1);


%% Q2

Testf1 = f1(1.2);
Testf1n2 = f1(5.7);


%% Q3

x = zeros(1000,1);
t=2;
while t<=1000
    x(t) = x(t-1) + randn;
    t= t+1;
end

figure(1)
plot(x)
title("Trajectoire d'une marche aléatoire")


%x = zeros(1000,1);
%for t=2:1000
%    x(t) = x(t-1)+randn;
%end
%plot(x)

%% Q4 Subplot

x = linspace(0,pi,30);
y1 = cos(x);y2 = sin(x); y3 = exp(-x);

%subplot()

figure(2)
subplot(3,1,1), plot(x,y1);
subplot(3,1,2), plot(x,y2);
subplot(3,1,3), plot(x,y3);

%% Q5 Fonction de Rosenbrock
e = randn(1000,1); x = zeros(1000,1);
rw =cumsum(e);
for t=2:1000
    x(t) = sum(e(1:t));
end


x0= randn(2,1);
options = optimset('MaxIter',100,'Tolfun',10^-10,'TolX',10^-10);
[x, fval,code,info,grad,hessien] = fminunc(@f2,x0,options)
% [x*,f(x*),>0,..., (Deltaf/deltax1, deltaf/deltax2), delta²f/
% deltax1*deltax2;delta²f/deltax2²)



%% Partie 2
% Q1 Simulation d'un processus AR(1)

clear all; clc;

%Yt = 0.5 + 0.8*Yt-1 + 0.9*Et
%Et üd N(0,1) => ut üd N(0,1)

c = 0.5;
phi1 = 0.8;
sigma1 = 0.9;

Xt =  zeros(500,1);

for t = 2:1:500
    Xt(t) = c + phi1*Xt(t-1)+sigma1*randn ;
end

figure(3)
plot(Xt)

%% Exerice 2 Identification de l ordre du processus

figure(4)
subplot(2,1,1);autocorr(Xt)
subplot(2,1,2); parcorr(Xt)


%% Exercice 3 Estimation du processus AR(1)

% Yt = c + phi1*Yt(t-1)+sigma1*Et <=> Y = X*Beta+Eps

T=500;
Y = Xt(2:T,1);
X = [ones(T-1,1) Xt(1:T-1,1)];
Beta=[c,phi1];

Bet = inv(X'*X)*(X'*y);
e = y - X * Bet;
Vare = e'*e/(T-1); size = sqrt(Vare);
Varbet = Vare *inv(X'*X);
%disp([[0.5;0.8;0.9] [Bet,size]])


%% Test de significativité des coefficients
sigbet = sqrt(diag(varbet));
tstat = bet./sigbet;
pval = 2 *(1-nomcdf(abs(tstat)));


%% Test de diagnostic des résidus (absence d'auto correlation AR)
[~,pval_ar] = lbqtest(e,'lags',[5;10;20]);

% h=0  indicates that there is not enough evidence to reject the null
% hypothesis that the residual
% H0 : p1=ph=0
% Ici si h= 20, pval = 34.75% > 5%, rejet de H0
% h_ar 0 si non rejet de h0 au seuil de alpha = 5%
% 1 sinon
disp('Test de Ljung Box sur les résidus');
disp(pValue)
[h_arch,pValue,stat,cValue] = lbqtest(resid.^2, 'lags',[5,10,20]);
disp('Test de Ljung Box sur les résidus au carré');
disp(pValue)
[h_norm,p,jbstat,critval] = jbtest(resid);
disp('test de Jarque Bera sur les résidus');

% Test de Jarque Bera :
% H0 : Epsilon suit une loi Normale vs H1:H0(barre)
% Stat de JB suit une loi de Chi(2)
% Conclusion : p-value = 46.41% > 5% => Non rejet de H0 au seuil de alpha =
% 5%




