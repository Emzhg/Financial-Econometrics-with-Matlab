%script Fiche 1 : Detection effet arch et Estimation d'un mod�le GARCG(1,1) M2 IEF

% step 1 : chargement des donn�es.
%clear;
clc;


[data]=xlsread('SP500_D.xlsx'); % Il faudra peut �tre changer le chemin d'acc�s


date = data(:,1); % matrice des dates - format num�rique
matlab_date = x2mdate(date); % conversion date format excel � matlab
indice = data(:,6); % selection de l'indice Adj Close

%%
%********************************************************
% Partie 1 Etude descriptive des rendements
%********************************************************


% Q1 : calcul des rendements + graphique

r=diff(log(indice));

figure(1)
plot(matlab_date(2:end),r)
datetick('x','mm/yyyy')
title('Rendement du S&P 500')
axis('tight')
ylim([-0.15 0.15])

%%
% Q2 : Calcul de statistiques descriptives
T = size(r,1);% nombre d'observations

% Calcul de la moyenne

m1 = mean(r);% fonction matlab
e = ones(size(r));
m2 = e'*r/size(r,1);% formule de calcul de la moyenne empirique

% Calcul de l'�cart type 

std1 = std(r); % fonction matlab
std2= sqrt((r-m2)'*(r-m2)/size(r,1));% formule de calcul de l'�cart type empirique


% Calcul du skewness

S1 = skewness(r); % fonction matlab
S2= ((r-m2).*(r-m2))'*(r-m2)/(size(r,1)*std2^3); % formule de calcul du skewness empirique
z_r = (r-m2)/std1; % rendements normalis�s standardis�s
S3 = ((z_r.*z_r)'*z_r)/(size(z_r,1)-1); % formule de calcul du skewness empirique


% Calcul du kurtosis
K1 = kurtosis(r);% fonction matlab
K2 = ((r-m2).*(r-m2))'*((r-m2).*(r-m2))/(size(r,1)*std2^4);  % formule de calcul du kurtosis empirique

% Tableau des statistiques descriptives
Stat={'Moyenne';'Ecart type';'Skewness';'Kurtosis'};
fonction_matlab = [m1;std1;S1;K1];
est_empirique = [m2;std2;S2;K2];

Stat_des =table(Stat,fonction_matlab,est_empirique);
disp(Stat_des);

% Histogramme
figure(2)
histogram(r)
title('Histogramme des rendements du S&P 500')

%%
% Q3 : Test de l'hypoth�se de normalit� de Jarques et Bera

[h,p,jbstat,critval] = jbtest(r); %fonction matlab

if h==1
    disp('Rejet hypothese loi normale')
else
    disp('Acceptation hypothese loi normale')
end

% statistique du test de Jarques et Bera

JB_test = (T)*S2^2/6+(T)*(K2-3)^2/24;
Crit_valJBtest = chi2inv(0.95,2);
Pval_JBtest = 1-chi2cdf(JB_test,2);

formatSpec='Statistique du test de Jarque et Bera = %4.2f et sa probabilit� critique = %4.2f ';
fprintf(formatSpec,JB_test,Pval_JBtest)
if Pval_JBtest<0.05
     disp('Rejet hypothese loi normale')
else
   disp('Acceptation hypothese loi normale')
end



%%
% Q4 : Calcul des rendements au carr� comme proxy de la variance � chaque
% date
% Objectif : montrer une alternance de p�riode de forte/faible variance
figure(3)
rd = (r-m1).^2;
plot(matlab_date(2:end),rdsq)
datetick('x','mm/yyyy')
title('Rendement au carr� du S&P 500')
axis('tight')
ylim([0 0.018])

%%
%%************************************************************************
% Partie 2 : test d'une effet arch
%***********************************************************************
q = 2; % nombre de retards dans la regression d'un effet ARCH
alpha = 0.05; %risque de premiere esp�ce

[LMstat,cValueLM,pvalLM,hLM,Fstat,cValueF,pvalF,hF] = ARCH_test(r,q,alpha); % construction d'une routine matlab
if hLM==1
    dLM= {'effet ARCH'};
else
    dLM ={'pas effet ARCH'};
end

if hF==1
    dF= {'effet ARCH'};
else
    dF={'pas effet ARCH'};
end

[h,pValue,stat,cValue] = archtest(r);% fonction matlab

% affichage des r�sultats

test = {'LMstat'; 'Fstat'};
valeur_stat=[LMstat ;Fstat];
seuil_rejet =[cValueLM ;cValueF];
pvalue = [pvalLM ; pvalF];
effet_arch = logical([hLM; hF]);

archtest_tab = table(test,valeur_stat,seuil_rejet,pvalue,effet_arch);

disp(archtest_tab);

%%
%**************************************************************************
% Partie 3 : Estimation des param�tres d'un mod�le GARCH(1,1)
%**************************************************************************

% M�thode d'estimation  : Estimateur du maximum de vraisemblance

% Q2 : valeurs initiales pour les parametres � estimer

alpha1 = 0.15;
beta1 = 0.80;
theta(1,1) =mean(r)*10^3;
theta(2,1)= alpha1;
theta(3,1)=alpha1+beta1;
theta(4,1)=var(r)*10^4;

%%
[LL] = LLgarch11(theta,r); %LLgarch11 : fonction de vraisemblance �valu�e pour les valeurs initiales des parametres et les rendements r
%cette ligne sert juste � v�rifier que la routine qui calcule (l'oppos� de)
%la log-vraisemblance fonctionne

%%
% Estimation des param�tres
% On cherche la valeur des param�tres qui maximise la log-vraisemblance :
% revient � chercher les param�tres qui minimise l'oppos� de la
% log-vraisemblance
options= optimoptions('fminunc','MaxIterations',500,'Display','iter'); %fminunc : algorithme de matlab pour la minimisation d'une fonction objectif sans contrainte sur les param�tres
f=@(x)LLgarch11(x,r);
[x,fval] = fminunc(f,theta,options);
% x= estimation des param�tres inconnus
% fval = valeur de la fonction objectif � l'optimum (= - max
% Log-vraisemblance)
%%
% Q4 : calcul des param�tres du GARCH(1,1) et des variances conditionnelles
% � partir de la solution optimale x
[theta_fin,sigmasqr] = fit_garch(x,r);

figure(4)
plot(matlab_date(2:end),sigmasqr)
datetick('x','mm/yyyy')
title('Variance conditionnelle de rendement du S&P 500')
% axis('tight')
% ylim([-0.15 0.15])

%%
clc;
formatSpec='Equation de la moyenne conditionnelle\n Constante = %1.4e\n';
fprintf(formatSpec,theta_fin(1));
formatSpec2 = 'Estimation du mod�le GARCH(1,1)\n constante = %1.4e\n alpha1 = %1.4f\n beta1 =%1.4f\n';
fprintf(formatSpec2,theta_fin(2),theta_fin(3),theta_fin(4));

%%
% Utilisation des fonctions de matlab pour estimer un mod�le GARCH(1,1)
Mdl=garch('GARCHLags',1,'ARCHLags',1,'offset',NaN);
[EstMdl,EstParamCov,logL,info] = estimate(Mdl,r);
summarize(EstMdl)

