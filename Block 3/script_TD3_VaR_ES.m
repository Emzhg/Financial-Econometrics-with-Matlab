%addpath('C:\Users\ylepen\Documents\MATLAB\cours\econometrie finance')

%%
% Partie I : présentation et analyse descriptive des données
[data,header] = xlsread('apple_exxon.xlsx');

date=data(:,1);
apple=data(:,2);
exxon=data(:,3);

m_date =x2mdate(date);
r_apple = diff(log(apple));
r_exxon=diff(log(exxon));
m_date = m_date(2:end);% perte de la premiere observation calcul rendement

r_port = 0.5*r_apple+0.5*r_exxon;

%%
% graphiques des rendements
figure(1)
subplot(3,1,1),plot(m_date,r_apple)
title('Rendement apple')
datetick('x','mmyyyy')
subplot(3,1,2),plot(m_date,r_exxon)
datetick('x','mmyyyy')
title('Rendement exxon')
subplot(3,1,3),plot(m_date,r_port)
datetick('x','mmyyyy')
title('Rendement portfeuille')

%%
figure(2)
subplot(3,1,1),histogram(r_apple)
title('Rendement apple')
subplot(3,1,2),histogram(r_exxon)
title('Rendement exxon')
subplot(3,1,3),histogram(r_port)
title('Rendement portfeuille')

% histogrammes montrent :
% distribution asymétrique
% excès de kurtosis
%%
R = [r_apple r_exxon r_port];
m_R =mean(R); 
median_R = median(R);
min_R= min(R);
max_R=max(R);
var_R =var(R);
S_R =skewness(R);
K_R= kurtosis(R);

%%
figure(3)
subplot(2,1,1)
autocorr(r_apple)
subplot(2,1,2)
parcorr(r_apple)

figure(4)
subplot(2,1,1)
autocorr(r_exxon)
subplot(2,1,2)
parcorr(r_exxon)


figure(5)
subplot(2,1,1)
autocorr(r_port)
subplot(2,1,2)
parcorr(r_port)

% Commentaire :
% les 3 séries de rendements ne présentent pas d'autocorrélation ni
% d'autocorrélation partielle
% Pas nécessaire de les modéliser à l'aide d'un modèle ARMA

%%
[h,pValue,stat,cValue] = lbqtest(r_port,'Lags',[5,10,15,20]);
varnames={'Retard','LBQtest','seuil critique','pvalue','Rejet H0?'};

Lag=[5,10,15,20];
T_lbqtest = table(Lag',stat',cValue',pValue',logical(h'),'VariableNames',varnames);
disp('Test absence autocorrelation LB pour le rendement du portefeuille')
disp(T_lbqtest)



%%
% 2 Estimation de la VaR et de l'ES à partir de l'historique

T= size(r_port,1);
N=500;
alpha1=0.01;
alpha2=0.05;
V=zeros(T-N+1,2);
ES =zeros(T-N+1,2);

for i=N:T
    r=r_port(i-N+1:i);
    r_sort = sort(r,'ascend');
    V(i,1)= - r_sort(alpha1*N,1);
    V(i,2)= - r_sort(alpha2*N,1);
    ES(i,1)= -sum(r_sort(1:alpha1*N,1))/(alpha1*N);
    ES(i,2)= -sum(r_sort(1:alpha2*N,1))/(alpha2*N);
end

subplot(2,1,1)
plot(m_date(N:end),V(N:end,:))
legend('1%','5%')
title('VaR à 1% et 5 %')
datetick('x','mm/yy')
subplot(2,1,2)
plot(m_date(N:end),ES(N:end,:))    
legend('1%','5%')
title('ES à 1% et 5 %')
datetick('x','mm/yy')


%%
% 3 Estimation de la VaR et de l'ES à partir d'un modèle GARCH

% Les rendements ne présentent pas d'autocorrélation
% on va juste les corriger de la moyenne

r_port_m = r_port-mean(r_port);

% estimation d'un modèle TARCH avec loi normale et effet d'asymétrie
P=1;
O=1;
Q=1;

[PARAM,LL,HT,VCVROBUST,VCV] = tarch(r_port_m,P,O,Q);
t_stat_param = PARAM./sqrt(diag(VCVROBUST));

%%
% residu standardises
z = r_port_m./sqrt(HT);
plot(z)
autocorr(z)
autocorr(z.^2)
[h,pval,LMstat,Cvalue]=archtest(z)

%%
% calcul de la VAR

% quantile de la loi normale

V_garch(:,1) = -(mean(r_port)+norminv(0.01,0,1)*sqrt(HT));
V_garch(:,2) = -(mean(r_port)+norminv(0.05,0,1)*sqrt(HT));

ES_garch(:,1) = 100*normpdf(norminv(0.01,0,1))*sqrt(HT)-mean(r_port);
ES_garch(:,2) = 20*normpdf(norminv(0.05,0,1))*sqrt(HT)-mean(r_port);



subplot(2,1,1)
plot(m_date,V_garch)
legend('1%','5%')
title('VaR à 1% et 5 %')
datetick('x','mm/yy')
subplot(2,1,2)
plot(m_date,ES_garch)    
legend('1%','5%')
title('ES à 1% et 5 %')
datetick('x','mm/yy')

%%
% Q4 Estimation de la VaR et de l'ES au niveau désagrégé

% Q4.1
r_apple_m=r_apple -mean(r_apple);
r_exxon_m = r_exxon-mean(r_exxon);

R = [r_apple_m r_exxon_m];

%%
% Q4.2 Modele Diagonal BEKK(1,1)

%[PARAMETERS,LL,HT,VCV,SCORES] = bekk(DATA,DATAASYM,P,O,Q,TYPE,STARTINGVALS,OPTIONS)

[param_bekk,LL_bekk,HT_bekk,VCV_bekk,] = bekk(R,[],P,O,Q,'Diagonal');

%%
w = [0.5 0.5];
var_port_bekk=zeros(T,1);
for t=1:T
    var_port_bekk(t,1)=w*HT_bekk(:,:,t)*w'; % calcul de la variance du portefeuille à partir des variances des actifs du portefeuille
end

V_bekk(:,1) = -(mean(r_port)+norminv(0.01,0,1)*sqrt(var_port_bekk));
V_bekk(:,2) = -(mean(r_port)+norminv(0.05,0,1)*sqrt(var_port_bekk));

ES_bekk(:,1) = 100*normpdf(norminv(0.01,0,1))*sqrt(var_port_bekk)-mean(r_port);
ES_bekk(:,2) = 20*normpdf(norminv(0.05,0,1))*sqrt(var_port_bekk)-mean(r_port);

subplot(2,1,1)
plot(m_date,V_bekk)
legend('1%','5%')
title('VaR à 1% et 5 %')
datetick('x','mm/yy')
subplot(2,1,2)
plot(m_date,ES_bekk)    
legend('1%','5%')
title('ES à 1% et 5 %')
datetick('x','mm/yy')


%%
% modele DCC

%[PARAMETERS,LL,HT,VCV,SCORES,DIAGNOSTICS] = dcc(DATA,DATAASYM,M,L,N,P,O,Q,...

% GJRTYPE,METHOD,COMPOSITE,STARTINGVALS,OPTIONS)
                                                    
M=1;% retard dcc
L=0; % asymetrie DCC
N=1;


[param_dcc,LL_dcc,HT_dcc,VCV_dcc] = dcc(R,[],M,L,N,P,O,Q);

%%
w = [0.5 0.5];
var_port_dcc=zeros(T,1);
for t=1:T
    var_port_dcc(t,1)=w*HT_dcc(:,:,t)*w';
end

V_dcc(:,1) = -(mean(r_port)+norminv(0.01,0,1)*sqrt(var_port_dcc));
V_dcc(:,2) = -(mean(r_port)+norminv(0.05,0,1)*sqrt(var_port_dcc));

ES_dcc(:,1) = 100*normpdf(norminv(0.01,0,1))*sqrt(var_port_dcc)-mean(r_port);
ES_dcc(:,2) = 20*normpdf(norminv(0.05,0,1))*sqrt(var_port_dcc)-mean(r_port);


subplot(2,1,1)
plot(m_date,V_dcc)
legend('1%','5%')
title('VaR à 1% et 5 %')
datetick('x','mm/yy')
subplot(2,1,2)
plot(m_date,ES_dcc)    
legend('1%','5%')
title('ES à 1% et 5 %')
datetick('x','mm/yy')
                                                     

%%
% Comparaison des méthodes d'estimation

compteur_hist_1 = sum(r_port(N:end)<-V(N:end,1))/size(r_port(N:end),1);
compteur_hist_2 = sum(r_port(N:end)<-V(N:end,2))/size(r_port(N:end),1);


compteur_garch_1 = sum(r_port<-V_garch(:,1))/size(r_port,1);
compteur_garch_5 = sum(r_port<-V_garch(:,2))/size(r_port,1);

compteur_bekk_1 = sum(r_port<-V_bekk(:,1))/size(r_port,1);
compteur_bekk_5 = sum(r_port<-V_bekk(:,2))/size(r_port,1);


compteur_dcc_1 = sum(r_port<-V_dcc(:,1))/size(r_port,1);
compteur_dcc_5 = sum(r_port<-V_dcc(:,2))/size(r_port,1);

compteur_1 = [compteur_hist_1;compteur_garch_1;compteur_bekk_1;compteur_dcc_1];
compteur_5 = [compteur_hist_5;compteur_garch_5;compteur_bekk_5;compteur_dcc_5];

nom={'1%','5%'};
row = {'hist','garch','bekk','dcc'};
resultat=table(compteur_1,compteur_5,'VariableNames',nom,'RowNames',row)


