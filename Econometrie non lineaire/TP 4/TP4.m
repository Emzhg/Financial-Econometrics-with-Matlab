 % TP4

clear; clc; close all hidden

% pib : Taux de croissance du PIB fran�ais (source : Insee) ;
% ipi : taux de croissance de la production manufacturi�re (source : Insee)
% L%indice de production industrielle a �t� converti en fr�quence trimestrielle par moyenne des valeurs mensuelles sur chaque trimestre. Nous reportons son taux de croissance trimestriel.

%% Q1 - Importation des donn�es
data = xlsread('data_TP4.xlsx'); % xlsread n'importe pas les dates naturellement sur windows
Y = data(:,1); % variable pr�vue : PIB
X = data(:,2); % variable explication : IPI
% n = le nb d'obs

ind_estim = (1:100)'; % Set d'observations d'entrainement
ind_test = (101:115)'; % Set d'observations de test


%% Q2 - Partition des donn�es
rng(1) % permet de fixer les donn�es pour une question de reproductivit�
ind_train = false(100,1);
ind_val = true(100,1);

tirage =randsample(100,round(0.5*100));
ind_train(tirage) = true;
ind_val(tirage) = false;

%% Q3 - Arbre de r�gression
RT = fitrtree(X(ind_train),Y(ind_train)); % fonction de fit du tree dans MATLAB
view(RT,'Mode','graph')


%% Q4 - Profondeur de l'arbre
m = max(RT.PruneList)

%% Q5 - Pr�vision pour diff�rentes profondeurs
pruneLevels = 0:1:m; % niveau d'�lagage vont de 0,1,...,9
Yfit = predict(RT,X(ind_val),'SubTrees',pruneLevels) % donne un vecteur ligne donnant pour chaque colonne les pr�dictions pour chaque �lagage
% fit sur la s�rie val et non train
% matrice de dimension 50 x10


%% Q6 - Choix du niveau d'�lagage optimal
mse = mean((Y(ind_val)- Yfit).^2) %(10x1)

min_mse = min(mse) % a= 0.1047 donc 3e colonne => a = 3-1e colone = 2
%L = loss(RT,X(ind_val),Y(ind_val),'Subtrees' ,pruneLevels) 
%Deux fonction �quivalentes pour faire le choix du nvieau d'�lagage optimal

disp([pruneLevels,mse])
[~,a0] = min(mse);
a = a0-1 % => a = 2

% 0 �lagage on se trouve tout en haut et c'est l'arbre le plus complexe;
% quand �lagage = 9 on est tout en bas de l'arbre

%% Q7 - Arbre final et pr�vision sur le test set : �lagage

RT_final = prune(RT,'Level',a); % a = 2
view(RT_final,'Mode','graph') % voir l'arbre
Yfit_RT = predict(RT_final,X(ind_test)); % commande predict(input sur les obs de l'arbre �lag� RT_final, sur ind_test)

% Rq : faire varier rng � la Q2 qui va faire varier les r�sultats





%% Q9 - For�t al�atoire

% ind_test = 101:115 data pour tester nos pr�visions sur cet ensemble l�
for taille = 5:5:15
hold on
RF = TreeBagger(100,X(1:100),Y(1:100),'Method','regression','InBagFraction',...
    0.5,'OOBPred','On','MinLeafSize',taille);
plot(oobError(RF))
end
xlabel("nb d'arbres");
ylabel('mse')
legend('taille=5','taille = 10','taille = 10')
hold off

% on cherche � prendre le MSE le plus petit; a partir du graph on voit que
% le MSE le plus petit est la courbe bleue de taille = 5
b = 5 % taille choisi

RF_final = TreeBagger(100,X(1:100),Y(1:100),'Method','regression','InBagFraction',0.5,...
    'MinLeafSize',b);
Yfit_RF = predict(RF_final,X(ind_test)); %ind_test : 101:115



%view(b.Trees{i},'Mode','graph')

% 'Method': 'regression (par d�faut, classification) ;
% 'InBagFraction': la fraction d�observations tir�es avec remise et utilis�e dans le training set pour la construction de chaque arbre (par d�faut, 1) ;
% 'NumPredictorsToSample': nombre entier donnant le nombre de pr�dicteurs retenus pour la construction des arbres (si �all�, bagging, sinon, random forest) ; par d�faut, 1/3 des pr�dicteurs sont tir�s au hasard parmi les p.
% 'MinLeafSize': taille minimale des feuilles
% OOBPred':'on' pour conserver l�information sur quelles observations sont �cart�es (out of bag) dans l�evaluation set pour chaque arbre.


%% Q10 - Pr�vision avec la for�t al�atoire sur le test set : Random Forest
% Obtention des r�sultats : 
% RT_Final = 0.1373
% RF_Final = 0.0911



%% Q11 - Comparaison des deux approches


mse_RT = mean((Y(ind_test)- Yfit_RT).^2) % = RT_Final = 0.1373
mse_RF = mean((Y(ind_test)- Yfit_RF).^2) % RF_Final = 0.0911 

% R�sultats bons
% Combinaison de mod�les => chercher la moyenne de nos mod�les
% Prix de l'�lectricit� : marche bien
















