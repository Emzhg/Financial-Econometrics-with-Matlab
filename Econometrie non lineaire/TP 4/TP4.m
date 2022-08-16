 % TP4

clear; clc; close all hidden

% pib : Taux de croissance du PIB français (source : Insee) ;
% ipi : taux de croissance de la production manufacturière (source : Insee)
% L%indice de production industrielle a été converti en fréquence trimestrielle par moyenne des valeurs mensuelles sur chaque trimestre. Nous reportons son taux de croissance trimestriel.

%% Q1 - Importation des données
data = xlsread('data_TP4.xlsx'); % xlsread n'importe pas les dates naturellement sur windows
Y = data(:,1); % variable prévue : PIB
X = data(:,2); % variable explication : IPI
% n = le nb d'obs

ind_estim = (1:100)'; % Set d'observations d'entrainement
ind_test = (101:115)'; % Set d'observations de test


%% Q2 - Partition des données
rng(1) % permet de fixer les données pour une question de reproductivité
ind_train = false(100,1);
ind_val = true(100,1);

tirage =randsample(100,round(0.5*100));
ind_train(tirage) = true;
ind_val(tirage) = false;

%% Q3 - Arbre de régression
RT = fitrtree(X(ind_train),Y(ind_train)); % fonction de fit du tree dans MATLAB
view(RT,'Mode','graph')


%% Q4 - Profondeur de l'arbre
m = max(RT.PruneList)

%% Q5 - Prévision pour différentes profondeurs
pruneLevels = 0:1:m; % niveau d'élagage vont de 0,1,...,9
Yfit = predict(RT,X(ind_val),'SubTrees',pruneLevels) % donne un vecteur ligne donnant pour chaque colonne les prédictions pour chaque élagage
% fit sur la série val et non train
% matrice de dimension 50 x10


%% Q6 - Choix du niveau d'élagage optimal
mse = mean((Y(ind_val)- Yfit).^2) %(10x1)

min_mse = min(mse) % a= 0.1047 donc 3e colonne => a = 3-1e colone = 2
%L = loss(RT,X(ind_val),Y(ind_val),'Subtrees' ,pruneLevels) 
%Deux fonction équivalentes pour faire le choix du nvieau d'élagage optimal

disp([pruneLevels,mse])
[~,a0] = min(mse);
a = a0-1 % => a = 2

% 0 élagage on se trouve tout en haut et c'est l'arbre le plus complexe;
% quand élagage = 9 on est tout en bas de l'arbre

%% Q7 - Arbre final et prévision sur le test set : élagage

RT_final = prune(RT,'Level',a); % a = 2
view(RT_final,'Mode','graph') % voir l'arbre
Yfit_RT = predict(RT_final,X(ind_test)); % commande predict(input sur les obs de l'arbre élagé RT_final, sur ind_test)

% Rq : faire varier rng à la Q2 qui va faire varier les résultats





%% Q9 - Forêt aléatoire

% ind_test = 101:115 data pour tester nos prévisions sur cet ensemble là
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

% on cherche à prendre le MSE le plus petit; a partir du graph on voit que
% le MSE le plus petit est la courbe bleue de taille = 5
b = 5 % taille choisi

RF_final = TreeBagger(100,X(1:100),Y(1:100),'Method','regression','InBagFraction',0.5,...
    'MinLeafSize',b);
Yfit_RF = predict(RF_final,X(ind_test)); %ind_test : 101:115



%view(b.Trees{i},'Mode','graph')

% 'Method': 'regression (par défaut, classification) ;
% 'InBagFraction': la fraction d’observations tirées avec remise et utilisée dans le training set pour la construction de chaque arbre (par défaut, 1) ;
% 'NumPredictorsToSample': nombre entier donnant le nombre de prédicteurs retenus pour la construction des arbres (si ‘all’, bagging, sinon, random forest) ; par défaut, 1/3 des prédicteurs sont tirés au hasard parmi les p.
% 'MinLeafSize': taille minimale des feuilles
% OOBPred':'on' pour conserver l’information sur quelles observations sont écartées (out of bag) dans l’evaluation set pour chaque arbre.


%% Q10 - Prévision avec la forêt aléatoire sur le test set : Random Forest
% Obtention des résultats : 
% RT_Final = 0.1373
% RF_Final = 0.0911



%% Q11 - Comparaison des deux approches


mse_RT = mean((Y(ind_test)- Yfit_RT).^2) % = RT_Final = 0.1373
mse_RF = mean((Y(ind_test)- Yfit_RF).^2) % RF_Final = 0.0911 

% Résultats bons
% Combinaison de modèles => chercher la moyenne de nos modèles
% Prix de l'électricité : marche bien
















