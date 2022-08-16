%--------------------------------------------------------------------------
% OBJET: Calcul des proba filtrées puis lissées - modèle MS de Vigfusson
%--------------------------------------------------------------------------
% DESCRIPTIF: [PL,PF] = Lissage (x);
% où: x = Valeur des Paramètres
% global endog : Variable endogene 
%--------------------------------------------------------------------------
% RESULTATS: 
%     PL = les probabilités lissées en x   (Tx1)
%     PF = les probabilités filtrées en x  (Tx1)
%--------------------------------------------------------------------------
% TP3 - 2021-22
%--------------------------------------------------------------------------


function [PL,PF] = lissage(x)

global endog

p11 =  normcdf(x(1));
[PF,PP]=filtrage(x);  % PF = Pr[St=1|It], PP = Pr(St=1|It-1)

T = size(endog,1);
PL = zeros(T,1);
PL(T,1) = PF(T,1);        % la derniere proba lissee=la derniere proba filtree

for t=T-1:-1:1
    PLi1 = p11 * PF(t) * PL(t+1) / PP(t+1);
    PLi2 = (1-p11) * PF(t) * (1-PL(t+1)) / (1-PP(t+1));
    PL(t,1)=PLi1+PLi2;
end

   
    

function [PR_STT1,PR_STL1]=filtrage(x)

global endog

   p11= normcdf(x(1)) ;  p22= normcdf(x(2));
   c= x(3:4);
   phi= x(5:6);
   VAR_L= x(7:8).^2;

% une matrice de probabilites de transition initiale 
   PR_TR=[p11 1-p22; 1-p11 p22];
   PR_TRF= [PR_TR(:,1); PR_TR(:,2)];
   
% calcul des probabilites ergodiques pour initialiser le filtre
   pi1 = (1-p22)/(2-p11-p22);   % probabilité inconditionnelle P(St=1)
   pi2 = (1-p11)/(2-p11-p22);   % probabilité inconditionnelle P(St=2)
     
   PROB__=   [pi1; pi1; pi2; pi2];  % vecteur des proba inconditionnelles (4x1)    

   T = length(endog);   % Taille de l'échantillon 
   LIKV=0.0; 
   PR_STL1 = zeros(T,1);  % vecteur de taille T qui contiendra P(St|It-1)
   PR_STT1 = zeros(T,1);  % vecteur de taille T qui contiendra P(St|It)
   
   for J_Iter=2:T  

     PROB_DD = PROB__ .* PR_TRF;   % Pr[St=i,St-1=j|It-1] 4x1  (étape 1 du filtre)

     PROB_DD = PROB_DD(1:2,1) + PROB_DD(3:4,1);  % Pr[St=i|It-1] 2x1  (étape 2 du filtre)

     PR_STL1(J_Iter,1) = PROB_DD(1,1);  % Pr(St=1|It-1) --- PP

     F_CAST = endog(J_Iter) - c - phi*endog(J_Iter-1);  % calcul des residus et|St (2x1)
     
     PR_VL = normpdf(F_CAST, 0, VAR_L.^(.5)) * PROB_DD';   % f(yt,St|It-1)  2x1 (étape 3 du filtre)

     PR_VAL = sum(PR_VL);    % f(yt|It-1), densite de yt sachant l'info passee (étape 4 du filtre)

     LIK=-1*log(PR_VAL);  % calcul de l'oppose de la log-vrais car fminunc mini les fonctions

     PROB__T = PR_VL ./ PR_VAL;   % Pr[St|It] 2x1 (étape 5 du filtre) 

     PR_STT1(J_Iter,1)=PROB__T(1,1);   % Pr(St=1|It) --- PF

     PROB__ = [PROB__T(1,1);PROB__T(1,1);PROB__T(2,1);PROB__T(2,1)];  
     % nv PROB__ utilise ds l'iteration suivante 4x1

     LIKV = LIKV + LIK;
     
   end