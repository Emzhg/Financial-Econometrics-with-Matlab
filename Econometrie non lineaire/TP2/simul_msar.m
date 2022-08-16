%--------------------------------------------------------------------------
%% OBJET: simulation du mod�le MS-AR(1)
%--------------------------------------------------------------------------
% DESCRIPTIF:  [endog,s] = simul_MSAR(x0,theta,P,T)
% INPUTS: 
% x0    : valeur initiale de X   (1x1)
% theta : coefficients (c1 c2 phi1 phi2 sig1 sig2)  (6x1)
% P : matrice de transition (2x2)
% T : taille de l'�chantillon
%--------------------------------------------------------------------------
% OUTPUTS: 
% x : variable simul�e suivant un MS avec switch de ts les coef : T x 1
% s : variable simul�e suivant une cha�ne de Markov d'ordre un  : T x 1
%--------------------------------------------------------------------------


function [x,s] = simul_msar(theta,P,T)

x = zeros(T,1); 

s = markov(P,T);

for t = 2:T
    eps = randn(1,1);
    x(t) = (theta(1) + theta(3)*x(t-1) + theta(5)*eps)*(s(t)==1) + ...
           (theta(2) + theta(4)*x(t-1) + theta(6)*eps)*(s(t)==2);
end



%--------------------------------------------------------------------------
%% OBJET: Simulation d'une chaine de Markov
%--------------------------------------------------------------------------
% DESCRIPTIF: s = markov(Ptrans,T)
% o�: Ptrans = matrice de transition qui somme a 1 en colonne
%     T      = taille de l'�chantillon
%--------------------------------------------------------------------------
% RESULTATS: 
% s = vecteur de taille nper (les �tats � chaque instant)
%--------------------------------------------------------------------------

function s = markov(Ptrans,T)

s = zeros(T,1);
s(1,1) = 1; % initialisation : �tat 1 en t=0

for t=2:T 
    
  u = rand(1,1);
  
  % si �tat 1 en t, �tat en t+1
  if s(t-1)==1 
      if u <= Ptrans(1,1)
        s(t) = 1;
      else
        s(t) = 2;
      end
      
  % si �tat 0 en t, �tat en t+1
  elseif s(t-1)==2
      if u <= Ptrans(2,2)
        s(t) = 2;
      else
        s(t) = 1;
      end
  end
    
end

