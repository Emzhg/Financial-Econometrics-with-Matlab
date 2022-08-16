%------------------------------------------------------------------
% OBJET: Estimation et représentation des PF et PL d'un MS-AR(1)
%------------------------------------------------------------------
% DESCRIPTIF: programme principal qui fait appel à
% 1. la fonction VRAIS
% 2. la fonction LISSAGE
%------------------------------------------------------------------
% TP3 - Année 2021-22
%------------------------------------------------------------------


clc; clear; 
global endog


%% a. Simulation du processus
T = 200; 
c1 = -0.5; c2 = 0.5; phi1 = 0.7; phi2 = 0.5; sig1 = 0.5; sig2 = 1;
P = [0.9 0.2;0.1 0.8]; p11 = P(1,1); p22 = P(2,2);
theta = [c1 phi1 sig1;c2 phi2 sig2];
[endog,s] = simul_msar(theta,P,T);


%% b. Choix du vecteur de parametres initiaux 
x_dgp = [norminv(p11) norminv(p22) c1 c2 phi1 phi2 sig1 sig2]';
x0 = x_dgp + randn(8,1)./100;  % on ajoute un vecteur de tirages aléatoires

% p11,p22,c1,c2,phi1,phi2,sig1,sig2
% rmq : attention, on donne l'ecart-type et non la variance 


%% c. Estimation par maximum de vraisemblance des paramètres

% Options de l'Optimisation
options=optimset('MaxIter', 10000, 'TolFun', 10^-10, 'TolX', 10^-10);

% Appel de l'algorithme d'optimisation numérique
[x, fval, code, info, g, H]=fminunc(@vrais, x0, options);



%% d. Affichage des résultats d'estimation

clc ; format compact ;
disp('RESULTATS D''ESTIMATION DU MODELE');

disp('le code de retour est');
disp(code);
disp('la valeur de la log-vraisemblance maximisee est');
disp(-fval);
t_stat = x ./ sqrt(diag(inv(H)));
p_val = 2*(1-normcdf(abs(t_stat)));
disp('--------------------------------------------');
disp('---   DGP     estim     tstat      p-val ---');
disp('--------------------------------------------');
disp([x_dgp, x, t_stat, p_val]);

disp('--------------------------------------------');

disp('Matrice de transition estimée')
mat_transition = [normcdf(x(1)) 1-normcdf(x(2)); 1-normcdf(x(1)) normcdf(x(2))];
disp(mat_transition);

%% e. Représentation des probabilités lissées
[PL,PF] = lissage(x);

subplot(3,1,1), plot(PL), title('Proba Lissée');
subplot(3,1,2), plot(PF), title('Proba Filtrée');
etat_1 = PL > (mat_transition(1,2)/(mat_transition(2,1) + mat_transition(1,2)));
subplot(3,1,3), bar(etat_1,1), title('Etat 1');

%% Test de diagnostic

eps = resid_gen(x);

[~, p_val_lbq] = lbqtest(eps, 'lags', [5 10 20]);
if p_val_lbq < 0.05, disp('Résidus non autocorrélés'), end
[~, p_val_arch] = lbqtest(eps.^2, 'lags', [5 10 20]);
if p_val_lbq < 0.05, disp('Résidus sans effet ARCH'), end
[~, p_val_jb] = jbtest(eps);
if p_val_lbq < 0.05, disp('Résidus normalement distribués'), end
%% Remarque : Interprétation du code de retour
% 1 Magnitude of gradient smaller than the specified tolerance. 
% 2 Change in x was smaller than the specified tolerance. 
% 3 Change in the objective function value was less than the specified tolerance. 
% 0 Number of iterations exceeded options.MaxIter or number of function evaluations exceeded options.FunEvals. 
% -1 Algorithm was terminated by the output function. 
% -2 Line search cannot find an acceptable point along the current search direction.