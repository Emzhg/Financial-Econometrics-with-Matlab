addpath(genpath('H:\MFE'))
addpath(genpath('H:\Econometrie de la volatilite'))
%%
[data,header] =xlsread('SP_Apple.xlsx','apple_sp');

date = data(:,1);
sp =data(:,2);
apple =data(:,3);

r_sp = diff(log(sp));
r_apple = diff(log(apple));

%%
%  Autocorrelogramme

figure (1)
autocorr(r_sp)
figure(2)
parcorr(r_sp)


figure (3)
autocorr(r_apple)
figure(4)
parcorr(r_apple)

%%
% Modelisation de l'espérance conditionnelle par un AR(1)
y= r_apple;
p=1;
q=0;
constant=1;
[parameters, LL, r_apple_f, SEregression, diagnostics, VCVrobust, VCV]=armaxfilter(y,constant,p,q);

autocorr(r_apple_f)

t_stat = parameters./sqrt(diag(VCVrobust))



%%
% Modélisation de l'espérance conditionnelle par un MA(2)
y= r_sp;
p=0;
q=[1 2];
constant=1;
[parameters, LL, r_sp_f, SEregression, diagnostics, VCVrobust, VCV]=armaxfilter(y,constant,p,q)

autocorr(r_sp_f)

t_stat = parameters./sqrt(diag(VCVrobust));
%% Question 2
R = [r_sp_f r_apple_f];

%% Question 4
% Estimation de la matrice de variance covariance par le modele bekk(1,1)
