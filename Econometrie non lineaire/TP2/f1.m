function [y] = X(t)
if t-1 <= 0
    y=C1+Phi1*X(t-1)+Sig1*Res1(t) ;
elseif -1< x && x<=3
    y = 2-exp(x^2/100);
else 
    y = 1;
end