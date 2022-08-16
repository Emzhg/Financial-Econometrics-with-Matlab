function [y] = func1(x)
if x<= -1
    y=0;
elseif -1< x && x<=3
    y = 2-exp(x^2/100);
else 
    y = 1;
end