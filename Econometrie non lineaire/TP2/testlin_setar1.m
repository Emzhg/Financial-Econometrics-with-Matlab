function wald = testlin_setar( x )
    T = length(x);
    y = x(2:T,1); cte = ones(T-1,1); x_1 = x(1:T-1);
	var_seuil = x_1; % variable de seuil 
    var_min = 10^10;
    rg_min = round(0.15*(T-1)); rg_max = round(0.85*(T-1));
    x_trie = sort(var_seuil);
    seuil_range = x_trie(rg_min:rg_max,1)'; 

    for seuil = seuil_range
    indic = var_seuil<=seuil;
	X = [cte.*indic x_1.*indic cte.*(1-indic) x_1.*(1-indic)];
	bet = inv(X'*X)*(X'*y);
	e = y - X*bet;
	var_res =e'*e/(T-1);

    if var_res<var_min 
        var_min = var_res;				
    % statistique de test de linearité 
        R=[1 0 -1 0;0 1 0 -1];
        r = zeros(2,1);
        wald = inv(var_min)*(R*bet-r)'*inv(R*inv(X'*X)*R')*(R*bet-r);
    end
    
    end

end

