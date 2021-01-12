%the leader model
function [Q_leader, c_leader] = leader(lambda_e, lambda_s, T, N)
	%length of z_t, y_t;
	len_zt = 2 ;
	y_t = [12,16] ;
    %variables
    z = 1 : len_zt * T;
    y = cell(sum(N),1);
    start = len_zt * T;
    for mm = 1 : length(N)      
            y{mm} = start + 1 : start + y_t(mm) * T;
            start = y{mm}(end); 
    end
    %sum of len_zt and y_t
    nu = len_zt * T;
    for mm = 1 : length(N)
            nu = nu  + y_t(mm) * T *N(mm);
    end
    %Q matrix
    Q_leader = zeros(nu * T, nu * T);
    for mm = 1 : N
        for nn = 1 : T
            Q_leader(z(1 + (nn - 1) * len_zt), y{mm}(3 +  (nn - 1) * y_t)) = -1 ;
            Q_leader(z(nn * len_zt), y{mm}(4 +  (nn - 1) * y_t)) = -1 ;
        end
    end
    %c matrix
    c_leader = zeros(nu * T, 1);
    for mm = 1 : N
        for nn = 1 : T
            c_leader(y{mm}(3 +  (nn - 1) * y_t), 1) = lambda_e(nn);
            c_leader(y{mm}(4 +  (nn - 1) * y_t), 1) = lambda_s(nn);
        end
    end
end