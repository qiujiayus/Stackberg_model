%the AES
function [Qi, ci, Ai, di, Ei, Mi, Bi, ri, nu, ineq, T, Ui] = AES(c, d, Ce, Cs, Eei_max, Eei_min, Pei_max, Psi_min, Psi_max, N, T)
	%M
    M = 1e10;
    %variables
    Ee = 1;
    Es = Ee + 1;
    Pe = Es + 1;
    Ps = Pe + 1;
    u1 = Ps + 1;
    u2 = u1 + 1;
    u3 = u2 + 1;
    u4 = u3 + 1;
    I1 = u4 + 1;
    I2 = I1 + 1;
    I3 = I2 + 1;
    I4 = I3 + 1;
    %the number of inequality constrains 
    nu = 4;
    ineq = 6;
    ineq_sym = {1,-1,-1,1,-1,1};
    ineq_x = {Ee,Ee,Pe,Pe,Ps,Ps};
    ineq_b = [Eei_max,Eei_min,0,Pei_max,0,Psi_max];
    nv = (nu +  ineq + ineq);
    %Q matrix
    Qi = zeros(nu * T, nu * T);
    for mm = 1 : T
        Qi(Ee + (mm - 1) * nu, Ee + (mm - 1) * nu) = d /2 ;
    end
    %c matrix
    ci = zeros(nu * T, 1);
    for mm = 1 : T
        ci(Ee + (mm - 1) * nu, 1) = -c;
        ci(Pe + (mm - 1) * nu, 1) = Ce;
        ci(Ps + (mm - 1) * nu, 1) = Cs;
    end
    %A matrix
    Ai = zeros(ineq * T, nu * T);
    for mm = 1 : T
        for nn = 1 : ineq
            Ai(nn + (mm - 1) * ineq, ineq_x{nn} + (mm - 1) * nu) = ineq_sym{nn};
        end
    end
    %d matrix
    di = zeros(ineq * T, 1);
    for mm = 1 : T
        for nn = 1 : ineq
            di(nn + (mm - 1) * ineq, 1) = ineq_b(nn);
        end
    end
    %Ei matrix
    Ei = diag(ones(ineq * T, 1));
    %Mi matrix
    Mi = M * Ei;
    %TODO
    
    %B matrix
    Bi = zeros((2 * nu + 2 * ineq + 2 * ineq) * T, nv);
    for mm = 1 :T
        Q_temp = Qi(1 + (mm - 1) * nu : mm * nu, 1 + (mm - 1) * nu : mm * nu);
        A_temp = Ai(1 + (mm - 1) * ineq : mm * ineq, 1 + (mm - 1) * nu : mm * nu);
        M_temp = Mi(1 + (mm - 1) * ineq : mm * ineq, 1 + (mm - 1) * ineq : mm * ineq);
        E_temp = Ei(1 + (mm - 1) * ineq : mm * ineq, 1 + (mm - 1) * ineq : mm * ineq);
        Bi_row = (2 * nu + 2 * ineq + 2 * ineq);
        Bi(1 + (mm - 1) * Bi_row : mm * Bi_row, 1 + (mm - 1) * nv : mm * nv) =...
           [2*Q_temp  A_temp' zeros(nu, ineq);
            -2*Q_temp  -A_temp' zeros(nu, ineq);
            A_temp zeros(ineq, ineq) zeros(ineq, ineq);
            -A_temp zeros(ineq, ineq) M_temp;
            zeros(ineq, nu) E_temp -M_temp;
            zeros(ineq, nu) -E_temp zeros(ineq, ineq)];
    end
    %Di matrix
    Di = repmat([-c; 0; Ce; Cs],[T, 1]);
    %ui matrix
    ri = zeros((2 * nu + 2 * ineq + 2 * ineq) * T, 1);
    for mm = 1 : T
        D_temp = Di(1 + (mm - 1) * nu : mm * nu,1);
        d_temp = di(1 + (mm - 1) * ineq : mm * ineq,1);
        M_temp = repmat(M, [ineq, 1]);
        ri(1 + (mm - 1) * (2 * nu + 2 * ineq + 2 * ineq) : mm * (2 * nu + 2 * ineq + 2 * ineq), 1) = ...
           [-D_temp;
            D_temp;
            d_temp;
            -d_temp + M_temp;
            zeros(ineq,1);
            zeros(ineq,1)];
    end
    %U matrix
    c1 = -1;
	c2 = -1;
    Ui = zeros(length(ri(:,1)), 2 *  T );
    for mm = 1 : T 
		Ui(3 + (mm - 1) * (2 * nu + 2 * ineq + 2 * ineq),1 + (mm - 1) * 2) = c1;
        Ui(4 + (mm - 1) * (2 * nu + 2 * ineq + 2 * ineq),mm * 2) = c2;  
		Ui(7 + (mm - 1) * (2 * nu + 2 * ineq + 2 * ineq),1 + (mm - 1) * 2) = c1;
        Ui(8 + (mm - 1) * (2 * nu + 2 * ineq + 2 * ineq),mm * 2) = c2;
    end
end