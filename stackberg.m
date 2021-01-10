clear;
N = [1];
a =[10,20,50];
b =[0.01,0.01,0.01];
c =[10,20,50];
d =[0.01,0.01,0.01];
Cs = [75, 72, 70];
Ce = [75, 72, 70];
Esi_max = [100, 60, 50];
Esi_min = [15, 10, 5];
Psi_max = [80, 100, 70];
Eei_max = [100, 60, 50];
Eei_min = [15, 10, 5];
Pei_max = [80, 100, 70];




B = cell(N,1);
U = cell(N,1);
r = [];
for mm = 1 : length(N)
    for nn = 1 : N(mm)
         [Qi, ci, Ai, di, Ei, Mi, Bi, ri, nu, ineq, T, Ui] = ship(a(nn), b(nn), Cs(nn), Esi_max(nn), Esi_min(nn), Psi_max(nn), N);
         %[Qi, ci, Ai, di, Ei, Mi, Bi, ri, nu, ineq, T, Ui] = AES(c(nn), d(nn), Ce(nn), Cs(nn), Esi_max(nn), Eei_min(nn), Pei_max(nn), Esi_min(nn), Psi_max(nn));
         B{nn} = Bi; 
         U{nn} = Ui;
         r = [r;ri];
    end
end


R = zeros(sum(N) * length(B{1}(:,1)),  length(U{1}(1,:)) + sum(N) * length(B{1}(1,:)));
for mm = 1 : sum(N)
    B_temp = B{mm};
    R(1 + (mm - 1) * length(B_temp(:,1)) : mm * length(B_temp(:,1)), 1:length(U{1}(1,:))) = U{mm};
    R(1 + (mm - 1) * length(B_temp(:,1)) : mm * length(B_temp(:,1)), 2 + 1 + (mm - 1) * length(B_temp(1,:)) : 2 + mm * length(B_temp(1,:))) = B{mm};
end
T = 2;
lambda_e = (1:24);
lambda_s = (1:24);
[Q_leader, c_leader] = leader(lambda_e, lambda_s, T, N);






%the leader model
% params.OutputFlag = 0;
params.NonConvex = 2;
clear models
modelx.modelsense = 'min';
modelx.obj = c_leader;
modelx.Q = sparse(Q_leader);
modelx.A = sparse(R);
modelx.rhs = r;
modelx.sense = repmat('<',length(modelx.rhs ),1);
modelx.vtype = [repmat('C',[2 * T, 1]);repmat([repmat('C',[8, 1]);repmat('B',[4, 1])], [T * N,1] )];
gurobi_write(modelx, 'modelx.lp')
resultx = gurobi(modelx, params);
x1 = resultx.x;


%the original model
params.OutputFlag = 0;
clear models
model.modelsense = 'min';
model.obj = ci;
model.Q = sparse(Qi);
model.A = sparse(Ai);
model.rhs = di;
model.sense = repmat('<',length(model.rhs ),1);
gurobi_write(model, 'mymodel.lp')
result = gurobi(model,params);
x0 = result.x;








%the lagrange model
params.OutputFlag = 0;
clear models
modelc.modelsense = 'min';
modelc.obj = 1*ones(length(Bi(1,:)),1);
modelc.A = sparse(Bi);
modelc.rhs = ri;
modelc.sense = repmat('<',length(modelc.rhs ),1);
modelc.vtype = repmat([repmat('C',[nu + ineq,1]);repmat('B',[ineq, 1])], [T,1] );
gurobi_write(modelc, 'modelc.lp')
resultc = gurobi(modelc, params);
x1 = resultc.x;









%%
%the leader model
function [Q_leader, c_leader] = leader(lambda_e, lambda_s, T, N)
    %variables
    z = 1 : 2 * T;
    y = cell(N);
    start = 2 * T;
    for mm = 1 : N
        y{mm} = start + 1 : start + 12 * T;
        start = y{mm}(end);
    end
    
    %the number of inequality constrains 
    nu = 14;
    %Q matrix
    Q_leader = zeros(nu * T, nu * T);
    for mm = 1 : N
        for nn = 1 : T
            Q_leader(z(1 + (nn - 1) * 2), y{mm}(3 +  (nn - 1) * 12)) = -1 ;
            Q_leader(z(nn * 2), y{mm}(4 +  (nn - 1) * 12)) = -1 ;
        end
    end
    %c matrix
    c_leader = zeros(nu * T, 1);
    for mm = 1 : N
        for nn = 1 : T
            c_leader(y{mm}(3 +  (nn - 1) * 12), 1) = lambda_e(nn);
            c_leader(y{mm}(4 +  (nn - 1) * 12), 1) = lambda_s(nn);
        end
    end
end


%the conventional ship
function [Qi, ci, Ai, di, Ei, Mi, Bi, ri, nu, ineq, T, Ui] = ship(a, b, Cs, Esi_max, Esi_min, Psi_max, N)
    T = 2;
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
    ineq = 4;
    ineq_sym = {1,-1,-1,1};
    ineq_x = {Es,Es,Ps,Ps};
    ineq_b = [Esi_max,Esi_min,0,Psi_max];
    nv = (nu +  ineq + ineq);
    %Q matrix
    Qi = zeros(nu * T, nu * T);
    for mm = 1 : T
        Qi(Es + (mm - 1) * nu, Es + (mm - 1) * nu) = b /2 ;
    end
    %c matrix
    ci = zeros(nu * T, 1);
    for mm = 1 : T
        ci(Es + (mm - 1) * nu, 1) = -a;
        ci(Ps + (mm - 1) * nu, 1) = Cs;
    end
    %A matrix
    Ai = zeros(ineq * T, nu * T);
    for mm = 1 : T
        for nn = 1 : ineq
            Ai(nn + (mm - 1) * ineq, ineq_x{nn} + (mm - 1) * ineq) = ineq_sym{nn};
        end
    end
    %b matrix
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
    %Ci matrix
    Di = repmat([0;-a;0;Cs],[T, 1]);
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
    c = 1;
    Ui = zeros(length(ri(:,1)), 2 * N);
    for mm = 1 : T
        Ui(4 + (mm - 1) * (2 * nu + 2 * ineq + 2 * ineq), mm * 2) = -c;  
        Ui(8 + (mm - 1) * (2 * nu + 2 * ineq + 2 * ineq), mm * 2) = c;
    end
end

%the AES
function [Qi, ci, Ai, di, Ei, Mi, Bi, ri, nu, ineq, T, Ui] = AES(c, d, Ce,Cs, Eei_max, Eei_min, Pei_max, Psi_min, Psi_max)
    T = 3;
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
    %b matrix
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
    %Ci matrix
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
    c = 1;
    Ui = zeros(length(ri(:,1)), 2);
    for mm = 1 : T
        Ui(4 + (mm - 1) * (2 * nu + 2 * ineq + 2 * ineq),2) = -c;  
        Ui(8 + (mm - 1) * (2 * nu + 2 * ineq + 2 * ineq),2) = c;
    end
end



% % % % % % function [Qi, ci, Ai, di, Ei, Mi, Bi, ri, nu, ineq, T] = AES()
% % % % % %     T = 1;
% % % % % %     c = 10;
% % % % % %     d = 41;
% % % % % %     Ce = 50;
% % % % % %     Cs = 50;
% % % % % %     eta = 0.7;
% % % % % %     Eei_min = 10;
% % % % % %     Eei_max = 100;
% % % % % %     Ee0 = 10;
% % % % % %     M = 1e10;
% % % % % %     %variables
% % % % % %     Ee = 1;
% % % % % %     Es = Ee + 1;
% % % % % %     Pe = Es + 1;
% % % % % %     Ps = Pe + 1;
% % % % % %     u1 = Ps + 1;
% % % % % %     u2 = u1 + 1;
% % % % % %     u3 = u2 + 1;
% % % % % %     u4 = u3 + 1;
% % % % % %     I1 = u4 + 1;
% % % % % %     I2 = I1 + 1;
% % % % % %     I3 = I2 + 1;
% % % % % %     I4 = I3 + 1;
% % % % % %     %the number of inequality constrains 
% % % % % %     nu = 4;
% % % % % %     ineq = 4;
% % % % % %     ineq_sym = {[1,-eta], [-1, eta], eta, -eta};
% % % % % %     ineq_x = {[Ee, Pe],[Ee, Pe],Pe,Pe};
% % % % % %     ineq_b = [Ee0, -Ee0, Eei_max - Ee0, -Eei_min + Ee0];
% % % % % %     nv = (nu +  ineq + ineq);
% % % % % %     %Q matrix
% % % % % %     Qi = zeros(nu * T, nu * T);
% % % % % %     for mm = 1 : T
% % % % % %         Qi(Ee + (mm - 1) * nu, Ee + (mm - 1) * nu) = d /2 ;
% % % % % %     end
% % % % % %     %c matrix
% % % % % %     ci = zeros(nu * T, 1);
% % % % % %     for mm = 1 : T
% % % % % %         ci(Ee + (mm - 1) * nu, 1) = -c;
% % % % % %         ci(Pe + (mm - 1) * nu, 1) = Ce;
% % % % % %         ci(Ps + (mm - 1) * nu, 1) = Cs;
% % % % % %     end
% % % % % %     %A matrix
% % % % % %     Ai = zeros(ineq * T, ineq * T);
% % % % % %     for mm = 1 : T
% % % % % %         for nn = 1 : ineq
% % % % % %             Ai(nn + (mm - 1) * ineq, ineq_x{nn} + (mm - 1) * ineq) = ineq_sym{nn};
% % % % % %         end
% % % % % %     end
% % % % % %     %b matrix
% % % % % %     di = zeros(ineq * T, 1);
% % % % % %     for mm = 1 : T
% % % % % %         for nn = 1 : ineq
% % % % % %             di(nn + (mm - 1) * ineq, 1) = ineq_b(nn);
% % % % % %         end
% % % % % %     end
% % % % % %     %Ei matrix
% % % % % %     Ei = diag(ones(ineq * T, 1));
% % % % % %     %Mi matrix
% % % % % %     Mi = M * Ei;
% % % % % %     %TODO
% % % % % %     %B matrix
% % % % % %     Bi = zeros((2 * nu + 2 * ineq + 2 * ineq) * T, nv);
% % % % % %     for mm = 1 :T
% % % % % %         Q_temp = Qi(1 + (mm - 1) * nu : mm * nu, 1 + (mm - 1) * nu : mm * nu);
% % % % % %         A_temp = Ai(1 + (mm - 1) * ineq : mm * ineq, 1 + (mm - 1) * nu : mm * nu);
% % % % % %         M_temp = Mi(1 + (mm - 1) * ineq : mm * ineq, 1 + (mm - 1) * ineq : mm * ineq);
% % % % % %         E_temp = Ei(1 + (mm - 1) * ineq : mm * ineq, 1 + (mm - 1) * ineq : mm * ineq);
% % % % % %         Bi_row = (2 * nu + 2 * ineq + 2 * ineq);
% % % % % %         Bi(1 + (mm - 1) * Bi_row : mm * Bi_row, 1 + (mm - 1) * nv : mm * nv) =...
% % % % % %            [2*Q_temp  A_temp' zeros(nu, ineq);
% % % % % %             -2*Q_temp  -A_temp' zeros(nu, ineq);
% % % % % %             A_temp zeros(ineq, ineq) zeros(ineq, ineq);
% % % % % %             -A_temp zeros(ineq, ineq) M_temp;
% % % % % %             zeros(ineq, nu) E_temp -M_temp;
% % % % % %             zeros(ineq, nu) -E_temp zeros(ineq, ineq)];
% % % % % %     end
% % % % % %     %Ci matrix
% % % % % %     Di = repmat([-c; 0; Ce; Cs],[T, 1]);
% % % % % %     %ui matrix
% % % % % %     ri = zeros((2 * nu + 2 * ineq + 2 * ineq) * T, 1);
% % % % % %     for mm = 1 : T
% % % % % %         D_temp = Di(1 + (mm - 1) * nu : mm * nu,1);
% % % % % %         d_temp = di(1 + (mm - 1) * ineq : mm * ineq,1);
% % % % % %         M_temp = repmat(M, [ineq, 1]);
% % % % % %         ri(1 + (mm - 1) * (2 * nu + 2 * ineq + 2 * ineq) : mm * (2 * nu + 2 * ineq + 2 * ineq), 1) = ...
% % % % % %            [-D_temp;
% % % % % %             D_temp;
% % % % % %             d_temp;
% % % % % %             -d_temp + M_temp;
% % % % % %             zeros(ineq,1);
% % % % % %             zeros(ineq,1)];
% % % % % %     end
% % % % % % %     %U matrix
% % % % % % %     Ui = zeros(length(ri(:,1)), 2);
% % % % % % %     for mm = 1 : T
% % % % % % %         Ui(3 + (mm - 1) * (2 * nu + 2 * ineq + 2 * ineq),1) = c1; 
% % % % % % %         Ui(4 + (mm - 1) * (2 * nu + 2 * ineq + 2 * ineq),1) = c1;  
% % % % % % %         Ui(7 + (mm - 1) * (2 * nu + 2 * ineq + 2 * ineq),2) = c2;
% % % % % % %         Ui(8 + (mm - 1) * (2 * nu + 2 * ineq + 2 * ineq),2) = c2;
% % % % % % %     end
% % % % % % end


%The leader model
% model.obj = zeros(nv * N *T, 1);
% model.A = B;
% model.rhs = b;

