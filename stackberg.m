
function stackberg(x,nv,N,T)
if nargin < 4
    T = 3;
    if nargin < 3
        N = 2;
        if nargin < 2
            nu = 4;
            nv = 12;
            if nargin < 1
                x = [1];
            end
        end
    end
end
a = 20;
b = 41;
Cs = 98;
Esi_max = 100;
Esi_min = 5;
Psi_max = 100;
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
ineq_sym = [1,-1,-1,1];
ineq_x = [Es,Es,Ps,Ps];
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
Ai = zeros(ineq * T, ineq * T);
for mm = 1 : T
    for nn = 1 : ineq
        Ai(nn + (mm - 1) * ineq, ineq_x(nn) + (mm - 1) * ineq) = ineq_sym(nn);
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
        A_temp zeros(ineq, ineq) M_temp;
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
resultc = gurobi(modelc,params);


%The leader model






% model.obj = zeros(nv * N *T, 1);
% model.A = B;
% model.rhs = b;

end