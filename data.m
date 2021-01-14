function [] = data(N) 
    a = repmat(10, [1, sum(N)]); 
    b = repmat(0.01, [1, sum(N)]); 
    c = repmat(20, [1, sum(N)]); 
    d = repmat(0.01, [1, sum(N)]); 
    Cs = repmat(75, [1, sum(N)]); 
    Ce = repmat(75, [1, sum(N)]); 
    Esi_max = repmat(100, [1, sum(N)]); 
    Esi_min = repmat(15, [1, sum(N)]); 
    Psi_max = repmat(80, [1, sum(N)]); 
    Eei_max = repmat(100, [1, sum(N)]); 
    Eei_min = repmat(15, [1, sum(N)]); 
    Pei_max = repmat(80, [1, sum(N)]); 
    save data;
end