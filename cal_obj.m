function obj = cal_obj(Y, L, D)
    T = (Y' * D * Y ).^0.5;
    t = diag(T);
    t(t~=0) = t(t~=0).^-1;
    G = D^0.5 * Y * diag(t); 
    obj = trace(G' * L * G);
end