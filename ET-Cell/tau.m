for i = 1:length(V)
    for j = 1:length(Ca)
       X(i,j) = V(i);
       Y(i,j) = Ca(j);
       tau_matrix(i,j) = w_tau(V(i),Ca(j)); 
    end
end