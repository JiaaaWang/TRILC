% this function is the implementation of the paper [11], the parameters are tuned to achieve the best performance
function [U, E_SPSA_set_norm] = SPSA2ILC(j, R, U_ini, m, N)
aj = 1.2 / (j^0.2);
cj = 0.5 / (j^0.25);
delta = 2 * round(rand(m * N,1)) - 1;
Uplus = U_ini + cj * delta;
Uminus = U_ini - cj * delta;
U_SPSA_set = zeros(m * N, 2);
U_SPSA_set(:, 1) = Uplus;
U_SPSA_set(:, 2) = Uminus;
E_SPSA_set = [];
E_SPSA_set_norm = [];
for i = 1 : 2
    Y = lifted_model(U_SPSA_set(:, i), N);
    E = R - Y;
    E_SPSA_set = [E_SPSA_set, E];
    E_SPSA_set_norm = [E_SPSA_set_norm; 0.5 * norm(E)^2];
end
ghat = (E_SPSA_set(:, 1) - E_SPSA_set(:, 2))./(2 * cj * delta);
U = U_ini - aj * ghat;