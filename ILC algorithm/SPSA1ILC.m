% this function is the implementation of the paper [18]
function [U, E_SPSA_set_norm] = SPSA1ILC(R, U_ini, m, N)
aj = 1;
cj = 1;
delta = 2 * round(rand(m * N,1)) - 1;
Uplus = U_ini + cj * delta;
U_SPSA_set = zeros(m * N, 1);
U_SPSA_set(:, 1) = Uplus;
E_SPSA_set = [];
E_SPSA_set_norm = [];
for i = 1 : 1
    Y = lifted_model(U_SPSA_set(:, i), N);
    E = R - Y;
    E_SPSA_set = [E_SPSA_set, E];
    E_SPSA_set_norm = [E_SPSA_set_norm; 0.5 * norm(E)^2];
end
ghat = E_SPSA_set(:, 1)./(2 * cj * delta);
U = U_ini - aj * ghat;
