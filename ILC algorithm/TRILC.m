function [U, E_model_diff, E_SI_set_norm, Delta_new] = TRILC(R, U_ini, E_ini, Delta, m, N, ep, mu)
E_SI_set_norm_temp_store = []; % store the additional tracking errors due to the criticality phase
reduce_para = 0.9; % criticality radius reduction
while 1
    U_SI_set = zeros(m * N, m * N); % initialize the additioanl trial inputs
    U_SI_set_diff = zeros(m * N, m * N); % initialize the difference between additioanl and main trial inputs
    for i = 1 : m * N
        U_SI_set(:, i) = U_ini; % the additional trial input is generated based on the main trial input
        U_SI_set(i, i) = U_SI_set(i, i) + (floor(rand * 2) * 2 - 1) * Delta; % generate the additional trial inputs by (8)
        U_SI_set_diff(:, i) = U_SI_set(:, i) - U_ini; % calculate the matrix W in (7)
    end
    E_SI_set = []; 
    E_SI_set_norm = []; 
    for i = 1 : m * N
        Y = lifted_model(U_SI_set(:, i), N); % calculate the output trajectory from the additioanl input
        E = R - Y; % calculate the additioanl trial tracking error
        E_SI_set = [E_SI_set, E]; % store the additioanl trial tracking errors
        E_SI_set_norm = [E_SI_set_norm; 0.5 * norm(E)^2]; % store the norm of additioanl trial tracking errors
    end
    E_SI_set_diff = [];
    for i = 1 : m * N % calculate the transposed right-hand side matrix in (7)
        E_temp = E_SI_set(:, i) - E_ini;
        E_SI_set_diff = [E_SI_set_diff, E_temp];
    end
    J = inv(U_SI_set_diff') * E_SI_set_diff'; % solve the linear interpolation problem in (7)
    J = J'; % obtain the estimated Jacobian
    if norm(J' * E_ini) <= ep && Delta <= mu * norm(J' * E_ini) % the criticality phase is not triggered, no more additional trials are required
        break;
    elseif norm(J' * E_ini) > ep % the criticality phase is not triggered, no more additional trials are required
        break;
    elseif norm(J' * E_ini) <= ep && Delta > mu * norm(J' * E_ini) % the criticality phase is triggered and run steps 5-9 in algorithm 1
        Delta = reduce_para * Delta; % reduce the trust-region radius
        E_SI_set_norm_temp_store = [E_SI_set_norm_temp_store; E_SI_set_norm]; % store the additional trial information
    end   
end
Delta_new = Delta; % update the trust-region radius due to the criticality phase 
if size(E_SI_set_norm_temp_store, 1) > 0
    E_SI_set_norm = E_SI_set_norm_temp_store; % if the criticality phase is triggered, store the additional trial information
end
[delta_U, E_delta_U] = trustregprob(J' * J,  -J' * E_ini, Delta, 'false'); % solve the trust-region subproblem
U = U_ini + delta_U; % update the main trial input
E_model_diff = 0.5 * norm(E_ini)^2 - (0.5 * norm(E_ini)^2 + E_delta_U); % calculate the denominator of the ratio in (20)