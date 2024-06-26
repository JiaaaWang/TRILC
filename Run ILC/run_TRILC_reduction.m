%% Algorithm parameters
Delta = 10; % trust region initial radius
Delta_max = 30; % maxiaml trust region radius
eta_1 = 0.01; % acceptence threshold
eta_2 = 0.9; % acceptence threshold
w_inc = 1.5; % trust region increase para
w_inc_max = 5; % trust region increase para
w_dec = 0.9; % trust region decrease para
ep = 0.01; % gradient threshold 
mu = 1; % gradient scaling 
tracking_error_th = 0.01; % accuracy tolerance
%% Algorithm initialization
U_all = []; % input storage, each column is one trial trajectory
Y_all = []; % output storage, each column is one trial trajectory
E_all = []; % tracking error storage, each column is one trial trajectory
U = U_ini_ini; % input trajectory of the first trial
E_SI_set_norm_all_reduction = []; % all tracking errors (main + additional trials)
%% Start experiment
for j = 1 : M
    Y = lifted_model(U, N); % calculate the output trajectory of each main trial 
    E = R - Y; % generate the main trial tracking error
    if 0.5 * norm(E)^2 <= tracking_error_th % stopping criteria
        E_SI_set_norm_all_reduction = [E_SI_set_norm_all_reduction; 0.5 * norm(E)^2]; % store the tracking errors (main + additioanl trials)
        break;
    end
    U_all = [U_all, reshape(U, [], 1)]; % store the main trial input trajectory 
    Y_all = [Y_all, reshape(Y, [], 1)]; % store the main trial output trajectory 
    E_all = [E_all, reshape(E, [], 1)]; % store the main trial tracking error trajectory 
    if j > 1
        ratio = ( 0.5 * norm(E_all(:, j - 1))^2 - 0.5 * norm(E_all(:, j))^2 ) / E_model_diff; % calculate the ratio in (20)
        if ratio >= eta_2 % successful input increment
            Delta = min(w_inc * Delta, Delta_max); % enlarge the trust-region radius
        elseif ratio < eta_2 && ratio >= eta_1 % acceptable input increment
            Delta = w_dec * Delta; % reduce the trust-region radius
        elseif ratio < eta_1 % unacceptable input increment
            Delta = w_dec * Delta; % reduce the trust-region radius
            U = U_all(:, j - 1); % discard the input increment
        end  
    end
    U_ini = U; % update the main trial input
    E_ini = E; % update the main trial tracking error 
    if j == 1 % do the same thing as the TRILC withou trial reduction in the first main trial
        U_SI_set = zeros(m * N, m * N);
        E_SI_set = [];
    else  
        E_SI_set = E_SI_set_new; % input the last main trial data set to the function after the first main trial
    end
    [U, E_model_diff, E_SI_set_norm, Delta_new, delta_U, U_SI_set_new, E_SI_set_new]...
        = TRILC_reduction(R, U_ini, E_ini, Delta, m, N, ep, mu, U_SI_set, j, E_SI_set); % run TRILC (with trial reduction) at this main trial, feed previous additional input and tracking error sets to the function
    E_SI_set_norm_all_reduction = [E_SI_set_norm_all_reduction; [0.5 * norm(E_ini)^2; E_SI_set_norm]]; % store the additional trial tracking errors
    Delta = Delta_new; % update the trust-region radius due to the criticality phase 
    U_SI_set = U_SI_set_new;  
end