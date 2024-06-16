%% Algorithm initialization
U_DDOILC_all = [];
Y_DDOILC_all = [];
E_DDOILC_all = [];
U_DDOILC = U_ini_ini;  
%% Start experiment
for j = 1 : size(E_SI_set_norm_all, 1) % DDOILC is terminated with the trial length of TRILC (withou trial reduction)
    Y_DDOILC = lifted_model(U_DDOILC, N); % calculate the output trajectory of each main trial 
    E_DDOILC = R - Y_DDOILC; % generate the main trial tracking error
    U_DDOILC_all = [U_DDOILC_all, reshape(U_DDOILC, [], 1)]; % store the main trial input trajectory 
    Y_DDOILC_all = [Y_DDOILC_all, reshape(Y_DDOILC, [], 1)]; % store the main trial output trajectory 
    E_DDOILC_all = [E_DDOILC_all, reshape(E_DDOILC, [], 1)]; % store the main trial tracking error trajectory
    if j == 1 
        Xi_DDOILC = 0.1 * ones(p * N, m * N); % initialize the estimated Jacobian with the initial condition 0.1
        for i = 2 : size(Xi_DDOILC, 2) 
            Xi_DDOILC(1 : (i - 1), i) = zeros(i - 1, 1); % generate a low triangular matrix due to the causality
        end
        Xi_DDOILC_ini = Xi_DDOILC; % input the initial estimated Jacobin matrix to DDOILC function (as the reset condition)
    end
    [U_DDOILC, Xi_DDOILC] = DDOILC(j, E_DDOILC, U_DDOILC, Y_DDOILC, U_DDOILC_all, Y_DDOILC_all, Xi_DDOILC_ini, Xi_DDOILC, p, m, N); % run DDOILC at this trial
end
E_DDOILC_norm = zeros(size(E_SI_set_norm_all, 1), 1);
for i = 1 : size(E_SI_set_norm_all, 1)
    E_DDOILC_norm(i) = 0.5 * norm( E_DDOILC_all(:, i) )^2; % store the norm of the each trial tracking error
end