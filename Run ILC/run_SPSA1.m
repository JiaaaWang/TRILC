%% Algorithm initialization
U_SPSA1_all = [];
Y_SPSA1_all = [];
E_SPSA1_all = [];
U_SPSA1 = U_ini_ini;  
%% Start experiment
for j = 1 : (size(E_SI_set_norm_all, 1)) / 2 % SPSA2ILC is terminated with the trial length of TRILC (withou trial reduction), note that each main trial is equipped one additional trial
    Y_SPSA1 = lifted_model(U_SPSA1, N); % calculate the output trajectory of each main trial 
    E_SPSA1 = R - Y_SPSA1; % generate the main trial tracking error
    U_SPSA1_all = [U_SPSA1_all, reshape(U_SPSA1, [], 1)]; % store the main trial input trajectory
    E_SPSA1_all = [E_SPSA1_all; 0.5 * norm( E_SPSA1 )^2]; % store the norm of the main trial tracking error trajectory 
    U_ini = U_SPSA1;
    [U_SPSA1, E_SPSA1_set_norm] = SPSA1ILC(R, U_ini, m, N); % run SPSA1ILC at this trial
    E_SPSA1_all = [E_SPSA1_all; E_SPSA1_set_norm]; % store the norm of the each trial tracking error
end