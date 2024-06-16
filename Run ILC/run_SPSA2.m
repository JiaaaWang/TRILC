%% Algorithm initialization
U_SPSA2_all = [];
Y_SPSA2_all = [];
E_SPSA2_all = [];
U_SPSA2 = U_ini_ini;  
%% Start experiment
for j = 1 : (size(E_SI_set_norm_all, 1)) / 3 % SPSA2ILC is terminated with the trial length of TRILC (withou trial reduction), note that each main trial is equipped two additional trials
    Y_SPSA2 = lifted_model(U_SPSA2, N); % calculate the output trajectory of each main trial 
    E_SPSA2 = R - Y_SPSA2; % generate the main trial tracking error
    U_SPSA2_all = [U_SPSA2_all, reshape(U_SPSA2, [], 1)]; % store the main trial input trajectory
    E_SPSA2_all = [E_SPSA2_all; 0.5 * norm( E_SPSA2 )^2]; % store the norm of the main trial tracking error trajectory 
    U_ini = U_SPSA2;
    [U_SPSA2, E_SPSA2_set_norm] = SPSA2ILC(j, R, U_ini, m, N); % run SPSA2ILC at this trial
    E_SPSA2_all = [E_SPSA2_all; E_SPSA2_set_norm]; % store the norm of the each trial tracking error
end