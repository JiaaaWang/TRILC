m = 1; % input number
p = 1; % output number
N = 100; % data length in one trial
M = 5; % Maximal main trial number
U_ini_ini = 5 * ones(N, m); % input trajectory of the first trial
Y_ini = zeros(N, p); % output trajectory of the first trial, initialize from 0
R = 1.96 * ones(p * N, 1); % reference trajectory