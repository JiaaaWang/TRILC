%% Algorithm parameters
Delta = 10; % trust region initial radius
Delta_max = 30; % maxiaml trust region radius
eta_1 = 0.01; % acceptence threshold
eta_2 = 0.9; % acceptence threshold
w_inc = 1.5; % trust region increase para
w_inc_max = 5; % trust region increase para
w_dec = 0.5; % trust region decrease para
ep = 0.01; % gradient threshold
mu = 1; % gradient scaling
tracking_error_th = 0.01; % accuracy tolerance
%% Remark
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% The parameters in Algorithm 1 (TRILC) can be adjusted in a range, which provide 
% some design freedom to manage the algorithm behaviour and improve the control performance. 
% In fact, different selections for the parameters of trust-region methods are discussed 
% comprehensively in Yuan (2000), Gould et al. (2005) and Chapter 17.1 of Conn et al. (2000), 
% which provide guidelines and foundations for applying Algorithm 1 to other nonlinear systems. 
% It is worth mentioning that the parameters eta_1, eta_2, w_dec and w_inc are important for
% managing the trust-region radius Delta_j, we discuss their selections from engineering point
% of view as follows. If r_j>=eta_2, the trust-region radius is enlarged from Delta_j to w_inc*Delta_j. 
% In practice, the nonlinearity around different inputs may vary greatly, therefore, it is safe 
% to not enlarge the trust-region radius easily. The radius should be enlarged only if the 
% data-driven representation produces a similar tracking error descent as the true nonlinear 
% system (describes the local nonlinearity well), which implies that eta_2 should be close to 1. 
% We choose eta_2=0.9 here. The radius should also not be enlarged too much, hence, we use w_inc=1.5, 
% instead of 2 or 2.5 for the general optimization problems (Cartis & Roberts, 2019; Conn et al., 2000). 
% On the other hand, the input increment will be discarded if r_j<eta_1. Even if the updated input 
% cannot improve the performance sufficiently, the tracking error is still possibly reduced, 
% which can be meaningful in practice, therefore, we choose eta_1=0.01 here, instead of the common used
% values 0.05 or 0.1. Note that the small value of eta_1 is also advocated in Conn et al. (2009). 
% If r_j<eta_2, the trust-region radius is reduced from Delta_j to w_dec*Delta_j, we select 
% w_dec=0.5 in TRILC (without trial reduction), which is consistent with the selection in Cartis & Roberts (2019).

% Reference:
% Yuan, Y. (2000). A review of trust region algorithms for optimization. In Proceedingss of the Fourth International Cogress on Industrial and Applied Mathematics (p. 271每282). Oxford, UK.
% Conn, A. R., Gould, N. I. M., & Toint, P. L. (2000). Trust region methods. Philadelphia: MPS/SIAM.
% Gould, N. I. M., Orban, D., Sartenaer, A., & Toint, P. L. (2005). Sensitivity of trust-region algorithms to their parameters. 4OR, 3, 227每241.
% Cartis, C., & Roberts, L. (2019). A derivative-free Gauss每Newton method. Mathematical Programming Computation, 11, 631每674.
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%% Algorithm initialization
U_all = []; % input storage, each column is one trial trajectory
Y_all = []; % output storage, each column is one trial trajectory
E_all = []; % tracking error storage, each column is one trial trajectory
U = U_ini_ini; % input trajectory of the first trial
E_SI_set_norm_all = []; % all tracking errors (main + additional trials)
%% Start experiment
for j = 1 : M
    Y = lifted_model(U, N); % calculate the output trajectory of each main trial 
    E = R - Y; % generate the main trial tracking error
    if 0.5 * norm(E)^2 <= tracking_error_th % stopping criteria
        E_SI_set_norm_all = [E_SI_set_norm_all; 0.5 * norm(E)^2]; % store the tracking errors (main + additioanl trials)
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
            Delta = w_dec * Delta;  % reduce the trust-region radius
            U = U_all(:, j - 1); % discard the input increment
        end       
    end  
    U_ini = U; % update the main trial input
    E_ini = E; % update the main trial tracking error 
    [U, E_model_diff, E_SI_set_norm, Delta_new] = TRILC(R, U_ini, E_ini, Delta, m, N, ep, mu); % run TRILC at this main trial
    E_SI_set_norm_all = [E_SI_set_norm_all; [0.5 * norm(E_ini)^2; E_SI_set_norm]]; % store the additional trial tracking errors
    Delta = Delta_new; % update the trust-region radius due to the criticality phase 
end