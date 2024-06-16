% this function is the implementation of the paper [5], the parameters are tuned to achieve the best performance
function [U_DDOILC, Xi_DDOILC] = DDOILC(j, E_DDOILC, U_DDOILC, Y_DDOILC, U_DDOILC_all, Y_DDOILC_all, Xi_DDOILC_ini, Xi_DDOILC, p, m, N)
rho = 0.12; 
lambda = 0.06;  
eta = 0.4; 
mu = 2; 
e_min = 0.0001;
e_max = 1000;
if j > 1
    Y_DDOILC_new = Y_DDOILC;
    Y_DDOILC = Y_DDOILC_all(:, j - 1)';
    U_DDOILC_new = U_DDOILC;
    U_DDOILC = U_DDOILC_all(:, j - 1)';
    Xi_DDOILC_new = Xi_DDOILC;
    for i = 1 : (p * N) % update the Jacobian matrix
    Xi_DDOILC_new(i, 1 : i) = Xi_DDOILC(i, 1 : i) ...
        + ( eta * ( Y_DDOILC_new(i) - Y_DDOILC(i) - Xi_DDOILC(i, 1 : i) * ( U_DDOILC_new(1 : i) - U_DDOILC(1 : i) )' )...
        / ( mu + norm( U_DDOILC_new(1 : i) - U_DDOILC(1 : i) )^2 ) ) ...
        * ( U_DDOILC_new(1 : i) - U_DDOILC(1 : i) );
        for ii = 1 : i % detect whether the each element sign is changed, or each element is too large/small
            if sign( Xi_DDOILC_new(i, ii) ) ~= sign( Xi_DDOILC_ini(i, ii) ) || norm(Xi_DDOILC_new(i, ii)) <= e_min || norm(Xi_DDOILC_new(i, ii)) >= e_max
                Xi_DDOILC_new(i, 1 : i) = Xi_DDOILC_ini(i, 1 : i); % if they are detected, reset the element to the initial value
                break;
            end
        end
    end
Xi_DDOILC = Xi_DDOILC_new; % store the updated the Jacobian matrix
end
U_DDOILC_new = reshape(U_DDOILC, [], 1) + (rho / ( lambda + norm(Xi_DDOILC)^2 ) ) * Xi_DDOILC' * reshape(E_DDOILC, [], 1); % update the input trajectory
U_DDOILC_new = reshape(U_DDOILC_new, m, N);     
U_DDOILC = U_DDOILC_new;