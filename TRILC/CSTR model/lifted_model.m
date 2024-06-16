function Y = lifted_model(U, N)
x_ini_1 = 0.57;
x_ini_2 = 0.3;
x = [x_ini_1; x_ini_2];
Y = [];  
u_temp = U;
for k = 0 : N - 1
    u = u_temp(k + 1);
    xdot = cstr_discrete(x, u);
    x_1 = xdot(1, 1);
    x_2 = xdot(2, 1);
    x = xdot;
    Y = [Y; x_2];
end