function xdot = cstr_discrete(x, u)
%% parameter
Ts = 0.1; 
theta = 1;
beta = 0.3;
gamma = 20;
B = 1;
Da = 0.072;
%% dynamics
xdot(1,1) = (1 - Ts * theta) * x(1,1) +...
    Ts * Da * (1 - x(1,1)) * exp(x(2,1) / (1 + (x(2,1) / gamma)));
xdot(2,1) = (1 - Ts * theta) * x(2,1) +...
    Ts * B * Da * (1 - x(1,1)) * exp(x(2,1) / (1 + (x(2,1) / gamma))) -...
    Ts * beta * x(2,1) + Ts * beta * u;