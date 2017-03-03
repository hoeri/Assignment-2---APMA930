%% Driver for advectDiffuse.m
% Use the "Run Section" button in the editor to run these examples 
% individually.
%
%%
% Fully Explicit - Small cell Reynolds number, no decay of high frequency modes
    clc
    mu = 1; disp(['mu = ', num2str(mu)])
    M = 128;
    dx = 2*pi/M;
    c = 1;, disp(['c = ', num2str(c)])
    dt = min(mu*dx/c^2, 0.5*dx^2/mu)
    T_final = 2*pi
    N = ceil(T_final/dt)
    u0 = @(x) 2*(.5 - heaviside(x-pi));
    u_exact = @(x,t) exp(-mu*t)*sin(x - c*t);
    advectDiffuse(c, mu, T_final, M, N, u0, 1, 't')
%%
% Fully Explicit - Cell Reynolds Number = 1, step function
    clc
    M = 256;
    dx = 2*pi/M;
    c = 1;, disp(['c = ', num2str(c)])
    mu = 0.5*c*dx; disp(['mu = ', num2str(mu)])
    dt = dx;
    T_final = 2*pi
    N = ceil(T_final/dt)
    u_exact = @(x,t) exp(-mu*t)*sin(x - c*t);
    u0 = @(x) 2*(.5 - heaviside(x-pi));
    u_exact = @(x,t) exp(-mu*t)*sin(x - c*t);
    advectDiffuse(c, mu, T_final, M, N, u0, 1, 't')
%%
% IMEX Euler - Small cell Reynolds number, this produces a nice solution
    clc
    mu = 1; disp(['mu = ', num2str(mu)])
    M = 256;
    dx = 2*pi/M;
    c = 1;, disp(['c = ', num2str(c)])
    dt = min(mu*dx/c^2, 0.5*dx^2/mu)
    T_final = 2*pi
    N = ceil(T_final/dt)
    u0 = @(x) 2*(.5 - heaviside(x-pi));
    u_exact = @(x,t) exp(-mu*t)*sin(x - c*t);
    advectDiffuse(c, mu, T_final, M, N, u0, 2, 't')
%%
% IMEX Euler - Cell Reynolds Number = 1, step function -  a number of modes
% don't damp out
    clc
    M = 256;
    dx = 2*pi/M;
    c = 1;, disp(['c = ', num2str(c)])
    mu = 0.5*c*dx; disp(['mu = ', num2str(mu)])
    dt = 0.9*dx;
    T_final = 2*pi
    N = ceil(T_final/dt)
    u_exact = @(x,t) exp(-mu*t)*sin(x - c*t);
    u0 = @(x) 2*(.5 - heaviside(x-pi));
    u_exact = @(x,t) exp(-mu*t)*sin(x - c*t);
    advectDiffuse(c, mu, T_final, M, N, u0, 2, 't')
%%
% SBDF  - Cell Reynolds Number = 1, step function - nice solution
    clc
    M = 256;
    dx = 2*pi/M;
    c = 1;, disp(['c = ', num2str(c)])
    mu = 0.5*c*dx; disp(['mu = ', num2str(mu)])
    dt = 0.5*dx;
    T_final = 2*pi
    N = ceil(T_final/dt)
    u_exact = @(x,t) exp(-mu*t)*sin(x - c*t);
    u0 = @(x) 2*(.5 - heaviside(x-pi));
    u_exact = @(x,t) exp(-mu*t)*sin(x - c*t);
    advectDiffuse(c, mu, T_final, M, N, u0, 3, 't')
%%
% CNLF  - Cell Reynolds Number = 1, step function, high frequency modes
% don't Decay
    clc
    M = 256;
    dx = 2*pi/M;
    c = 1;, disp(['c = ', num2str(c)])
    mu = 0.5*c*dx; disp(['mu = ', num2str(mu)])
    dt = dx;
    T_final = 2*pi
    N = ceil(T_final/dt)
    u_exact = @(x,t) exp(-mu*t)*sin(x - c*t);
    u0 = @(x) 2*(.5 - heaviside(x-pi));
    u_exact = @(x,t) exp(-mu*t)*sin(x - c*t);
    advectDiffuse(c, mu, T_final, M, N, u0, 4, 't')
%%

% CNAB - Cell Raynolds Number = 1,

  clc
    M = 256;
    dx = 2*pi/M;
    c = 1;, disp(['c = ', num2str(c)])
    mu = 0.5*c*dx; disp(['mu = ', num2str(mu)])
    dt = 0.5*dx; 
    T_final = 2*pi
    N = ceil(T_final/dt)
    u_exact = @(x,t) exp(-mu*t)*sin(x - c*t);
    u0 = @(x) 2*(.5 - heaviside(x-pi));
    u_exact = @(x,t) exp(-mu*t)*sin(x - c*t);
    advectDiffuse(c, mu, T_final, M, N, u0, 5, 't')




