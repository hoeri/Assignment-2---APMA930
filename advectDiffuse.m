%% Advection Diffusion Equation
% Solve the advection diffusion equation
%   u_t + c u_x = mu u_xx, 0 < x < 2 pi 
% with periodic boundary conditions via a variety of time stepping methods,
% all derivatives centred in space.
%
% advectDiffuse(c, mu, T_final, M, N, init_cond, imethod)
%
% Input variables:
%   c, mu = advection speed, viscosity
%   T_final = final time
%   M = number of grid points in x 
%   N = number of time steps
%   init_cond = function defining u(x,0)
%   imethod =
%       1:  forward Euler
%       2:  imex Euler
%       3:  imex SBDF
%       4:  imex CNLF
%   debug = 't' if there is an exact solution
%   u_exact = function defining exact solution u(x,t)
%
function advectDiffuse(c, mu, T_final, M, N, init_cond, imethod, ...
                       debug)

    close all
    
    % Store c and mu for later
   
%
% Define computational grid
    dx = 2*pi/M; x_j = 0: dx: 2*pi-dx;
    dt = T_final/N;
    disp(['dx = ', num2str(dx, 4),'  dt = ', num2str(dt, 4)])
    disp(' ')
    
% 
% Calculate numerical parameters
    lambda = c * dt / dx;
    disp(['lambda = ', num2str(lambda, 4)])
    sigma = mu * dt / dx^2;
    disp(['sigma = ', num2str(sigma, 4)])
    disp(['Cell Reynolds number = ', num2str(0.5*lambda/sigma)])
    disp(' ')
    disp(['lambda^2 = ', num2str(lambda^2),' 2 sigma = ', num2str(2*sigma)])
    disp(['sigma - lambda/2 = ', num2str(sigma - lambda/2)])
    disp(' ')
    
%
% Plot stability region and location within stability region
    alpha_dt = dt*2*mu*(cos(dx*(1:M/2))-1)/dx^2;
    beta_dt = -dt*c*sin(dx*(1:M/2))/dx;
    figure(1)
    subplot(1, 2, 1)
    stabilityAdvDiff (imethod, alpha_dt, beta_dt)
    
%
% Get Fourier components of initial condition for exact solution
    uHat_0 = fft(init_cond(x_j))/M;
    kMode = [0:M/2-1 -M/2:-1];
    omega = c*kMode - 1i*mu*kMode.^2;
    
%
% Plotting
    nplot = ceil(N/250);
%    nplot = ceil(N/N)
    
%
% Allocate space for solution array
    usave = zeros(N+1, M);
    tgrid = zeros(N+1, 1);

%
% Set up system matrix according to imethod.
    switch imethod
        case 1   % forward Euler
            diagonal = (1-2*sigma)*ones(1, M);
            up = (sigma - 0.5*lambda)*ones(1, M-1);
            low = (sigma + 0.5*lambda)*ones(1, M-1);
            sysMatrix = diag(diagonal) + diag(up, 1) + diag(low, -1);
            sysMatrix(1, end) = sigma + 0.5*lambda;
            sysMatrix(end, 1) = sigma - 0.5*lambda;
        case 2   % imex Euler
            diagonal = (1+2*sigma)*ones(1, M);
            up = -sigma*ones(1, M-1);
            low = -sigma*ones(1, M-1);
            sysMatrix = diag(diagonal) + diag(up, 1) + diag(low, -1);
            sysMatrix(1, end) = -sigma;
            sysMatrix(end, 1) = -sigma;
        case 3   % imex SBDF
          %
          % System matrix for first time step
            diag1 = (1+2*sigma)*ones(1, M);
            up1 = -sigma*ones(1, M-1);
            low1 = -sigma*ones(1, M-1);
            sysMatrix1 = diag(diag1) + diag(up1, 1) + diag(low1, -1);
            sysMatrix1(1, end) = -sigma;
            sysMatrix1(end, 1) = -sigma;
          %
          % System matrix for remaining time steps
            diagonal = (3+4*sigma)*ones(1, M);
            up = -2*sigma*ones(1, M-1);
            low = -2*sigma*ones(1, M-1);
            sysMatrix = diag(diagonal) + diag(up, 1) + diag(low, -1);
            sysMatrix(1, end) = -2*sigma;
            sysMatrix(end, 1) = -2*sigma;
        case 4   % imex CNLF
          %
          % System matrix for first time step
            diag1 = (1+2*sigma)*ones(1, M);
            up1 = -sigma*ones(1, M-1);
            low1 = -sigma*ones(1, M-1);
            sysMatrix1 = diag(diag1) + diag(up1, 1) + diag(low1, -1);
            sysMatrix1(1, end) = -sigma;
            sysMatrix1(end, 1) = -sigma;
          %
          % System matrix for remaining time steps
            diagonal = (1+2*sigma)*ones(1, M);
            up = -sigma*ones(1, M-1);
            low = -sigma*ones(1, M-1);
            sysMatrix = diag(diagonal) + diag(up, 1) + diag(low, -1);
            sysMatrix(1, end) = -sigma;
            sysMatrix(end, 1) = -sigma;
            
            
        case 5 % Imex CNAB   
            
            diag1 = (1+sigma)*ones(1,M) ;
            up1 = -sigma/2*ones(1,M-1) ;
            low1 = -sigma/2*ones(1,M-1) ;
            sysMatrix = diag(diag1) + diag(up1, 1) + diag(low1, -1) ;
            sysMatrix(1,end) = -sigma/2 ;
            sysMatrix(end,1) = -sigma/2 ;
    
            
            
    end

%
% Initialize
    u = init_cond(x_j)';
    switch imethod
        case {1, 2}
            usave(1,:) = u';
            tgrid(1) = 0;
            kstart = 1;
        case {3, 4}     % need first two time steps
            u0 = u;
            uHat = uHat_0.*exp(-omega*1i*dt);
            u1 = M*real(ifft(uHat)); u1 = u1';
            usave(1:2,:) = [u0'; u1'];
            tgrid(1:2) = [0;dt];
            kstart = 2;
            
        case 5
            u0 = u;
            uHat = uHat_0.*exp(-omega*1i*dt);
            u1 = M*real(ifft(uHat)); u1 = u1';
            usave(1:2,:) = [u0'; u1'];
            tgrid(1:2) = [0;dt];
            kstart = 2;
            u_exact = @(x,t) exp(-mu*t)*sin(x - c*t);
    end
            
%
% Plot initial condition
    figure(1)
    subplot(1, 2, 2)
    disp('Adjust figure, hit return')
    pause
    shg
    plot(x_j, init_cond(x_j), 'r')
   % hold on
    xlabel('x')
    ylabel('y')
    title('Advection Diffusion')

    for k = kstart : N
        t = k*dt;
        tgrid(k+1) = t;
        
        switch imethod
%
%  solve for next timestep
            case 1   % forward Euler
                rhs = u;
                u = sysMatrix*rhs;
            case 2   % imex Euler
                rhs = u + 0.5*lambda*([u(end); u(1:end-1)] ...
                                       - [u(2:end); u(1)]);
                u = sysMatrix\rhs;
            case 3   % imex SBDF
                rhs = 4*u1 - u0  + 2*lambda*([u1(end); u1(1:end-1)] ...
                                           - [u1(2:end); u1(1)]) ...
                                 - lambda*([u0(end); u0(1:end-1)] ...
                                         - [u0(2:end); u0(1)]); 
                u = sysMatrix\rhs;
                u0 = u1;
                u1 = u;
            case 4   % imex CNLF
                rhs = u0  + lambda*([u1(end); u1(1:end-1)] ...
                                  - [u1(2:end); u1(1)]);
                u = sysMatrix\rhs;
                u0 = u1;
                u1 = u;
                
            case 5 % Imex CNAB
                rhs = (1-sigma)*u1  + (sigma/2 - 3*lambda/4)*[u1(2:end) ; u1(1)] +  (sigma/2 + 3*lambda/4)*[u1(end) ; u1(1:end-1)] ...
                + (lambda/4)*[u0(2:end) ; u0(1)] - lambda/4*[u0(end) ; u0(1:end-1)] ;
                
                u = sysMatrix\rhs ;
                u0 = u1;
                u1 = u;
            
        end
        usave(k+1, :) = u';
        
%   
% interval plotting
            if mod(k, nplot) == 0
                plot(x_j,u,'k', 'LineWidth', 2);
                hold on
                uHat = uHat_0.*exp(-omega*1i*t);
                u_exact = M*real(ifft(uHat));
                plot(x_j, u_exact, 'r', 'LineWidth', 2)
                axis([0 2*pi -1 1])
                hold off
                err = norm(u'-u_exact)/sqrt(M);
            axis([0 2*pi -1.5 1.5])
            title(['Time = ', num2str(t), ' Error = ', num2str(err)])
            xlabel('x')
            ylabel('u')
            drawnow
            end
%        end
    end
    if debug == 't'
        err = norm(u'-u_exact)/sqrt(M)
    end
end

% Stability for Advection Diffusion

function stabilityAdvDiff (imethod, alphaPoint, betaPoint)

    [alpha_dt, beta_dt] = meshgrid([-5:.05:5], [-5:.05:5]);
    
    switch imethod
        case 1 % forward Euler
            A = 1 + alpha_dt + 1i*beta_dt;
        case 2 % imex Euler
            A = (1 + 1i*beta_dt)./(1 - alpha_dt);
        case 3 % spdf
            a = 3-2*alpha_dt;
            b = -4-4*1i*beta_dt;
            c = 1+2*1i*beta_dt;
            xi_1 = -b./(2*a) + (1./(2*a)).*sqrt(b.^2 - 4*a.*c);
            xi_2 = -b./(2*a) - (1./(2*a)).*sqrt(b.^2 - 4*a.*c);
            A = max(abs(xi_1),abs(xi_2));
        case 4 % cnlf
            a = 1-alpha_dt;
            b = -2*1i*beta_dt;
            c = -1-alpha_dt;
            xi_1 = -b./(2*a) + (1./(2*a)).*sqrt(b.^2 - 4*a.*c);
            xi_2 = -b./(2*a) - (1./(2*a)).*sqrt(b.^2 - 4*a.*c);
            A = max(abs(xi_1),abs(xi_2));
            
        case 5 % CNAB
            M = 256 ;
            dx = 0.5*pi/M
            c = 1 ;
            mu = 0.5*c*dx;
             tau1 = (1/2*1i)*(-1i*alpha_dt*mu+3*beta_dt*c-2*1i+sqrt(-mu^2*alpha_dt.^2-(10*1i)*c*mu*alpha_dt*beta_dt+9*c^2*beta_dt.^2-4*alpha_dt*mu-(4*1i)*beta_dt*c-4))./(mu*alpha_dt-2) ;
            tau2 = (1/2)*(-(3*1i)*beta_dt*mu-alpha_dt*c-2+sqrt(-9*c^2*beta_dt.^2+(10*1i)*c*mu*alpha_dt.*beta_dt+alpha_dt.^2*mu^2+(4*1i)*beta_dt*c+4*alpha_dt*mu+4))./(mu*alpha_dt-2) ;
            A = abs(tau2) ;
            
    end
    pcolor(alpha_dt, beta_dt, abs(A))
    shading flat 
    colormap(jet)
    colorbar
    hold on
    contour (alpha_dt, beta_dt, abs(A), [1 1],'k') % plot c=1 contour line
    plot (alphaPoint, betaPoint, 'k*')
    xlabel('\alpha \Delta t')
    ylabel('\beta \Delta t')
    caxis([0 2])
    switch imethod
        case 1
            axis([-3 0 -2 0])
            title('Stability - Forward Euler')
        case 2
            title('Stability - Backwards Euler')
        case 3
            title('Stability - SBDF')
        case 4 
            title('Stability - CNLF')
        case 5 
            title('Stability - CNAB')
    end
end
