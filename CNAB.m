%% CNAB Implementation for the advection diffusion equation 



close all
clear all
%
% Define computational grid
    M = 100 ;
    N = 100 ;
    T_final = 2 ; 
    c = 1 ; 
    mu = 0.5 ; 
    dx = 2*pi/M; x_j = 0: dx: 2*pi-dx;
    dt = T_final/N;
    disp(['dx = ', num2str(dx, 4),'  dt = ', num2str(dt, 4)])
    disp(' ')
    
    nplot = ceil(N/250)
    
    
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
    
    
    
% Build the system Matrix

diag1 = (1+sigma)*ones(1,M) ;
up1 = -sigma/2*ones(1,M-1) ;
low1 = -sigma/2*ones(1,M-1) ;
sysMatrix = diag(diag1) + diag(up1, 1) + diag(low1, -1) ;
sysMatrix(1,end) = -sigma/2 ;
sysMatrix(end,1) = -sigma/2 ;
    
    
    
      
    
% Set the initial data 

u0 = zeros(1,M) ; 
u0(1:M/2) = 1 ;
u0 = sin(x_j) ;


 % Initialize
    u = u0;
    uHat_0 = fft(u0)/M ;
    kMode = [0:M/2-1 -M/2:-1];
    omega = c*kMode - 1i*mu*kMode.^2;
        % need first two time steps
            u0 = u;
            uHat = uHat_0.*exp(-omega*1i*dt);
            u1 = M*real(ifft(uHat)); 
            
            
            
% Solve the IMEX CNAB scheme system

for k = 2:N
    
    % Update time
    t = k*dt;
        tgrid(k+1) = t ;
    
    
    
    
    % Construct the right hand side
    rhs = (1-sigma)*u1 + (sigma/2 - 3*lambda/4)*[u1(2:end)  u1(1)] +  (sigma/2 + 3*lambda/4)*[u1(end)  u1(1:end-1)] ...
       + (lambda/4)*[u0(2:end)  u0(1)] - lambda/4*[u0(end)  u0(1:end-1)] ;
   
    % Use the system Matrix to solve for the next time step
    u = sysMatrix/rhs ;
    u = u';
    % Update the previous time steps
    u0 = u1 ;
    u1 = u ;
    
    
    
    if mod(k, nplot) == 0
                plot(x_j,u,'k', 'LineWidth', 2);
                hold on
                uHat = uHat_0.*exp(-omega*1i*t);
                u_exact = M*real(ifft(uHat));
                plot(x_j, u_exact, 'r', 'LineWidth', 2)
                axis([0 2*pi -1 1])
                hold off
                err = norm(u'-u_exact')/sqrt(M);
            axis([0 2*pi -1.5 1.5])
            title(['Time = ', num2str(t),' Error = ', num2str(err)])
            xlabel('x')
            ylabel('u')
            drawnow
            end
   
        
end 
            
            
            
            
            
            
 
            
   
    