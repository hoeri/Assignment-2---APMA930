%% One Way Wave Equation - LF
% Solve u_t + c u_x = 0, periodic on [0, 2 pi] using
% Leap Frog time stepping, centred in space.
%
    close all; clear;
%
% Wave speed
    c = 1;


% Define computational grid for M = 64 and lambda = 0.5
    M = 64    ;   dx = 2*pi/M ;    x_j = 0:dx:2*pi-dx ; 
    T = 4*pi ;  N = 256 ;   dt = T/N  ;
    disp(['\Delta x = ', num2str(dx, 4),'  \Delta t = ', num2str(dt, 4)])
    lambda = c * dt / dx;
    disp(['\lambda = ', num2str(lambda, 4)])
    
% Define computational grid for M = 64 and lambda = 1
%     M = 64    ;   dx = 2*pi/M ;    x_j = 0:dx:2*pi-dx ; 
%     T = 4*pi ;  N = 128 ;   dt = T/N  ;
%     disp(['\Delta x = ', num2str(dx, 4),'  \Delta t = ', num2str(dt, 4)])
%     lambda = c * dt / dx;
%     disp(['\lambda = ', num2str(lambda, 4)])    
   
    
    
     

%
% Initial condtions
    f0 = @(x) exp(-4*(x-pi).^2);
    u0 = f0(x_j);
    
    % Initiali conditions as step function
    u0 = [zeros(1,M/2) ones(1,M/2)];
    IC = u0 ;
    fHat = fft(u0);
    
% first time step
    IV1_eul  = u0 - 0.5*lambda*([u0(2:end) u0(1)] - [u0(end) u0(1:end-1)]) ;
    IV1_exact = f0(x_j - c*dt);
    IV1_wrong = f0(x_j - 2*c*dt); % wrong first time step, 
                                         % useful for triggering comp. mode
    u1 = IV1_eul;
    

%

    plot(x_j,u0,'r');  hold on;  plot(x_j,u0,'kx');
    axis([0 2*pi 0 1.25]);
    xlabel('x');  ylabel('u')
    title('Crank-Nicolson, Leap-Frog')

% Plot initial information

    for k = 2 : N
        t = k*dt;
        u = u0 - lambda*([u1(2:end) u1(1)] - [u1(end) u1(1:end-1)]) ;
          
        
        
      
        
      u0 = u1;
      u1 = u;
    end

    hold off
            plot(x_j, IC, 'r', 'LineWidth', 2)
            hold on
            plot(x_j, u, 'k', 'LineWidth', 2);  
            axis([0 2*pi -.2 1.2])
            xlabel('x');  ylabel('u')
            title(['One-Way Wave, CNLF: t = ', num2str(t, '%5.3f') ' and lambda = ' num2str(lambda)])
            shg
             pause(1.d-8)

