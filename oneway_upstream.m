%% One Way Wave Equation
% Solve u_t + c u_x = 0, periodic on [0, 2 pi] using
% forward differences in time, upstream weighting in space.
%
    close all; clear;
    
    movie = 1;
%
% Wave speed
    c = 1;
%


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

    %IVs  = exp(-4*(x_j-pi).^2);
    IVs = sin(5*x_j);
    IVs = [zeros(1,32) ones(1,32)] ;
    ymax = max(IVs);
    ymin = min(IVs);

%
% Plot initial information

    u = IVs;  
    plot(x_j,u,'r');  hold on;  plot(x_j,u,'kx');
    axis([0 2*pi 0 1.25]);
    xlabel('x');  ylabel('u')
    title('One-Way Wave, Upstream Weighting')


    for k=1 : 5*N/5
        t = k*dt;
        u = u - lambda*(u - [u(end) u(1:end-1)]) ;
       
        
           
    end
    
      hold off
            plot(x_j, IVs, 'r', 'LineWidth', 2)
            hold on
            plot(x_j, u, 'k', 'LineWidth', 2);  
            axis([0 2*pi ymin - 0.1  ymax + 0.1])
            xlabel('x');  ylabel('u')
            title(['One-Way Wave, Upstream Weighting, Time =' num2str(t) ', Lambda =' num2str(lambda)])
        pause(0.0001)

        
        
    
    