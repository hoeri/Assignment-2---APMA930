%% Lax Wendroff Method for solving the advection diffusion equation 


close all ; clear ;

% Wave speed 
c = 1;



% Define computational grid for M = 64 and lambda = 0.5
    M = 128    ;   dx = 2*pi/M ;    x_j = 0:dx:2*pi-dx ; 
    T = 2*pi ;  N = 128 ;   dt = T/N  ;
    disp(['\Delta x = ', num2str(dx, 4),'  \Delta t = ', num2str(dt, 4)])
    lambda = c * dt / dx;
    disp(['\lambda = ', num2str(lambda, 4)])
    
% Define computational grid for M = 64 and lambda = 1
%     M = 128    ;   dx = 2*pi/M ;    x_j = 0:dx:2*pi-dx ; 
%     T = 4*pi ;  N = 256 ;   dt = T/N  ;
%     disp(['\Delta x = ', num2str(dx, 4),'  \Delta t = ', num2str(dt, 4)])
%     lambda = c * dt / dx;
%     disp(['\lambda = ', num2str(lambda, 4)])  

% Define computational grid for M = 64 and lambda = 1/sqrt(2)
%     M = 128    ;   dx = 2*pi/M ;    x_j = 0:dx:2*pi-dx ; 
%     T = 4*pi ;  N = 256*sqrt(2) ;   dt = T/N  ;
%     disp(['\Delta x = ', num2str(dx, 4),'  \Delta t = ', num2str(dt, 4)])
%     lambda = c * dt / dx;
%     disp(['\lambda = ', num2str(lambda, 4)]) 



    
% Initial condtions
    f0 = @(x) exp(-4*(x-pi).^2);
    IVs = f0(x_j);
     IVs = [zeros(1,M/2) ones(1,M/2)] ;
    
    ymax = max(IVs);
    ymin = min(IVs);
    
    
    
 % Initnialize u   
    u = IVs;
    
    
 
    
for k = 1:N
    
    t = k*dt ;
    
    u = u  + lambda/2*([u(end) u(1:end-1)] - [u(1:end-1) u(1)]) + lambda^2/2 *([u(1:end-1) u(1)] - 2*u + [u(end) u(1:end-1)] )  ;
    
   
    
    
    
end


hold off
            plot(x_j, IVs, 'r', 'LineWidth', 2)
            hold on
            plot(x_j, u, 'k', 'LineWidth', 2);  
            axis([0 2*pi ymin - 0.1  ymax + 0.1])
            xlabel('x');  ylabel('u')
            title(['Lax Wendroff, Time =' num2str(t) ', Lambda =' num2str(lambda)])
        pause(0.0001)