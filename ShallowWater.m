%% Linearized shallow water equaiotns in one space dimension

clear all;
close all;

% Pick the parameter G
G = 2;



% Define computational grid for M = 64 and lambda = 1
%     M = 128    ;   dx = 2*pi/M ;    x_j = 0:dx:2*pi-dx ; 
%     T = 4*pi ;  N = 512 ;   dt = T/N  ;
%     disp(['\Delta x = ', num2str(dx, 4),'  \Delta t = ', num2str(dt, 4)])
%     lambda = (1+ sqrt(G)) * dt / dx;
%     disp(['\lambda = ', num2str(lambda, 4)])
    
% Define computational grid for M = 64 and lambda = 0.5
    M = 128    ;   dx = 2*pi/M ;    x_j = 0:dx:2*pi-dx ; 
    T = 4*pi ;  N = round(256*(1+sqrt(2))) ;   dt = T/N  ;
    disp(['\Delta x = ', num2str(dx, 4),'  \Delta t = ', num2str(dt, 4)])
    lambda = (1+sqrt(G)) * dt / dx;
    disp(['\lambda = ', num2str(lambda, 4)]) 
%     
    
    
% Define computational grid for M = 64 and lambda = 1
%     M = 128    ;   dx = 2*pi/M ;    x_j = 0:dx:2*pi-dx ; 
%     T = 4*pi ;  N = 256*G ;   dt = T/N  ;
%     disp(['\Delta x = ', num2str(dx, 4),'  \Delta t = ', num2str(dt, 4)])
%     lambda = G * dt / dx;
%     disp(['\lambda = ', num2str(lambda, 4)])    


% Initial data 
f0 = @(x) exp(-4*((x)-pi).^2);
    IVu = f0(x_j);
    IVh = f0(x_j);
   
 
    ymax = max(IVu);
    ymin = min(IVu);  
    
u0 = IVu ;
h0 = IVh ;

% Compute the first step using an euler method. 
IV1_eul_u  = u0 - (dt/dx)*([u0(2:end) u0(1)] - [u0(end) u0(1:end-1)]) - G*(dt/dx)*([h0(2:end) h0(1)] - [h0(end) h0(1:end-1)]) ;

IV1_eul_h  = h0 - (dt/dx)*([h0(2:end) h0(1)] - [h0(end) h0(1:end-1)]) - (dt/dx)*([u0(2:end) u0(1)] - [u0(end) u0(1:end-1)]);

u1 = IV1_eul_u ;
h1 = IV1_eul_h ; 


for k = 2:4*N/5
    
    % Update the time and the solutions
    t = k*dt ;  
    u = u0 - (dt/dx)*([u1(2:end) u1(1)] - [u1(end) u1(1:end-1)]) - G*(dt/dx)*([h1(2:end) h1(1)] - [h1(end) h1(1:end-1)]) ; 
    h = h0 - (dt/dx)*([h1(2:end) h1(1)] - [h1(end) h1(1:end-1)]) - (dt/dx)*([u1(2:end) u1(1)] - [u1(end) u1(1:end-1)]) ;
    

    % Compute the exact solutions
    Indexu = round(mod(k*lambda*(G),numel(u))) ;
    ExSolu = [IVu(end-Indexu+1:end) IVu(1:end-Indexu) ] ; 
    
    Indexh = round(mod(k*lambda,numel(u))) ;
    ExSolh = [IVu(end-Indexh+1:end) IVu(1:end-Indexh) ] ; 
    
    % Plot the results
     hold off
            plot(x_j, IVh, 'r', 'LineWidth', 2)
            hold on
            plot(x_j, h, 'k', 'LineWidth', 2);
            plot(x_j, u , 'b*','LineWidth', 2);
            plot(x_j,ExSolu,'g-','LineWidth',2);
            plot(x_j,ExSolh,'c--','LineWidth',2);
            axis([0 2*pi -.2 1.2])
            legend('Initial conditions','h','u','Exact u','Exact h')
            xlabel('x');  ylabel('u')
            title(['Shallow water, CNLF: t = ', num2str(t, '%5.3f') ' and  lambda = ' num2str(lambda) ', G =' num2str(G)])
            shg
             pause(1.d-8)
       
             
    u0 = u1 ;
    u1 = u;
    
    h0 = h1;
    h1 = h ;
    
    
end




    
    
    
    
    
    