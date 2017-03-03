%% Stokes Eddies
%
% Solve the driven cavity problem for Stokes flow in a wedge
% using the streamfunction/vorticity formulation.
%
% Builds a big system matrix
%
    close all;  clear all; clc;
%
% Problem parameters:
    U = -1;
    Rmax = 1;
    alpha = pi/2;
    Re = 150 ;
    
% Time parameters
Tfinal = 0.5*pi ;
T = 1000 ;
dt = Tfinal/T ;
    
    
%
% Set up finite difference grid
    M =100; dr = Rmax/(M-1);
    N =100; dth = alpha/(N-1);
    [rg, thg] = meshgrid(0: dr :Rmax, ...
                         alpha: -dth: 0);                     
%
% Unknowns and numbering
    numUn = M*N;
    nP = reshape(1:numUn, size(rg));
    nO = reshape(numUn+1:2*numUn, size(rg));
    numUn = 2*numUn;
%
% Build system matrices and rhs
    PsiOmSys = SystemMatComp(numUn, nP, nO, M, N, alpha, dr, dth);
    PsiOmSys = sparse(PsiOmSys);
    % Modify the code so that the Omega term includes the time term (do
    % this in a sparse way)
%     TimeOm = sparse(numUn,numUn) ;
%     for n = (numUn/2+1):numUn
%     TimeOm(n,n) = (3/2) ;
%     end
    PsiOmSys(numUn/2+1:end,numUn/2+1:end) = sparse(eye(numUn/2,numUn/2)) + (2*dt/(3*Re))*PsiOmSys(numUn/2+1:end,numUn/2+1:end) ;
    
    A = PsiOmSys  ;
  
    [LL, UU, PP, QQ, RR] = lu(A);

      % write out information about sparse algorithms
    %spy(PsiOmSys)
    %drawnow
    rhs = ConstructRhs(numUn, nP, nO, M, N, Rmax, dr, U);
%%
%  % Create a 1/r grid
    R1 = zeros(N,M) ;
        for n = 2:M ;
            R1(:,n) = 1/((n-1)*dr) ;
        end
        R1 = reshape(R1,N*M,1) ;




% Start with zero initial conditions
    psivort0 = zeros(size(rhs)) ;
    psivort1 = zeros(size(rhs)) ;


Time = 0;
tic

    Jac0 = [zeros(numUn/2,1) ; R1.*jaco(psivort0(1:numUn/2,1),psivort0(numUn/2+1:end,1),M,N,dr,dth,U)] ;
    Jac1 = [zeros(numUn/2,1) ; R1.*jaco(psivort1(1:numUn/2,1),psivort1(numUn/2+1:end,1),M,N,dr,dth,U)] ;

for k = 1:T

    Time = k*dt ;
    % Construct the jacobians
    
   

    % Solve  
    %psivort = A \ ( 2*psivort1 - 0.5*psivort0 -2*dt*Jac1 + 1*dt*Jac0 +   rhs );  % sparse solve
    
    RHS = (4*psivort1 - 1*psivort0 -4*dt*Jac1 + 2*dt*Jac0)/3 + rhs ;
    RHS(19902:19999) = U/Rmax + 3*U/dr; ;
    
    psivort = (QQ*(UU\(LL\(PP*(RR\(RHS)))))) ;
    
    % Update the previous steps
    psivort0 = psivort1 ;
    psivort1 = psivort ;
    
    Jac0 = Jac1;
    Jac1 = [zeros(numUn/2,1) ; R1.*jaco(psivort(1:numUn/2,1),psivort(numUn/2+1:end,1),M,N,dr,dth,U)] ;
    

end
t = toc;
disp(['Time taken for linear system solve = ', num2str(t)]);

% Plot
    psi = reshape(psivort(1:numUn/2), size(rg));
    omega = reshape(psivort(numUn/2+1:numUn), size(rg));
    figure()
    subplot(1, 2, 1)
        pcolor(rg.*cos(thg), rg.*sin(thg), psi); colorbar;
        shading flat;  colormap(jet);  
        xlabel('x')
        ylabel('y')
        title('Streamfunction')
        axis([0 Rmax 0 Rmax])
        axis square
    subplot(1, 2, 2)
        pcolor(rg.*cos(thg), rg.*sin(thg), omega); colorbar;
        shading flat;  colormap(jet);  
        xlabel('x')
        ylabel('y')
        title('Vorticity')
        axis([0 Rmax 0 Rmax])
        axis square
%
% Look for Eddies!
    figure()
    subplot(1, 2, 1)
        contour(rg.*cos(thg), rg.*sin(thg), psi, [0 0],'k','LineWidth',2); 
        hold on;
        contour(rg.*cos(thg), rg.*sin(thg), psi, 40); 
%         c = linspace(0, 4.4d-4, 40);
%         contour(rg.*cos(thg), rg.*sin(thg), psi, c); 
        shading flat;  colormap(jet);  
        xlabel('x')
        ylabel('y')
        title('Streamfunction')
        axis([0 Rmax 0 Rmax])
        axis square
    subplot(1, 2, 2)
        c = linspace(-20, 20, 40);
        contour(rg.*cos(thg), rg.*sin(thg), omega, c);
        shading flat;  colormap(jet);  
        xlabel('r')
        ylabel('\theta')
        title(['Vorticity, at time = ' num2str(Time)])
        axis([0 Rmax 0 Rmax])
        axis square
