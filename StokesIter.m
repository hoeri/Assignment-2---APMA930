%% Stokes Eddies - Iteration version
%
% Solve the driven cavity problem for Stokes flow in a wedge
% using the streamfunction/vorticity formulation.
%
% Solves two separate BVPs for Psi and Omega & iterates.
%
    close all;  clear; clc;
%
% Problem parameters:
    U = -1;
    Rmax = 1;
    alpha = pi/2;
    Re = 50 ;
%
% Set up finite difference grid
    M = 60; dr = Rmax/(M-1);
    N = 50; dth = alpha/(N-1);
    [rg, thg] = meshgrid(0: dr :Rmax, ...
                         alpha: -dth: 0);                     
%
% Unknowns and numbering
    numUn = M*N;
    nP = reshape(1:numUn, size(rg)); % numbering of psi unknowns
    nO = reshape(numUn+1:2*numUn, size(rg)); % numbering of omega unknowns
%
% Build system matrices and rhs
    PsiOmSys = SystemMat(2*numUn, nP, nO, M, N, alpha, dr, dth);
    PsiOmSys = sparse(PsiOmSys);
    % 
    % Block Matrices
        AP = PsiOmSys(1:numUn, 1:numUn);
        BO = PsiOmSys(1:numUn, numUn+1:2*numUn);
        BP = PsiOmSys(numUn+1:2*numUn, 1:numUn);
        AO = PsiOmSys(numUn+1:2*numUn, numUn+1:2*numUn);
        
        % Construct 1/r matrix (it is a grid of the wedge of 1 divided by the r value at
        % that point on the grid)
        R1 = zeros(N,M) ;
        for n = 2:M ;
            R1(:,n) = 1/(n*dr) ;
        end
        R1 = reshape(R1,numUn,1) ;
        
        
    %
    % LU Decomposition
    tic;
    [LP, UP, PP, QP, RP] = lu(AP);
    [LO, UO, PO, QO, RO] = lu(AO);
    tlu = toc;
    disp(['Time for LU Decomposition = ', num2str(tlu)])
    
        
    spy(PsiOmSys)
    drawnow
    spparms('spumoni', 0)
    rhs = ConstructRhs(2*numUn, nP, nO, M, N, Rmax, dr, U);
        rhsP = rhs(1:numUn); rhsO = rhs(numUn+1:2*numUn);
%
% Initialize psi and omega
    newPsi = zeros(numUn, 1);
    tilPsi = zeros(numUn, 1);
    newOm = zeros(numUn, 1);
    tilOm = zeros(numUn, 1);
    tempPsi = zeros(numUn, 1);
    tempOm = zeros(numUn, 1) ;
    
     
%
% Relaxation parameters
    alpha = 0.99; beta = 0.01;
%
% Solve via iteration
    
    tic;
    itmax = 1000;
    normPsi = zeros(1, itmax);
    normOm = zeros(1, itmax);
    tol = 1.d-9;
    normPsi(1) = 1;
    normOm(1) = 1;
    iter = 2;
    
    while (normPsi(iter-1) > tol || normOm(iter-1) > 1.d-2*tol) && ...
            iter < itmax
% 
% Calculate \psi^\tilde_{m+1}, \omega^\tilde_{m+1}

%  System matrix not factorized
          Jac = jaco(newPsi,newOm,M,N,dr,dth,U) ;
          Jac = R1.*Jac ;
        tilPsi = AP\(rhsP - BO*newOm);
        tilOm = AO\(rhsO - BP*newPsi + Re*Jac);
        
%  System maqtrix factorized
%         cP = PP * (RP \ (rhsP - BO*newOm) ) ;
%         tilPsi = QP * (UP \ (LP \ cP));
%         cO = PO * (RO \ (rhsO - BP*newPsi) ) ;
%         tilOm = QO * (UO \ (LO \ cO) );
        
    %
    % Relax - Calculate \psi_{m+1}, \omega_{m+1}
        tempPsi = alpha*tilPsi + (1-alpha)*newPsi;
        tempOm = beta*tilOm + (1-beta)*newOm;
    
    % Check - relative nomrs of  \psi_{m+1}-\psi_m, \omega_{m+1}-\omega_m
        normPsi(iter) = norm(tempPsi-newPsi)/norm(tempPsi);
        normOm(iter) = norm(tempOm - newOm)/norm(tempOm);
        newPsi = tempPsi;
        newOm = tempOm;
       disp(['Iteration: ', num2str(iter), ...
              ' Residual Psi = ', num2str(normPsi(iter)), ...
              ' Residual Omega = ', num2str(normOm(iter))])
        iter = iter + 1;
    end
    iter = iter - 1;
    t = toc;
    disp(' ')
    disp(['Time taken for solve = ', num2str(t)]);
    spy(reshape(Jac,N,M))
%
% Plot convergence
    figure()
    loglog(1:iter, normPsi(1:iter), 'r')
    hold on
    loglog(1:iter, normOm(1:iter), 'b')
    xlabel('iter')
    ylabel('|| x^{m+1} - x^{m} ||/ || x^{m+1} ||')
    legend('\psi', '\omega', 'location', 'NorthEast')
%
% Plot
    psi = reshape(newPsi, size(rg));
    omega = reshape(newOm, size(rg));
    figure()
    subplot(1, 2, 1)
        pcolor(rg.*cos(thg), rg.*sin(thg), psi); colorbar;
        shading flat;  colormap(jet);  
        xlabel('x')
        ylabel('y')
        title('Streamfunction')
        axis([0 1 0 1])
        axis square
    subplot(1, 2, 2)
        pcolor(rg.*cos(thg), rg.*sin(thg), omega); colorbar;
        shading flat;  colormap(jet);  
        xlabel('x')
        ylabel('y')
        title('Vorticity')
        axis([0 1 0 1])
        axis square
%%
% Look for Eddies!
    figure()
    subplot(1, 2, 1)
        contour(rg.*cos(thg), rg.*sin(thg), psi, [0 0],'k','LineWidth',2); 
        hold on;
        contour(rg.*cos(thg), rg.*sin(thg), psi, 40);         
%         c = linspace(-5d-7, 0, 40);
%         contour(rg.*cos(thg), rg.*sin(thg), psi, c); 
        hold on
        shading flat;  colormap(jet);  
        xlabel('x')
        ylabel('y')
        title('Streamfunction')
        axis([0 1 0 1])
        axis square
    subplot(1, 2, 2)
        c = linspace(-20, 20, 40);
        contour(rg.*cos(thg), rg.*sin(thg), omega, c);
        shading flat;  colormap(jet);  
        xlabel('x')
        ylabel('y')
        title('Vorticity')
        axis([0 1 0 1])
        axis square
    
