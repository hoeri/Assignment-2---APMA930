%% Stokes Eddies
%
% Solve the driven cavity problem for Stokes flow in a wedge
% using the streamfunction/vorticity formulation.
%
% Returns the r derivative of psi at r=1. 

function   psi_r = PsyStokesEddies(M,N) 

%
    close all;   clc;
%
% Problem parameters:
    U = -1;
    Rmax = 1;
    alpha = pi/2;
%
% Set up finite difference grid
    M ; dr = Rmax/(M-1);
    N ; dth = alpha/(N-1);
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
    A = PsiOmSys;
    [LL, UU, PP, QQ, RR] = lu(A);
 %   spparms('spumoni', 1)  % write out information about sparse algorithms

    rhs = ConstructRhs(numUn, nP, nO, M, N, Rmax, dr, U);
%%
%
% Solve
    tic
    psivort = A \ rhs ;  % sparse solve
    t = toc;
 %   disp(['Time taken for linear system solve = ', num2str(t)]);
    
%
% Plot
    psi = reshape(psivort(1:numUn/2), size(rg));
    omega = reshape(psivort(numUn/2+1:numUn), size(rg));
    
    
    psi1 = psi(:,end) ;
    psi2 = psi(:,end-1) ;
    
    psi_r = (psi2-psi1)/dr ;
 
    
    