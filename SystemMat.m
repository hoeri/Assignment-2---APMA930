function [PsiOm] = SystemMat(numUn, nP, nO, M, N, alpha, dr, dth)
    PsiOm = zeros(numUn, numUn);
%
%  Boundary conditions at origin 
    for jrow = 1: N
        ij = nP(jrow, 1);
        ijp1 = nP(jrow, 2);
        PsiOm(ij, ij) = (N-2);
        PsiOm(ij, ijp1) = -1;
        ijO = nO(jrow, 1);
        PsiOm(ij, ijO) = -alpha*(dr/2)^2/dth;
        ijOp1 = nO(jrow, 2);
        PsiOm(ijO, ijO) = 1/dr;
        PsiOm(ijO, ijOp1) = -1/dr;
    end
    for icol = 2: M-1
        r_i = (icol-1)*dr;
%
%  Boundary condition at top - Dirichlet 
        ij = nP(1, icol);
        ijm1 = nP(2, icol);
        ijm2 = nP(3, icol);
        PsiOm(ij, ij) = 1;
        ijO = nO(1, icol);
        PsiOm(ijO, ijO) = 1;
        PsiOm(ijO, ijm1) = 4/(r_i*dth)^2;
        PsiOm(ijO, ijm2) = -0.5/(r_i*dth)^2;
%
%  Interior points - 5 point stencil
        for jrow = 2: N-1
            ij = nP(jrow, icol);
            ip1j = nP(jrow, icol+1);
            im1j = nP(jrow, icol-1);
            ijp1 = nP(jrow-1, icol);
            ijm1 = nP(jrow+1, icol);
            PsiOm(ij, im1j) = -(r_i - 0.5*dr)/(r_i*dr^2);
            PsiOm(ij, ip1j) = -(r_i + 0.5*dr)/(r_i*dr^2);
            PsiOm(ij, ijm1) = -1/(r_i*dth)^2;
            PsiOm(ij, ijp1) = -1/(r_i*dth)^2;
            PsiOm(ij, ij) = 2/dr^2 + 2/(r_i*dth)^2;
            ijO = nO(jrow, icol);
            PsiOm(ij, ijO) = -1;   % lapl(psi) = -Omega;
            ij = ijO;
            ip1j = nO(jrow, icol+1);
            im1j = nO(jrow, icol-1);
            ijp1 = nO(jrow-1, icol);
            ijm1 = nO(jrow+1, icol);
            PsiOm(ij, im1j) = -(r_i - 0.5*dr)/(r_i*dr^2);
            PsiOm(ij, ip1j) = -(r_i + 0.5*dr)/(r_i*dr^2);
            PsiOm(ij, ijm1) = -1/(r_i*dth)^2;
            PsiOm(ij, ijp1) = -1/(r_i*dth)^2;
            PsiOm(ij, ij) = 2/dr^2 + 2/(r_i*dth)^2;
        end
%
%  Boundary condition at bottom - Dirichlet
        ij = nP(N, icol);
        ijp1 = nP(N-1, icol);
        ijp2 = nP(N-2, icol);
        PsiOm(ij, ij) = 1;
        ijO = nO(N, icol);
        PsiOm(ijO, ijO) = 1;
        PsiOm(ijO, ijp1) = 4/(r_i*dth)^2;
        PsiOm(ijO, ijp2) = -0.5/(r_i*dth)^2;
    end
%
%  Boundary conditions at right - Dirichlet
    for jrow = 1: N
        ij = nP(jrow, M);
        im1j = nP(jrow, M-1);
        im2j = nP(jrow, M-2);
        PsiOm(ij, ij) = 1;
        ijO = nO(jrow, M);
        PsiOm(ijO, ijO) = 1;
        PsiOm(ijO, im1j) = 4/dr^2;
        PsiOm(ijO, im2j) = -0.5/dr^2;
    end  
end