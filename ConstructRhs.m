function [r] = ConstructRhs(numUn, nP, nO, M, N, R, dr, U)
    r = zeros(numUn, 1);
%
%  Boundary conditions at right: omega given
    for jrow = 2: N-1
        ijO = nO(jrow, M);
        r(ijO) = U/R + 3*U/dr;
    end  
end