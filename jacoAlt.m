%% Alternate way of gettin the jacobian. 


function   Jac = jacoAlt(newPsi,newOm,M,N,dr,dth,U)

Jac = zeros(M*N,1) ;


for n = 2:N-1
    
    
    
    for m = 2:M-1
        
        
        psi_r = (newPsi(N*n + m + 1) - newPsi(N*n + m - 1))/(2*dr) ;
        om_r = (newOm(N*n + m + 1) - newOm(N*n + m - 1))/(2*dr) ;
        
        psi_th = (newPsi(M*m + n + 1 ) - newPsi(M*m + n - 1 ))/(2*dth) ;
        om_th = (newOm(M*m + n + 1 ) - newOm(M*m + n - 1 ))/(2*dth) ;
        
        Jac(N*(n-1) + m) = om_r.*psi_th - om_th.*psi_r ;
        
    end
    
end