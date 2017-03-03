%% Construct the Jacobian Matrix for Psi and Omega (note that the 1/r term is not included)


function       Jac  =  jaco(newPsi,newOm,M,N,dr,dth,U) 



Jac = zeros(N,M) ;

% Create grids out of the Matrices ;
Psi = reshape(newPsi,N,M) ;
Om = reshape(newOm,N,M) ;

% Compute the r derivatives seperatively
psi_th = zeros(N,M) ;
om_th = zeros(N,M) ;

psi_th(2:end-1,:) = ( Psi(3:end,:) - Psi(1:end-2,:) )./(2*dth) ;
om_th(2:end-1,:) = ( Om(3:end,:) - Om(1:end-2,:) )./(2*dth) ;


% Compute the theta derivatives
psi_r = zeros(N,M) ;
om_r = zeros(N,M) ;

psi_r(:,2:end-1) = (Psi(:,3:end) - Psi(:,1:end-2))./(2*dr) ;
om_r(:,2:end-1) =  (Om(:,3:end) - Om(:,1:end-2))./(2*dr) ;


% Fill in the interior of the Jacobian. 
Jac = -(om_r.*psi_th - om_th.*psi_r) ; 



% disp(size(psi_r1))
% disp(size(om_r1))
% disp(size(psi_th1))
% disp(size(om_th1))



% Compute the boundary conditions at the origin, we use a one sided
% difference for this

%One sided r derivatives 
% psi_r1 = (Psi(:,2) - Psi(:,1))/dr ;
% om_r1 = (Om(:,2) - Om(:,1))/dr ; 
% 
% psi_th1 = [(Psi(1,1)-Psi(2,1))/dth ; -(Psi(3:N,1) - Psi(1:N-2,1))/(2*dth); (Psi(N-1,1)-Psi(N,1))/dth] ;
% om_th1 = [ (Om(1,1)-Om(2,1))/dth; (Om(3:N,1) - Om(1:N-2,1))/(2*dth); (Om(N-1,1)-Om(N,1))/dth] ;
% 
% 
% Jac(:,1) = om_r1.*psi_th1 - om_th1.*psi_r1 ;

% % Compute the boundary conditions at theta = 0
% 
% psi_th1 = (Psi(1,:) - Psi(2,:))/dth ;
% om_th1  = (Om(1,:) - Om(2,:))/dth ; 
% 
% psi_r1 = [ (Psi(1,2) -Psi(1,1))/dr (Psi(1,3:M) - Psi(1,1:M-2))/(2*dr) (Psi(1,M) -Psi(1,M-1))/dr] ;
% om_r1 = [ (Om(1,2) -Om(1,1))/dr   (Om(1,3:M) - Om(1,1:M-2))/(2*dr)   (Om(1,M) -Om(1,M-1))/dr] ;

% disp(size(psi_r1))
% disp(size(om_r1))
% disp(size(psi_th1))
% disp(size(om_th1))


%Jac(1,:) = om_r1.*psi_th1 - om_th1.*psi_r1 ;


% 
% % Compute the boundary conditions at theta = 0
% 
% psi_th = (Psi(2,:) - Psi(1,:))/dth ;
% om_th  = (Om(2,:) - Om(1,:))/dth ;
% 
% psi_r = [(Psi(1,2:M)-Psi(1,1:M-1))/dr   Psi(1,M)];
% om_r = [(Om(1,2:M)-Om(1,1:M-1))/dr  Om(1,M) ] ;
% 
%         % Fill in the top of the Jacobian
% Jac(1,1:M) = om_r.*psi_th - om_th.*psi_r ;
% 
% 
% % Compute the boundary conditions at theta = alpha
% 
% psi_th = (Psi(N,:) - Psi(N-1,:))/dth ;
% om_th  = (Om(N,:) - Om(N-1,:))/dth ;
% 
% psi_r = [(Psi(N,2:M)-Psi(N,1:M-1))/dr  Psi(N,M)];
% om_r = [(Om(N,2:M)-Om(N,1:M-1))/dr  Om(N,M) ] ;
% 
% 
%         % Fill in the bottom of the jacobian
% Jac(N,1:M) = om_r.*psi_th - om_th.*psi_r ;
% 
% 
% 
% Compute the boundary conditions at r = Rmax, keeping in mind that the
% theta derivative of psi is zero there, and the r derivative of psi is -U

% om_th = [-(Om(1:N-2,M) - Om(3:N,M))/(2*dth) ] ;
% 
% disp(size(om_th))
% 
% %om_th = [om_th ; (Om(M,N) - CornerOm)/dth ];
% 
% 
% 
% Jac(2:N-1,M) = -U*om_th ;

%psi_th = [(Psi(2:N,M) - Psi(1:N-1,M))/dth ; (Psi(N,M)-Psi(N-1,M))/dth ]  ;

% om_th1 = [ (Om(1,M)-Om(2,M))/dth; (Om(3:N,M) - Om(1:N-2,M))/(2*dth); (Om(N-1,M)-Om(N,M))/dth] ;
% 
%  psi_r1 = (Psi(:,M) - Psi(:,M-1))/dr ; 
% % om_r =  (Om(:,M) - Om(:,M-1))/dr ;
% 
% Jac(:,M) =  - om_th1.*psi_r1 ;




% Compute the (rescaled) Jacobian

Jac = reshape(Jac,N*M,1) ;



