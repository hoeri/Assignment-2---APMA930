%% Revovered boundary conditions for Stokes Flow in a wedge

% We see how well the boundary conditions that the theta derivative equals
% one, is recovered from the artificial boundary conditions that was
% imposed at r=Rmax. 


% We keep the grid points in r and theta equal, and increase them in
% increments of 25, from 50 to 300

% Set number of grid points
clf

M = 100 ;
N = 100 ;

% Set angle
alpha = pi/2
Rmax = 1 ;
    
   
AvgPsi = zeros(5,5) ;


for R = 1:5
    for th = 1:5    
        AvgPsi(R,th) = mean(PsyStokesEddies(100*R,100*th)) ;       
    end
end

AvgPsi1 = AvgPsi ;
AvgPsi2 = AvgPsi1 ;

% for n = 1:5
%     
%     AvgPsi2(:,n) = AvgPsi1(:,6-n) ;
%     AvgPsi(n,:) = AvgPsi2(6-n,:) ;
%     
% end

clf
hold on 

surf(1:5,1:5,AvgPsi) 

title('Average Psi derivative with respect to r, for varying dr and dth')

xlabel('dth')
ylabel('dr')
zlabel('Average derivative')

hold off




