function StabilityAdvecDiff (imethod)
    close all; 
    
    [alpha_dt, beta_dt] = meshgrid([-3:.025:0], [-3:.025:0]);
    
    switch imethod
        case 1 % forward Euler
            A = 1 + alpha_dt + 1i*beta_dt;
        case 2 % imex Euler
            A = (1 + 1i*beta_dt)./(1 - alpha_dt);
        case 3 % spdf
            a = 3-2*alpha_dt;
            b = -4-4*1i*beta_dt;
            c = 1+2*1i*beta_dt;
            xi_1 = -b./(2*a) + (1./(2*a)).*sqrt(b.^2 - 4*a.*c);
            xi_2 = -b./(2*a) - (1./(2*a)).*sqrt(b.^2 - 4*a.*c);
            A = max(abs(xi_1),abs(xi_2));
        case 4 % cnlf
            a = 1-alpha_dt;
            b = -2*1i*beta_dt;
            c = -1-alpha_dt;
            xi_1 = -b./(2*a) + (1./(2*a)).*sqrt(b.^2 - 4*a.*c);
            xi_2 = -b./(2*a) - (1./(2*a)).*sqrt(b.^2 - 4*a.*c);
            A = max(abs(xi_1),abs(xi_2));
            
        case 5 % CNAB
            mu = 0.1;
            C = 10;
            tau1 = (1/2*1i)*(-1i*alpha_dt*mu+3*beta_dt*C-2*1i+sqrt(-mu^2*alpha_dt.^2-(10*1i)*C*mu*alpha_dt*beta_dt+9*C^2*beta_dt.^2-4*alpha_dt*mu-(4*1i)*beta_dt*C-4))./(mu*alpha_dt-2) ;
            tau2 = (1/2)*(-(3*1i)*beta_dt*mu-alpha_dt*C-2+sqrt(-9*C^2*beta_dt.^2+(10*1i)*C*mu*alpha_dt.*beta_dt+alpha_dt.^2*mu^2+(4*1i)*beta_dt*C+4*alpha_dt*mu+4))./(mu*alpha_dt-2) ;
            A = abs(tau2) ; 
            
    end
    pcolor(alpha_dt, beta_dt, abs(A))
    shading flat 
    colormap(jet)
    colorbar
    hold on
    contour (alpha_dt, beta_dt, abs(A), [1 1],'k', 'LineWidth', 2) % plot c=1 contour line
    xlabel('\alpha \Delta t')
    ylabel('\beta \Delta t')
    caxis([0 2])
    switch imethod
        case 1
%            axis([-3 0 -2 0])
            title('Stability - Forward Euler')
        case 2
            title('Stability - Backwards Euler')
        case 3
            title('Stability - SBDF')
        case 4
            title('Stability - CNLF')
        case 5
            title('Stability - CNAB')
    end
end