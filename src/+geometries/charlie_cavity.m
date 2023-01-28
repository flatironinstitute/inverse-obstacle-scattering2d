function [src_info] = charlie_cavity(j,n)
%  j must be a real number between 0 and 10

    if(j<0 || j>11)
        fprintf('j must be between 0 and 10');
        return
    end
    
    src_info = [];
    ts = (0:(n-1))/n*2*pi;
    rx=(1-(j-1)/11)+0.03*(j-1);
    ry=(1-(j-1)/11)+0.8*(j-1);
    
    
    ct = cos(ts);
    st = sin(ts);
    
    xt = rx*ct + 20;
    dxdt = -rx*st;
    d2xdt2 = -rx*ct;
    yt = ry*st;
    dydt = ry*ct;
    d2ydt2 = -ry*st;
    
    rmax = 400*rx*rx/(ry*ry-rx*rx) + 400 + ry^2;
    
    a = (1+ 0.7*(j-1))/2;
    C = 1/rmax^a;
    rho_t = C*(xt.^2 + yt.^2).^(a);
    drho_t_dt = 2*C*a*(xt.^2 + yt.^2).^(a-1).*(xt.*dxdt + yt.*dydt);
    d2rho_t_dt2 = C*((a-1)*a*(xt.^2+yt.^2).^(a-2).*(2*xt.*dxdt + 2*yt.*dydt).^2 + ...  
         a*(xt.*xt + yt.*yt).^(a-1).*(2*xt.*d2xdt2 + 2*dxdt.*dxdt + 2*yt.*d2ydt2 + 2*dydt.*dydt));
    
    
    
    rfac = 2*a;
    theta_t = rfac*atan2(yt,xt);
    dtheta_t_dt = rfac*(xt.*dydt - yt.*dxdt)./(xt.*xt + yt.*yt);
    d2theta_t_dt2 = rfac*(xt.*yt.*(2*dxdt.*dxdt + yt.*d2ydt2 -2*dydt.*dydt) - ...
            xt.*xt.*(yt.*d2xdt2 + 2*dxdt.*dydt) + ...
             yt.*yt.*(2*dxdt.*dydt - yt.*d2xdt2) + xt.*xt.*xt.*d2ydt2)./ ...
         (xt.*xt + yt.*yt).^2;
    
    
    
    src_info.xs = rho_t.*cos(theta_t);
    src_info.ys = rho_t.*sin(theta_t);
    src_info.dxs = drho_t_dt.*cos(theta_t) - rho_t.*sin(theta_t).*dtheta_t_dt;
    src_info.dys = drho_t_dt.*sin(theta_t) + rho_t.*cos(theta_t).*dtheta_t_dt;
    
    d2xs = d2rho_t_dt2.*cos(theta_t) - 2*drho_t_dt.*sin(theta_t).*dtheta_t_dt -... 
              rho_t.*(cos(theta_t).*dtheta_t_dt.^2 + sin(theta_t).*d2theta_t_dt2);
    d2ys = d2rho_t_dt2.*sin(theta_t) + 2*drho_t_dt.*cos(theta_t).*dtheta_t_dt +... 
              rho_t.*(-sin(theta_t).*dtheta_t_dt.^2 + cos(theta_t).*d2theta_t_dt2);      
    
    
    src_info.ds = sqrt(src_info.dxs.^2 + src_info.dys.^2);
    src_info.H = (src_info.dxs.*d2ys - src_info.dys.*d2xs)./(src_info.ds).^3; 
    h = 2*pi/n;
    src_info.h = h;
    src_info.L = sum(src_info.ds)*h;
    src_info.paramL = 2*pi;
    rsc = 2*pi/src_info.paramL;
    src_info.Der_param = specdiffmat_ds(n,ones(n,1))*rsc;
    src_info.Der = specdiffmat_ds(n,src_info.ds)*rsc;
end