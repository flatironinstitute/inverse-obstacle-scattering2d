function [src_info] =  cone(alpha,beta,gamma,n)

    src_info = [];
    ts = -pi/2 + (0:(n-1))/n*pi;
    src_info.xs = sqrt(alpha.^2 + sin(ts).^2);
    src_info.ys = gamma*sin(ts).^2 - beta*cos(ts).*sin(ts);
    src_info.dxs = sin(ts).*cos(ts)./sqrt(alpha^2 + sin(ts).^2);
    src_info.dys = -beta*(cos(ts).^2 - sin(ts).^2) + 2*gamma*sin(ts).*cos(ts);
    src_info.ds = sqrt(src_info.dxs.^2 + src_info.dys.^2);
    d2xdt2 = (-cos(ts).^2.*sin(ts).^2 + ...
      (cos(ts).^2-sin(ts).^2).*(alpha^2 + ...
      sin(ts).^2))./(alpha^2 + sin(ts).^2).^1.5;
    d2ydt2 = 4*beta*cos(ts).*sin(ts)+2*gamma*(cos(ts).^2 - sin(ts).^2);
    src_info.H = (src_info.dxs.*d2ydt2 - src_info.dys.*d2xdt2)./(src_info.ds).^3; 
    h = pi/n;
    src_info.h = h;
    src_info.L = sum(src_info.ds)*h;
    src_info.paramL = pi;
    rsc = 2*pi/src_info.paramL;
    src_info.Der_param = specdiffmat_ds(n,ones(n,1))*rsc;
    src_info.Der = specdiffmat_ds(n,src_info.ds)*rsc;
end