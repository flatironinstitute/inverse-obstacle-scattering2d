function [src_info] = starn(coefs,nc,n)  
    src_info = [];
    ts = (0:(n-1))/n*2*pi;
    
    coefs_use = coefs(:);
    
    if (nc == 0)
        r = coefs(1)*ones(size(ts));
        drdt = zeros(size(r));
        d2rdt2 = zeros(size(r));
    else
        r = (cos(ts'*(0:nc))*coefs_use(1:(nc+1)) + sin(ts'*(1:nc))*coefs_use((nc+2):end)).';
        drdt = (-sin(ts'*(0:nc))*(coefs_use(1:(nc+1)).*(0:nc).') + cos(ts'*(1:nc))*(coefs_use((nc+2):end).*(1:nc).')).';
        d2rdt2 = (-cos(ts'*(0:nc))*(coefs_use(1:(nc+1)).*(0:nc).^2.') - sin(ts'*(1:nc))*(coefs_use((nc+2):end).*(1:nc).^2.')).';
    end
    src_info.xs = r.*cos(ts);
    src_info.ys = r.*sin(ts);
    src_info.dxs = drdt.*cos(ts) - r.*sin(ts);
    src_info.dys = drdt.*sin(ts) + r.*cos(ts);
    src_info.ds = sqrt(src_info.dxs.^2 + src_info.dys.^2);
    d2xdt2 = -2*drdt.*sin(ts) + d2rdt2.*cos(ts) - r.*cos(ts);
    d2ydt2 = 2*drdt.*cos(ts) + d2rdt2.*sin(ts) - r.*sin(ts);
    src_info.H = (src_info.dxs.*d2ydt2 - src_info.dys.*d2xdt2)./(src_info.ds).^3; 
    h = 2*pi/n;
    src_info.h = h;
    src_info.L = sum(src_info.ds)*h;
    src_info.paramL = 2*pi;
    rsc = 2*pi/src_info.paramL;
    src_info.Der_param = specdiffmat_ds(n,ones(n,1))*rsc;
    src_info.Der = specdiffmat_ds(n,src_info.ds)*rsc;
end