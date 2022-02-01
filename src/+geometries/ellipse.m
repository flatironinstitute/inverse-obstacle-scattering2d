function [src_info] = ellipse(a,b,n)
   
    
    src_info = [];
    ts = (0:(n-1))/n*2*pi;
    src_info.xs = a*cos(ts);
    src_info.ys = b*sin(ts);
    src_info.ds = sqrt((b*cos(ts)).^2 + (a*sin(ts)).^2);
    src_info.dxs = -a*sin(ts);
    src_info.dys = b*cos(ts);
    d2xdt2 = -a*cos(ts);
    d2ydt2 = -b*sin(ts);
    src_info.H = (src_info.dxs.*d2ydt2 - src_info.dys.*d2xdt2)./(src_info.ds).^3; 
    h = 2*pi/n;
    src_info.h = h;
    src_info.L = sum(src_info.ds)*h;
    src_info.paramL = 2*pi;
    rsc = 2*pi/src_info.paramL;
    src_info.Der_param = specdiffmat_ds(n,ones(n,1))*rsc;
    src_info.Der = specdiffmat_ds(n,src_info.ds)*rsc;
end