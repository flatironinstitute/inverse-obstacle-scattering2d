function [src_info] = larrycup(a,b,n,m)
    if(nargin == 3)
        m = n
    end
    nhalf = ceil(m/2);
    s = ((1:nhalf)-0.5)/nhalf * pi;  % note half-offset, needed for easy reflection abt z
    r = 1 - a*erf((s-pi/2)/a);  % radius: starts at 1+a, ends at 1-a
    c = a; %*(1-b/pi);  % is theta rounding scale
    sabs = @(x) exp(-(x/c).^2)*c/sqrt(pi)+x.*erf(x/c); % c-smoothed absval
    th = b-a + 2*(1-(b-a)/pi)*sabs(s-pi/2);
    rho = r.*sin(th); z = r.*cos(th);  % theta down from z axis as in 3D cyl coords
    z = z*1.2;  % vert stretch! makes ellipse cavity
    Z = [rho -rho(end:-1:1)] + 1i*[z z(end:-1:1)]; % complex coords of full curve
    zhat = fft(Z(:))/m;
    t1 = (0:(n-1))/n;
    h = 1.0/n;
    xy = fourierZ(zhat,t1);
    dxydt = fourierZp(zhat,t1);
    src_info = [];
    src_info.xs = real(xy);
    src_info.ys = imag(xy);
    src_info.ds = abs(dxydt);
    src_info.dxs = real(dxydt);
    src_info.dys = imag(dxydt);
    src_info.h = h;
    src_info.L = sum(src_info.ds)*h;
    src_info.paramL = 1.0;
    rsc = 2*pi/src_info.paramL;
    src_info.Der_param = specdiffmat_ds(n,ones(n,1))*rsc;
    src_info.Der = specdiffmat_ds(n,src_info.ds)*rsc;
    src_info.H = rla.get_curvature(src_info);

end
