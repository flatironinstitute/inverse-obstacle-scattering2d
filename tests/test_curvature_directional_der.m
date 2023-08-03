%test_curvature_directional_der
%


clear all; clc; close all;

n  = 300;
path_to_ios2d = '../../inverse-obstacle-scattering2d/';
addpath(path_to_ios2d);
addpath(genpath_ex(path_to_ios2d));

% max number of wiggles
nc = 3;

% parameters 'a'
coefs = zeros(2*nc+1,1);
coefs(1) = 1;
coefs(nc+1) = 0.3;

src_info = geometries.starn(coefs,nc,n);
xs = src_info.xs; xs = xs(:);
ys = src_info.ys; ys = ys(:);
diffmat = src_info.Der_param;

ts = linspace(0,2*pi,n+1); ts = ts(1:end-1).'; ts = ts(:);

h = cos(3*ts) + 0.5*sin(4*ts);

ntimes = 6;

der = curvature_directional_der(src_info,h);

for ii = 1:ntimes
    delta = 10^(-ii);

    ds = sqrt(src_info.dxs.^2 + src_info.dys.^2);
    nx = src_info.dys./ds; ny = -src_info.dxs./ds;
    nx = nx(:); ny = ny(:);

    xplus = xs + delta*h.*nx;
    yplus = ys + delta*h.*ny;
    xminus = xs - delta*h.*nx;
    yminus = ys - delta*h.*ny;

    dxplus = diffmat*xplus;
    dyplus = diffmat*yplus;
    d2xplus = diffmat*dxplus;
    d2yplus = diffmat*dyplus;

    dxminus = diffmat*xminus;
    dyminus = diffmat*yminus;
    d2xminus = diffmat*dxminus;
    d2yminus = diffmat*dyminus;

    Hplus = (dxplus.*d2yplus - dyplus.*d2xplus)./ ...
        ( sqrt(dxplus.^2+dyplus.^2)).^3;

    Hminus = (dxminus.*d2yminus - dyminus.*d2xminus)./ ...
        ( sqrt(dxminus.^2+dyminus.^2)).^3;

    derapprox = (Hplus - Hminus)/(2*delta);

    norm(der-derapprox)
end
    

