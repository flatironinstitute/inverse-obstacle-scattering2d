function [der] = curvature_directional_der(src_info,h)
%CURVATURE_DIRECTIONAL_DER compute the directional derivative of the 
% curvature in the specified direction 
%
% Let H( [x(t),y(t)] ) be the signed curvature for a curve, i.e. 
%
% H( [x(t),y(t)] ) = (x'(t)y''(t)-x''(t)y'(t))/( x'(t)^2+y'(t)^2 )^(3/2)
%
% This code evaluates the directional derivative
%
% lim delta -> 0 of 
%
% (H( [x(t),y(t)] + delta*h(t)*[nx(t),ny(t)] ) - H([x(t),y(t)])) / delta
%
% where h(t) is a supplied function on the curve and [nx(t),ny(t)] is the
% normal to the curve, i.e. 
%
% nx(t) = y'(t)/(x'(t)^2+y'(t)^2)^(1/2)
% ny(t) = -x'(t)/(x'(t)^2+y'(t)^2)^(1/2)
%
% This is done using a direct formula (obtained via Mathematica).
% 
% -h*H^2 + h'*(x'*x''+y'*y'')/(x'^2+y'^2)^2 - h''/(x'^2+y'^2)
%
% input:
%
% src_info - the curve parameterization in the standard format
% h - the variable direction 
%
% output:
% 
% the directional derivative of the curvature in the direction h 
%

xs = src_info.xs; xs = xs(:);
ys = src_info.ys; ys = ys(:);

Dmat = src_info.Der_param;

dxs = Dmat*xs;
dys = Dmat*ys;
d2xs = Dmat*dxs;
d2ys = Dmat*dys;

dh = Dmat*h;
d2h = Dmat*dh;

der = -h.*(src_info.H(:)).^2 + ...
    dh.*(dxs.*d2xs + dys.*d2ys)./(dxs.^2+dys.^2).^2 - ...
    d2h./(dxs.^2+dys.^2);

