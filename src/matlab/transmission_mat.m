function [xmat] = transmission_mat(zks,a,b,norder,h,srcinfo)
%
%
%  Representation
%    u_{i}  = -(1/b_{i}) S_{k_{i}}[\sigma] + \frac{1}{b_{i}} D_{k_{i}}[\mu]
%    
%  PDE 
%
%  [au]/q = f/q, [b du/dn] = g
%  q = 0.5*(a_{1}/b_{1} + a_{2}/b_{2})
%
%  Unknown ordering \sigma_{1},\mu_{1},\sigma_{2},\mu_{2}....
%  Output ordering [au]/q_{1}, [b du/dn]_{1}, [au]/q_{2}, [b du/dn]_{2}...
%
%  Input: 
%    zks(2) - complex
%       Helmholtz parameters
%    a(2) - complex
%       scaling for jump in u
%    b(2) - complex
%      scaling for jump in dudn
%    norder - Alpert quadrature rule order
%    srcinfo(6,n) - source info
%       srcinfo(1:2,:) - locations
%       srcinfo(3:4,:) - normals
%       srcinfo(5,:) - dsdt
%       srcinfo(6,:) - curvature
%
%    h - step size in parameter space
%
%  Output:
%    xmat - complex(2*n,2*n)
%       transmission matrix

  [m,n] = size(srcinfo);
  assert(m==6,'srcinfo must be of shape (5,n)');
  xmat = complex(zeros(2*n),0);
  nsys = 2*n;
  mex_id_ = 'trans_mat(i int[x], i int[x], i double[x], i double[xx], i dcomplex[x], i dcomplex[x], i dcomplex[x], io dcomplex[xx])';
[xmat] = kern_mats(mex_id_, n, norder, h, srcinfo, zks, a, b, xmat, 1, 1, 1, 6, n, 2, 2, 2, nsys, nsys);
end
  
