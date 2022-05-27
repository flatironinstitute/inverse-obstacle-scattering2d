function [xmat] = ddiff0_neu_mat(zpars,norder,h,srcinfo)
%
%  Representation:
%    u = D_{k1} - D_{0} [\sigma] 
%
%  Data returned:
%    Neumann data (du/dn)
%
%
%  Input: 
%    zpars(1) - 
%      zpars(1) - Helmholtz paramter (k1)
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
%    xmat - complex(n,n)
%      d_{k}' - D_{0}'
%       

  [m,n] = size(srcinfo);
  assert(m==6,'srcinfo must be of shape (5,n)');
  xmat = complex(zeros(n),0);
  mex_id_ = 'ddiff0_neu_mat(i int[x], i int[x], i double[x], i double[xx], i dcomplex[x], io dcomplex[xx])';
[xmat] = kern_mats(mex_id_, n, norder, h, srcinfo, zpars, xmat, 1, 1, 1, 6, n, 1, n, n);
end
%  
%
%
%
%
