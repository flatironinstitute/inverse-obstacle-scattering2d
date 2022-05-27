function [xmat] = comb_ext_mat(zpars,norder,h,srcinfo)
%
%  Representation:
%    u = \alpha S_{k} [\sigma] + \beta D_{k}[\sigma]
%
%  Data returned:
%    Dirichlet data (u)
%
%
%  Input: 
%    zpars(3) - 
%      zpars(1) - Helmholtz paramter
%      zpars(2) - single layer strength
%      zpars(3) - double layer strength
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
%       combined field matrix

  [m,n] = size(srcinfo);
  assert(m==6,'srcinfo must be of shape (5,n)');
  xmat = complex(zeros(n),0);
  mex_id_ = 'comb_ext_mat(i int[x], i int[x], i double[x], i double[xx], i dcomplex[x], io dcomplex[xx])';
[xmat] = kern_mats(mex_id_, n, norder, h, srcinfo, zpars, xmat, 1, 1, 1, 6, n, 3, n, n);
end
%  
%
%
%  
%
%
