function [ux,uy] = helm_c_g(zk,src,targ)
%
%   This subroutine evaluates the kernel i/4 nx \cdot \nabla_{x} H_{0}(k|x-y|)
%   x = targ(1:2), nx = targ(3:4), y = src
%  
%
  [n0,n] = size(src);
  [n1,m] = size(targ);

  src0 = src(1:2,:);
  targ0 = targ(1:2,:);

  ux = complex(zeros(m,n));
  uy = complex(zeros(m,n));
  mex_id_ = 'helm_c_g(i int[x], i double[xx], i int[x], i double[xx], i dcomplex[x], io dcomplex[xx], io dcomplex[xx])';
[ux, uy] = helm_kernels(mex_id_, n, src0, m, targ0, zk, ux, uy, 1, 2, n, 1, 2, m, 1, m, n, m, n);
end
%
%------------------------------------------------------
