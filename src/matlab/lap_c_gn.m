function [u] = lap_c_gn(src,targ)
%
%   This subroutine evaluates the kernel -1/2pi nx \cdot \nabla_{x} log(|x-y|)
%   x = targ(1:2), nx = targ(3:4), y = src
%  
%
  [n0,n] = size(src);
  [n1,m] = size(targ);

  src0 = src(1:2,:);
  targ0 = targ(1:4,:);

  u = complex(zeros(m,n));
  mex_id_ = 'lap_c_gn(i int[x], i double[xx], i int[x], i double[xx], io dcomplex[xx])';
[u] = lap_kernels(mex_id_, n, src0, m, targ0, u, 1, 2, n, 1, 4, m, m, n);
end
%
%
%------------------------------------------------------
