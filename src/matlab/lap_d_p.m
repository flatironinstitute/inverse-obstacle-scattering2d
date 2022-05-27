function [u] = lap_d_p(src,targ)
%
%   This subroutine evaluates the kernel -1/2\pi ny \cdot \nabla_{y} log(|x-y|)
%   x = targ, y = src(1:2), ny = src(3:4)
%
  [n0,n] = size(src);
  [n1,m] = size(targ);

  src0 = src(1:4,:);
  targ0 = targ(1:2,:);

  u = complex(zeros(m,n));
  mex_id_ = 'lap_d_p(i int[x], i double[xx], i int[x], i double[xx], io dcomplex[xx])';
[u] = lap_kernels(mex_id_, n, src0, m, targ0, u, 1, 4, n, 1, 2, m, m, n);
end
