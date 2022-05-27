function [u] = helm_d_gn(zk,src,targ)
%
%   This subroutine evaluates the kernel i/4 ny nx \cdot \nabla_{y} \nabla_{x} H_{0}(k|x-y|)
%   x = targ, y = src(1:2), ny = src(3:4)
%
  [n0,n] = size(src);
  [n1,m] = size(targ);

  src0 = src(1:4,:);
  targ0 = targ(1:4,:);

  u = complex(zeros(m,n));
  mex_id_ = 'helm_d_gn(i int[x], i double[xx], i int[x], i double[xx], i dcomplex[x], io dcomplex[xx])';
[u] = helm_kernels(mex_id_, n, src0, m, targ0, zk, u, 1, 4, n, 1, 4, m, 1, m, n);
end
