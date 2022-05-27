function [u] = helm_c_p(zk,src,targ)
%
%     
%   This suborutine evaluates the kernel i/4 H_{0} (k|x-y|) 
%     x = targ, y = src
%
  [n0,n] = size(src);
  [n1,m] = size(targ);

  src0 = src(1:2,:);
  targ0 = targ(1:2,:);

  u = complex(zeros(m,n));
  mex_id_ = 'helm_c_p(i int[x], i double[xx], i int[x], i double[xx], i dcomplex[x], io dcomplex[xx])';
[u] = helm_kernels(mex_id_, n, src0, m, targ0, zk, u, 1, 2, n, 1, 2, m, 1, m, n);
end
%
%
%------------------------------------------------------
