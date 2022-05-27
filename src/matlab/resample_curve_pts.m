function [varargout] = resample_curve_pts(xys,nb,eps,nuse)
%
%  This subroutine resamples a given curve through a collection
%  of points xys and attempts to pass a bandlimited curve through
%  them.
%
%  Input arguments:
%    xys(2,n) -
%      x,y coordinates of the input points
%    nb - band limit of output curve
%    eps - 
%      (optional) desired accuracy for interpolation of
%      output curve thorugh input points. Default value is
%      1e-13.
%    nuse - 
%      max number of doublings in sampling points to achieve
%      desired accuracy
%      
%  
%  Output arguemnts:
%    rltot - Length of the curve
%    wsave(lsave) - array which stores info for resampling points
%      later at desired values
%
  [~,n] = size(xys);
  if( nargout < 2 || nargout > 3)
    fprintf('invalid number of output arguments\n');
    fprintf('out arguments must be 2,3\n');
    varargout(1:nargout) = {0};
    return;
  end

  if (nargin == 2)
    epsuse = 1e-13;
    nmax = 3;
  elseif (nargin == 3)
    epsuse = eps;
    nmax = 3;
  elseif (nargin == 4)
    epsuse = eps;
    nmax = nuse;
  else
    fprintf('invalid number of arguments\n');
  end

  nlarge = 0;
  nout = 0;
  lsave = 0;
  lused = 0;
  ierm = 0;
  mex_id_ = 'simple_curve_resampler_mem(i int[x], i double[xx], i int[x], i double[x], i int[x], io int[x], io int[x], io int[x], io int[x], io int[x])';
[nlarge, nout, lsave, lused, ierm] = curve_resampler(mex_id_, n, xys, nb, epsuse, nmax, nlarge, nout, lsave, lused, ierm, 1, 2, n, 1, 1, 1, 1, 1, 1, 1, 1);

  if(ierm == 4)
    fprintf('nb too small resulting curve self intersecting');
    fprintf('change nb or eps');
    fprintf('returning without resampling curve');
    varargout(1:nargout) = {0};
  else if(ierm == 2)
    fprintf('warning: desired interpolation accurcy not reached');
    fprintf('try changing nb');
  end

  wsave = zeros(lsave,1);
  nnn = n+1;
  tts = zeros(nnn,1);
  sinfo = zeros(6,nout);
  hout = 0.0;
  rltot = 0.0;
  mex_id_ = 'simple_curve_resampler_guru(i int[x], i double[xx], i int[x], i int[x], i int[x], i int[x], i int[x], io double[xx], io double[x], io double[x], io double[x], io double[x], io int[x])';
[sinfo, hout, rltot, wsave, tts, ierm] = curve_resampler(mex_id_, n, xys, nb, nlarge, lsave, lused, nout, sinfo, hout, rltot, wsave, tts, ierm, 1, 2, n, 1, 1, 1, 1, 1, 6, nout, 1, 1, lsave, nnn, 1);
  
  varargout{1} = rltot;
  varargout{2} = wsave;
  if(nargout>=3)
    varargout{3} = tts(1:n);
  end

end
%
%
%
%
