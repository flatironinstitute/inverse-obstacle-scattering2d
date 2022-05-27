function [varargout] = resample_curve(srcinfo,rl,nh,hcoefs,nout,eps)
%
%  This subroutine resamples a given curve whose xy 
%  and attempts to pass a bandlimited curve thorugh them.
%
%  Makes no assumption about the input being equispaced to begin with
%
%  Input arguments:
%    srcinfo(6,n) -
%      srcinfo(1:2,:) - x,y coordinates of the input points
%      srcinfo(3:4,:) - x and y component of normals
%      srcinfo(5,:) - dst at input points 
%      srcinfo(6,:) curv at input points
%    h - 
%       Length of parameter space/n
%    rl - The input curve is assumed to be a map from
%       [0,rl] \to R^{2}
%    nh - number of non-zero coeffs in h
%    hcoefs - the coeffs
%       h = \sum_{j=0}^{nh} c_{j} \cos(2*pi*t/rl) +
%          \sum_{j=1}^{nh} s_{j} \sin(2*pi*t/rl)
%    where hcoefs(1:nh+1) = c_{j} and hcoefs(nh+2:2*nh+1) = s_{j}
%    nout -
%      (optional) number of points at which output curve
%      should be sampled
%    eps - 
%      (optional) desired accuracy for interpolation of
%      output curve thorugh input points. Default value is
%      1e-13.
%  
%  Output arguemnts:
%    srcinfoout(6,nout) - 
%      resampled curve
%    hout - 
%      spacing between points in parameter space
%    rltot -
%      (optional, if requested) Length of the curve
%    ier - 
%       (optional, if requested) error code for 
%        obtainining arc length parametrization
%    tts(1:nout) -
%       (optional, if requested) location between
%        [0,rl] corresponding to equispaced nodes
%
  [~,n] = size(srcinfo);
  if( nargout < 2 || nargout > 5)
    fprintf('invalid number of output arguments\n');
    fprintf('out arguments must be 3,4,5\n');
    varargout(1:nargout) = {0};
    return;
  end

  if (nargin == 4)
    epsuse = 1e-13;
    nuse = n + 4*nh;
  elseif (nargin == 5)
    epsuse = 1e-13;
    nuse = nout;
  elseif (nargin == 6)
    epsuse = eps;
    nuse = nout;
  else
    fprintf('invalid number of arguments\n');
  end
  xfft = fft(srcinfo(1,1:n))/n;
  yfft = fft(srcinfo(2,1:n))/n;
  dxdt = -srcinfo(4,:).*srcinfo(5,:);
  dydt = srcinfo(3,:).*srcinfo(5,:);
  dxfft = fft(dxdt)/n;
  dyfft = fft(dydt)/n;

  nx = n;
  nhuse = 2*nh+1;
  npar1 = 4*nx + nhuse;
  a = zeros(npar1,1);
  par1 = complex(a,0);
  par1(1:n) = xfft + 1j*yfft;
  par1((n+1):(2*n)) = dxfft + 1j*dyfft;

 
  hfft = zeros(nhuse,1);
  hfft(1) = hcoefs(1);
  hfft(2:(nh+1)) = hcoefs(2:(nh+1))/2;
  hfft((nh+2):nhuse) = flip(hcoefs(2:(nh+1)))/2;
  hfft(2:(nh+1)) = hfft(2:(nh+1))+hcoefs((nh+2):nhuse)/2/1j;
  hfft((nh+2):nhuse) = hfft((nh+2):nhuse)-flip(hcoefs((nh+2):nhuse))/2/1j;
  dhfft = zeros(nhuse,1);
  dhfft(2:(nh+1)) = hfft(2:(nh+1))*1j.*(1:nh)'*2*pi/rl;
  dhfft((nh+2):nhuse) = -hfft((nh+2):nhuse)*1j.*flip(1:nh)'*2*pi/rl;
  par1((2*n+1):(2*n+nhuse)) = hfft;
  par1((2*n+nhuse+1):(2*(n+nhuse))) = dhfft;

  
  nn = nuse+1;
  lw = 10000;
  work = zeros(lw,1);
  tts = zeros(nn,1);
  ier = 0;
  sinfo = zeros(4,nuse);
  hout = 0.0;
  rltot = 0.0;
  lsave = 0;
  
  mex_id_ = 'curve_resampler_guru(io int[x], i int[x], i int[x], i dcomplex[x], i double[x], i int[x], i double[x], io double[x], io double[xx], io double[x], io double[x], io double[x], i int[x], io int[x])';
[ier, tts, sinfo, hout, rltot, work, lsave] = curve_resampler(mex_id_, ier, n, nhuse, par1, rl, nuse, epsuse, tts, sinfo, hout, rltot, work, lw, lsave, 1, 1, 1, npar1, 1, 1, 1, nn, 4, nuse, 1, 1, lw, 1, 1);
  
  srcinfoout = zeros(6,nuse);
  srcinfoout(1,:) = sinfo(1,:);
  srcinfoout(2,:) = sinfo(2,:);
  srcinfoout(3,:) = sinfo(4,:);
  srcinfoout(4,:) = -sinfo(3,:);
  srcinfoout(5,:) = ones(1,nuse);
  
  varargout{1} = srcinfoout;
  varargout{2} = hout;
  if(nargout>=3)
     varargout{3} = rltot;
  end

  if(nargout>=4)
     varargout{4} = ier;
  end

  if(nargout>=5)
     varargout{5} = tts(1:nuse);
  end
  
end
%
%
%
%
%
%
