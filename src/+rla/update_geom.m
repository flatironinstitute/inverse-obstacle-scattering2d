function [src_out,varargout] = update_geom(src_info,nh,hcoefs,opts)
%
%  This subroutine updates a given source discretization given an 
%  update set in fourier coeffs. It checks for self intersection, and
%  if the updated curve satsfies a certain criterion for the curvature
% 
%
%  The source curve is assumed to be equispaced discretization
%  on [0,src_info.paramL] with n nodes at 0:paramL/n:paramL*(1-1/n)
%
%  The updated curve is given by
%
%  xs_new(t) = xs(t) + (\sum_{j=0}^{nh} hcoefs(j+1)*cos(2*pi*t/L) + 
%                          \sum_{j=1}^{nh} hcoefs(j+nh+1)*sin(2*pi*t/L))
%                          nx(t)
%  ys_new(t) = ys(t) + (\sum_{j=0}^{nh} hcoefs(j+1)*cos(2*pi*t/L) + 
%                          \sum_{j=1}^{nh} hcoefs(j+nh+1)*sin(2*pi*t/L))
%                          ny(t)
%   where  nx(t), ny(t) are the x and y components of the unit normal
%    to the curve (xs(t),ys(t)), and L is the length of the original curve
%   
%  
%  If the curve is self intersecting, the routine returns the input curve
%  with ier = 1;
%
%  If eps_curv < Inf; then suppose H_{j} are the fourier coeffs of
%   the curvature of the updated curve, given n_curv = N; 
%   if (\sum_{j=-N}^{N} |H_{j}|^2)^{1/2}/(\sum_{j} |H_{j}|^2|)^{1/2} <
%   1-eps_curv; then
%  the input curve is returned with ier = 2. The above condition implies
%  that the Fourier energy of the tail of curvature of the updated curve 
%  in an arc-length parametrization is too high.
%
%  If n_curv is unspecified, but eps_curv is, then
%  n_curv = max(0.1*n,30) where L is the length of the updated
%  curve, and n is the number of points in the output curve.
%
%  Input arguments:
%   src_info - source info struct;
%      src_info.xs = x coordinates;
%      src_info.ys = y coordinates;
%      src_info.dxs = dxdt;
%      src_info.dys = dydt;
%      src_info.ds = sqrt(dxdt^2 + dydt^2);
%      src_info.h = h in trapezoidal parametrization;
%      src_info.lambda - imepdance value at discretization nodes 
%           (optional if solving impedance boundary value problem);
%   nh - max fourier content of normal update to curve
%   hcoefs(2*nh+1,1) or (1,2*nh+1) - fourier coeffs of normal update to the
%       curve;
%   opts - options struct (optional)
%      opts.eps_curv - precision for testing tail of fourier coeffs of
%         curvature of the updated curve in arc length parameterization
%         (Inf)
%      opts.n_curv - index determining the tail of the curvature test 
%         above
%
%  Output arguments:
%    src_info_out - updated source info struct
%    ier - error code (optional)
%      ier = 0, implies successful curve update
%      ier = 1, implies updated curve was self-intersecting
%      ier = 2, implies updated curve had high fourier tail in curvature
%                   of updated curve
%           


    if(nargin < 4)
        opts = [];
    end
    n = length(src_info.xs);
    paramL = src_info.paramL;
    t = 0:(2*pi/n):2*pi*(1.0-1.0/n);
    hcoefs_use = hcoefs(:);
    if(nh>0)
       h_upd = (cos(t'*(0:nh))*hcoefs_use(1:(nh+1)) + sin(t'*(1:nh))*hcoefs_use((nh+2):end)).';
    else
       h_upd = hcoefs(1);
    end
    
    x_upd = src_info.xs + h_upd.*src_info.dys./src_info.ds;
    y_upd = src_info.ys - h_upd.*src_info.dxs./src_info.ds;
    
    
    check_curv = false;
    if isfield(opts,'eps_curv')
        check_curv = true;
        eps_curv = opts.eps_curv;
        
    end

    if (~rla.issimple(x_upd,y_upd))
        src_out = src_info;
        ier = 1;
        varargout{1} = ier;
        return
    else
        dx_upd_dt = (src_info.Der_param*x_upd')';
        dy_upd_dt = (src_info.Der_param*y_upd')';
        
        ds_upd_dt = sqrt(dx_upd_dt.^2 + dy_upd_dt.^2);
        
        L = sum(ds_upd_dt)*src_info.h;
        
        
        nout = n;
        if (isfield(opts,'nppw'))
            nppw = opts.nppw;
            rlam = opts.rlam;
            Nw = L/rlam;
            nout = max(ceil(nppw*Nw),300);
            if mod(nout,2)
                nout = nout + 1;
            end
            nout = max(nout,n);
        end
        if(isfield(opts,'n_curv'))
           n_curv = opts.n_curv;
        else
           n_curv = max(floor(0.1*nout),30);
        end
        srcin = zeros(6,n);
        srcin(1,:) = src_info.xs;
        srcin(2,:) = src_info.ys;
        srcin(3,:) = src_info.dys./src_info.ds;
        srcin(4,:) = -src_info.dxs./src_info.ds;
        srcin(5,:) = src_info.ds;
        srcin(6,:) = src_info.H;
        [srctmp,hout,Lout,~,tts] = resample_curve(srcin,paramL,nh,hcoefs_use,nout);
        ier = 0;
        src_out = [];
        src_out.xs = srctmp(1,:);
        src_out.ys = srctmp(2,:);
        src_out.dxs = -srctmp(4,:);
        src_out.dys = srctmp(3,:);
        src_out.ds = srctmp(5,:);
        src_out.h = hout;
        src_out.L = Lout;
        src_out.paramL = Lout;
        rsc = 2*pi/Lout;
        src_out.Der = specdiffmat_ds(nout,src_out.ds)*rsc;
        src_out.Der_param = src_out.Der;
        src_out.H = rla.get_curvature(src_out);
        
        if(check_curv)
            curvhat = fft(src_out.H);

            curvhat_bounded = curvhat(n_curv+2:end-n_curv);        
            ratio_curv = norm(curvhat_bounded)/norm(curvhat);
            if(ratio_curv > eps_curv)
                ier = 2;
                varargout{1} = ier;
                src_out = src_info;
                return
            end
        end
        
        if isfield(src_info,'lambda')
            lambda_hat = fft(src_info.lambda);
            tt_use = tts(:)/paramL;
            
            
            kk = [(0:(n/2)) ((-n/2+1):-1)];
            src_out.lambda = sum(exp(1i*2*pi*tt_use*kk).*((lambda_hat(:)).'),2)/n;
        end
    end
    varargout{1} = ier;
    varargout{2} = tts(:);
    
    

end
