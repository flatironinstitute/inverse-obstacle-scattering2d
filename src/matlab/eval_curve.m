function [binfo] = eval_curve(ts,wsave)
 %
 %  This subroutine evaluates the precomputed
 %  bandlimited curve generated using resample_curve
 %  at a collection of points in parameter space.
 %
 %  Note it is the user's responsibility to ensure
 %  that 0<ts(i)<lenght of curve, otherwise junk
 %  values would be returned for ts which do not 
 %  satisfy the criterion above
 %
 %  Input arguments:
 %    ts -
 %      coordinates where curve needs to be evaluated 
 %    wsave(1:lsave) -
 %      precomputed saved array from resample_curve
 %
 %  Output arguments:
 %    binfo(6,n) -
 %      x,y,rnx,rny,dst,curvature at requested points
 %
   tsuse = ts(:);
   n = length(tsuse);
   binfo = zeros(6,n);
   lsave = length(wsave);

   mex_id_ = 'eval_curve_multi(i int[x], i double[x], i int[x], i double[x], io double[xx])';
[binfo] = curve_resampler(mex_id_, n, tsuse, lsave, wsave, binfo, 1, n, 1, lsave, 6, n);
end
