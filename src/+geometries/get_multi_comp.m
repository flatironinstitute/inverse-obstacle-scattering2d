function [src_out] = get_multi_comp(src_info,centers,scales)
   src_out = [];
   [~,nc] = size(centers);
   if(nargin == 2)
       scales = ones(nc,1);
   end
   n = length(src_info.xs);
   
   nout = nc*n;
   xs = zeros(1,nout);
   ys = zeros(1,nout);
   dxs = zeros(1,nout);
   dys = zeros(1,nout);
   ds = zeros(1,nout);
   H = zeros(1,nout);
   src_out.h = src_info.h;
   for i=1:nc
       istart = (i-1)*n+1;
       iend = i*n;
       xs(istart:iend) = centers(1,i) + scales(i)*src_info.xs;
       ys(istart:iend) = centers(2,i) + scales(i)*src_info.ys;
       dxs(istart:iend) = scales(i)*src_info.dxs;
       dys(istart:iend) = scales(i)*src_info.dys;
       ds(istart:iend) = scales(i)*src_info.ds;
       H(istart:iend) = src_info.H/scales(i);
   end
   src_out.xs = xs;
   src_out.ys = ys;
   src_out.dxs = dxs;
   src_out.dys = dys;
   src_out.ds = ds;
   src_out.H = H;
end