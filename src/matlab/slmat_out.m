function out = slmat_out(kh,h,src,tgt)

x_s   = src(1,:);
y_s   = src(2,:);
dx_s  = src(3,:);
dy_s  = src(4,:);

x_t   = tgt(1,:);
y_t   = tgt(2,:);

m = length(x_t);

rr = sqrt(bsxfun(@minus,x_t',x_s).^2 + bsxfun(@minus,y_t',y_s).^2);

ds = repmat(sqrt(dx_s.^2 + dy_s.^2),m,1);

out = (h)*(1i/4)*besselh(0,1,kh*rr).*ds;