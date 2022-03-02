function out = dlmat_out(kh,h,src,tgt)

x_s   = src(1,:);
y_s   = src(2,:);
dx_s  = src(3,:);
dy_s  = src(4,:);

x_t   = tgt(1,:);
y_t   = tgt(2,:);

m = length(x_t);

rr = sqrt(bsxfun(@minus,x_t',x_s).^2 + bsxfun(@minus,y_t',y_s).^2);

out = (h)*(1i*kh/4) * ...
    (bsxfun(@minus,x_t',x_s).*repmat(dy_s,m,1)-bsxfun(@minus,y_t',y_s).*repmat(dx_s,m,1)) .* ...
    besselh(1,1,kh*rr)./rr;