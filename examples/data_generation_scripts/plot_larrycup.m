n = 500;
m = 200;

src_info = geometries.larrycup(0.2,pi/12,n,m);

nh = 0;
hcoefs = zeros(1,1);
[src_out] = rla.update_geom(src_info,nh,hcoefs);

plot(src_out.xs,src_out.ys,'k.');
axis equal

L = src_info.L;
