function ipass = test_update_geom()
ipass = 1;
ifplot = false;


n  = 600;
a = 1.1;  
b = 1.1;  
src_info = geometries.ellipse(a,b,n);

if(ifplot)
    figure
    plot(src_info.xs,src_info.ys,'b.');
    axis equal
    drawnow;
end

lambda_fun = @(t) sin(t + 0.3) + 0.1*cos(3*t) + 0.3;
%lambda_fun = @(t) sin(t + 0.3);
n = length(src_info.xs);
ts = (0:(n-1))/n*2*pi;
src_info.lambda = lambda_fun(ts).';

if(ifplot)
    figure
    clf
    plot(ts',src_info.lambda,'k.');
    drawnow;
end

nh = 0;
hcoefs = 0;
[src_out,ier,tts] = rla.update_geom(src_info,nh,hcoefs);

if(ifplot)
    figure
    clf
    plot(src_out.xs,src_out.ys,'k.');
    axis equal
    drawnow
end


% test interpolated lambda
lambda_out = src_out.lambda;
xs_out = src_out.xs/a;
ys_out = src_out.ys/b;
thet = atan2(ys_out,xs_out);
lambda_test = lambda_fun(thet).';

err_lam = norm(lambda_out(:) - lambda_test(:))/norm(lambda_test(:));
if(err_lam>1e-15) 
    ipass = 0;
    fprintf('failed lambda interpolation in update_geom test: %d\n',err_lam);
end




% Test updated curve with
nh = 4;
hcoefs = zeros(1,2*nh+1);
hcoefs(nh+1) = 0.3;

src_info = rmfield(src_info,'lambda');


[src_out,ier] = rla.update_geom(src_info,nh,hcoefs);

if(ifplot)
    figure
    clf
    plot(src_out.xs,src_out.ys,'b.')
    axis equal
    drawnow
end

% test all quantities associated with resampled curve
xs = src_out.xs;
ys = src_out.ys;
rs = sqrt(xs.^2 + ys.^2);
thet = atan2(ys./rs,xs./rs);

rfun = @(t) a + hcoefs(nh+1)*cos(nh*t);
drfun = @(t) -nh*hcoefs(nh+1)*sin(nh*t);
d2rfun = @(t) -nh*nh*hcoefs(nh+1)*cos(nh*t);

xsout = rfun(thet).*cos(thet);
ysout = rfun(thet).*sin(thet);
dxsout = drfun(thet).*cos(thet) - rfun(thet).*sin(thet);
dysout = drfun(thet).*sin(thet) + rfun(thet).*cos(thet);

dsout = sqrt(dxsout.^2 + dysout.^2);
dxsouttest = dxsout./dsout;
dysouttest = dysout./dsout;

dstest = sqrt(src_out.dxs.^2 + src_out.dys.^2);

dxstest = src_out.dxs./dstest;
dystest = src_out.dys./dstest;

errdxs = norm(dxsouttest - dxstest)/norm(dxsouttest);
errdys = norm(dysouttest - dystest)/norm(dysouttest);


d2xsout = d2rfun(thet).*cos(thet) - 2*drfun(thet).*sin(thet) - rfun(thet).*cos(thet);
d2ysout = d2rfun(thet).*sin(thet) + 2*drfun(thet).*cos(thet) - rfun(thet).*sin(thet);

curv_out = (dxsout.*d2ysout - dysout.*d2xsout)./dsout.^3;

errcurv = norm(curv_out-src_out.H)/norm(src_out.H);
if(errdxs>1e-14) 
    ipass = 0;
    fprintf('failed dx/dt test in update_geom test: %d\n',errdxs);
end

if(errdys>1e-14) 
    ipass = 0;
    fprintf('failed dy/dt test in update_geom test: %d\n',errdys);
end

if(errcurv>1e-8) 
    ipass = 0;
    fprintf('failed curvature test in update_geom test: %d\n',errcurv);
end

end
