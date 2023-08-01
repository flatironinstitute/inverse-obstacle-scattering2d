src_info = geometries.ellipse(1.1,1.1,600);
nh = 4;
hcoefs = zeros(1,2*nh+1);
hcoefs(nh+1) = 0.3;
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

noveruse = 5;

nn = numel(x_upd);
x_ft = fft(x_upd);
y_ft = fft(y_upd);
x_big_ft = zeros([noveruse*nn,1]);
y_big_ft = zeros([noveruse*nn,1]);
x_big_ft(1:(nn/2)) = x_ft(1:(nn/2));
y_big_ft(1:(nn/2)) = y_ft(1:(nn/2));
x_big_ft(noveruse*nn -nn/2 +1:end) = x_ft((nn/2+1):end);
y_big_ft(noveruse*nn -nn/2 +1:end) = y_ft((nn/2+1):end);
xbig = real(ifft(x_big_ft))*(noveruse);
ybig = real(ifft(y_big_ft))*(noveruse);
%plot(xs,ys,'b.',xbig,ybig,'r.')

P = polyshape(xbig,ybig,'KeepCollinearPoints',true);