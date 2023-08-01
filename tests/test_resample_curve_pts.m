clear;
clc;
addpath('../src');

x = importdata('bunny10out');
x = x';
nb = 150;
% simple interface
tic, [rltot, wsave] = resample_curve_pts(x,nb); toc;

fprintf('Length of curve=%d\n',rltot);

% Full interface
eps = 1e-9;
tic, [rltot, wsave, tts] = resample_curve_pts(x,nb,eps,4); toc;


% Evalute curve at user specified points
tic, binfo = eval_curve(tts,wsave); toc;

figure(1);
clf();
plot(x(1,:),x(2,:),'rx'), hold on;
plot(binfo(1,:),binfo(2,:),'k-');

% Routine to resample curve at a given number of nodes
% currently being used to test interpolation error
err = norm(binfo(1,:) - x(1,:)).^2 + norm(binfo(2,:)-x(2,:)).^2;
rr = norm(x(:))^2;
err = sqrt(err/rr);
fprintf('error in interpolation=%d\n',err);


%% convert to src struct format

[~,n] = size(binfo);

srcin = zeros(6,n);
srcin(1,:) = binfo(1,:);
srcin(2,:) = binfo(2,:);
srcin(3,:) = binfo(3,:);
srcin(4,:) = binfo(4,:);
srcin(5,:) = binfo(5,:);
srcin(6,:) = 0;

nout = n;
nh = 0;
hcoefs = zeros(3,1);
[srctmp,hout,Lout,~,tts] = resample_curve(srcin,rltot,nh,hcoefs,nout);
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

