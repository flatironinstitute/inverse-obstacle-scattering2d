% clear
% addpath('./src')
BW = imread('SU571.png');
BW = im2bw(BW, graythresh(BW));
[B,L] = bwboundaries(BW);%maybe noholes option
x1=B{2};%or B{1}
xs=x1(:,1)'/100-2;
ys=x1(:,2)'/100-3.5;
xs1=xs(end:-1:1);
ys1=ys(end:-1:1);
A=[0 1;-1 0];
b=A*[xs1;ys1];
xs1=b(1,:);
ys1=b(2,:);

%calculating L
difx = diff(xs1);
dify = diff(ys1);
L    = sum(sqrt(difx.^2+dify.^2));

n_bd = length(xs1);
N_bd =50;
src = zeros(5,n_bd);
src(1,:) = xs1/2.+0.5;
src(2,:) = ys1/2.;
src(3,:) = ones(1,n_bd);
src(4,:) = ones(1,n_bd);
src(5,:) = ones(1,n_bd);
nterms   = 20;
[srcout,~,Lplane,~,~] = resample_curve(src,L,N_bd,zeros(2*N_bd+1,1),nterms);
h = Lplane/nterms;
tplane=0:h:Lplane;
xplane = csape(tplane,[srcout(1,:) srcout(1,1)],'periodic');
yplane = csape(tplane,[srcout(2,:) srcout(2,1)],'periodic');
dxplane.form   = xplane.form;
dyplane.form   = yplane.form;
dxplane.breaks = xplane.breaks;
dyplane.breaks = yplane.breaks;
dxplane.pieces = xplane.pieces;
dyplane.pieces = yplane.pieces;
dxplane.order  = xplane.order-1;
dyplane.order  = yplane.order-1;
dxplane.dim    = xplane.dim;
dyplane.dim    = yplane.dim;
dxplane.coefs  = bsxfun(@times,xplane.coefs(:,1:end-1),[3 2 1]); 
dyplane.coefs  = bsxfun(@times,yplane.coefs(:,1:end-1),[3 2 1]); 

N = 1000;
h = Lplane/N;
tplane = 0:h:Lplane;

xs  = fnval(tplane,xplane);
ys  = fnval(tplane,yplane);
dxs = fnval(tplane,dxplane);
dys = fnval(tplane,dyplane);
ds  = sqrt(dxs.^2 + dys.^2);
rnx = dys./ds;
rny = -dxs./ds;
bd_ref = struct('t_bd',tplane,'h_bd',h,'xs',xs,'ys',ys,'dxs',dxs,'dys',...
    dys,'ds',ds,'Lplane',L,'nterms',nterms);

%plot(xs,ys);hold on; quiver(xs,ys,rnx,rny)

% save('plane.mat','xplane','yplane','dxplane','dyplane','Lplane')

