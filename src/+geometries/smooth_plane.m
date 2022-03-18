function [src_info] = smooth_plane(nterms,n)
   
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
    N_bd = 50;
    src = zeros(6,n_bd);
    src(1,:) = xs1/2.+0.5;
    src(2,:) = ys1/2.;
    src(3,:) = ones(1,n_bd);
    src(4,:) = ones(1,n_bd);
    src(5,:) = ones(1,n_bd);
    
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

    
    h = Lplane/n;
    tplane = 0:h:(Lplane-h);

    xs  = fnval(tplane,xplane);
    ys  = fnval(tplane,yplane);
    dxs = fnval(tplane,dxplane);
    dys = fnval(tplane,dyplane);
    ds  = sqrt(dxs.^2 + dys.^2);
    rnx = dys./ds;
    rny = -dxs./ds;

    srcin = zeros(6,n);
    srcin(1,:) = xs;
    srcin(2,:) = ys;
    srcin(3,:) = rnx;
    srcin(4,:) = rny;
    srcin(5,:) = ds;
    
    
    nh = 3;
    hcoefs_use = zeros(2*nh+1,1);
    [srctmp,hout,Lout,~,~] = resample_curve(srcin,Lplane,nh,hcoefs_use,n);
    src_info = [];
    src_info.xs = srctmp(1,:);
    src_info.ys = srctmp(2,:);
    src_info.dxs = -srctmp(4,:);
    src_info.dys = srctmp(3,:);
    src_info.ds = srctmp(5,:);
    src_info.h = hout;
    src_info.L = Lout;
    src_info.paramL = Lout;
    rsc = 2*pi/Lout;
    src_info.Der = specdiffmat_ds(n,src_info.ds)*rsc;
    src_info.Der_param = src_info.Der;
    src_info.H = rla.get_curvature(src_info);

    
end