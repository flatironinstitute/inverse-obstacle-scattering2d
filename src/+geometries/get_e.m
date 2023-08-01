function [src_info] = get_e(n,iwid)
    if(n<384)
        fprintf('n must be >=384\n');
        return
    end

    if(iwid<1 || iwid > 3)
       fprintf('Invalid argument for iwid aborting\n');
       return
    end

    
    if(iwid==1), srcin = load('+geometries/letter_e1.dat'); end;
    if(iwid==2), srcin = load('+geometries/letter_e2.dat'); end;
    if(iwid==3), srcin = load('+geometries/letter_e3.dat'); end;

    srcin = srcin.';
    if(iwid==1), L = load('+geometries/letter_e1_length.dat'); end;
    if(iwid==2), L = load('+geometries/letter_e2_length.dat'); end;
    if(iwid==3), L = load('+geometries/letter_e3_length.dat'); end;
    
    
    nh = 3;
    hcoefs_use = zeros(2*nh+1,1);
    [srctmp,hout,Lout,~,~] = resample_curve(srcin,L,nh,hcoefs_use,n);
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
