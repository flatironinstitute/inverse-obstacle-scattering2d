% This script generates the data for a star shaped domain, for a fixed
% number of sensors and incident directions where data is available for all
% sensors at each incident direction
clear

run ../inverse-obstacle-scattering2d/startup.m
%addpath('/Users/borges/Desktop/Research/Current_Projects/Jeremy/rla-monograph-tests/cavity-tests')

warning('off');
ifgenerate = 0; % flag for generating the required data

% define geometry type
charlie_flag = 10.2;
dom_type = 10;

%ratio type for contrast
q_plus_1 = 1.5;

%domain jump
rho_d = 1;

%transmission parameters
a = [complex(1.0); complex(1.0)];
b = [complex(1.0); complex(rho_d)];


% define k0 (starting frequncy, dk, spacing in frequency and 
% number of frequencies (nk)
k0 = 1;
dkinv = 4;
dk = 1.0/dkinv;
%nk = 117;

khmax = 30;
nk = (khmax-1)*dkinv+1;


% setting incident waves and receptors
% inc_type = 1, 10 inc direction with 200 receptors
% inc_type = 2, 50 inc directions with 200 receptors
% inc_type = 3, 100 inc directions with 200 receptors
% inc_type = 4  2*k incident directions with 8*k receptors
% inc_type = 5, 10*k incident directions with 10*k receptors
inc_type = 5;


% choice of noise levels
% noise_type = 0; no noise
% noise_type = 1; additive
% noise_type = 2; multiplicative
noise_type = 0;
noise_lvl = 0.02;

% Boundary condition parameters
bc = [];
bc.type = 'Transmission';
bc.invtype = 'o';
% bc.transk = %this changes at every frequency
bc.q = q_plus_1;
bc.transa = a;
bc.transb = b;


% optimization parameters
optim_opts = [];
opts = [];
opts.verbose=false;%true
optim_opts.optim_type = 'sd';
optim_opts.filter_type = 'gauss-conv';
opts.store_src_info = true;
optim_opts.maxit = 100;
opts.use_lscaled_modes = true;
optim_opts.eps_upd = 1e-3;

ifcons = 1;
if(ifcons)
    optim_opts.eps_curv = 1e-1;
end

optim_opts.n_curv_min = 20;


% Data and solution directories
dir_data = './';%'./charlie-data/';
dir_sol = './';%'./charlie-sol/';
dir_diary = './';%'./charlie-diary/';


fname = [dir_data 'dom' num2str(dom_type) '_ik' num2str(k0) '_nk' int2str(nk) '_dk' ...
     num2str(dk) '_cflag' num2str(charlie_flag) '_contrast' num2str(q_plus_1) '_inctype' ...
     int2str(inc_type) ...
     '_noise' int2str(noise_type) 'noise_lvl' num2str(noise_lvl) ... 
     '_data_' bc.type '.mat'];

fname_sol = [dir_sol 'dom' num2str(dom_type) '_ik' num2str(k0) '_nk' int2str(nk) '_dk' ...
     num2str(dk) '_cflag' num2str(charlie_flag) '_contrast' num2str(q_plus_1) '_inctype' ...
     int2str(inc_type) ...
     '_noise' int2str(noise_type) 'noise_lvl' num2str(noise_lvl) ... 
     '_data_' bc.type '_optimtype_' optim_opts.optim_type '_filtertype_' ...
     optim_opts.filter_type '_ifcons' int2str(ifcons) '_ncurvmin' ...
     int2str(optim_opts.n_curv_min) '_epscurv' num2str(optim_opts.eps_curv) ... 
     '_lscaled.mat'];

fname_diary = [dir_diary 'dom' num2str(dom_type) '_ik' num2str(k0) '_nk' int2str(nk) '_dk' ...
     num2str(dk) '_cflag' num2str(charlie_flag) '_contrast' num2str(q_plus_1) '_inctype' ...
     int2str(inc_type) ...
     '_noise' int2str(noise_type) 'noise_lvl' num2str(noise_lvl) ... 
     '_data_' bc.type '_optimtype_' optim_opts.optim_type '_filtertype_' ...
     optim_opts.filter_type '_ifcons' int2str(ifcons) '_ncurvmin' ...
     int2str(optim_opts.n_curv_min) '_epscurv' num2str(optim_opts.eps_curv) ...
     '_lscaled.mat'];

diary(fname_diary);
 
%the first id for the cavity the second for charlie's
% src0 = [0.01;-1.1];
%flag 4
src0 = [0.85;0.0];
%flag 10.2 
% src0 = [0.6;0.0];

opts.test_analytic = true;
opts.src_in = src0;
opts.verbose=true;

% Set of frequencies (k_{i})
kh = 1:dk:(1+(nk-1)*dk);

if (ifgenerate)
     n  = 300;

    src_info = geometries.charlie_cavity(charlie_flag,n);    
    
    L = src_info.L;

    save(fname,'src_info');
    u_meas = cell(nk,1);

    nppw = 30;

    for ik=1:nk
        if(inc_type == 1)
            n_dir = 10;
            n_tgt = 200;
        end
        if(inc_type == 2)
            n_dir = 50;
            n_tgt = 200;
        end
        if(inc_type == 3)
            n_dir = 100;
            n_tgt = 200;
        end
        if(inc_type == 4)
            n_dir = ceil(2*kh(ik));
            n_tgt = ceil(8*kh(ik));
        end
        
        if(inc_type == 5)
            n_dir = ceil(10*kh(ik));
            n_tgt = ceil(10*kh(ik));
        end

        % set target locations
        %receptors (r_{\ell})
        r_tgt = 10;
        t_tgt = 0:2*pi/n_tgt:2*pi-2*pi/n_tgt;

        % Incident directions (d_{j})
        t_dir = 0:2*pi/n_dir:2*pi-2*pi/n_dir;

        [t_tgt_grid,t_dir_grid] = meshgrid(t_tgt,t_dir);
        t_tgt_grid = t_tgt_grid(:);
        t_dir_grid = t_dir_grid(:);
        xtgt = r_tgt*cos(t_tgt_grid);
        ytgt = r_tgt*sin(t_tgt_grid);
        tgt   = [ xtgt'; ytgt'];


        sensor_info = [];
        sensor_info.tgt = tgt;
        sensor_info.t_dir = t_dir_grid;


        n = ceil(nppw*L*abs(kh(ik))/2/pi);
        n = max(n,500);
        src_info = geometries.charlie_cavity(charlie_flag,n);         
        khi = kh(ik)*bc.q;        
        bc.transk = [khi kh(ik)];

       [mats,erra] = rla.get_fw_mats(kh(ik),src_info,bc,sensor_info,opts);
       fields = rla.compute_fields(kh(ik),src_info,mats,sensor_info,bc,opts);

       u_meas0 = [];
       u_meas0.kh = kh(ik);
       u_meas0.uscat_tgt = fields.uscat_tgt;
       u_meas0.tgt = sensor_info.tgt;
       u_meas0.t_dir = sensor_info.t_dir;
       u_meas0.err_est = erra;
       u_meas{ik} = u_meas0;
       fprintf('kh = %d    err= %d\n',kh(ik),u_meas0.err_est);
    end


    save(fname,'u_meas','-append','-v7.3');
    return
    
else
    S = load(fname);
    u_meas = S.u_meas;
end

% stop

%%%%%%%%%%%%%%%
diary off

niter = 1;

allxs = {};
allys = {};
errs  = [];
reso  = [];
resa  = [];
edict = {};
ik_lists = {};
src_init=[];

path_flag = 0.5;

% fname_solution = [dir_data 'cf' num2str(charlie_flag) '_result.mat'];
fname_result = [dir_data 'charlie_ik' num2str(k0) '_nk' int2str(nk) '_dk' ...
     num2str(dk) '_cflag' num2str(charlie_flag) '_contrast' num2str(q_plus_1) '_inctype' ...
     int2str(inc_type) ...
     '_noise' int2str(noise_type) 'noise_lvl' num2str(noise_lvl) ... 
     '_data_' bc.type '_optimtype_' optim_opts.optim_type '_filtertype_' ...
     optim_opts.filter_type '_ifcons' int2str(ifcons) '_ncurvmin' ...
     int2str(optim_opts.n_curv_min) '_epscurv' num2str(optim_opts.eps_curv) ...
     '_paths' int2str(niter) '_pflag' num2str(path_flag)  '_resultv3.mat'];

for ijk=1:niter
    
fprintf('Iteration=%d\n',ijk)

lam_init=[];
rla_path_opts = [];
llen = length(u_meas);

if (ijk >1)
iii = (rand(10000,1)).^(path_flag);
iii = max(cumsum(2*round(iii)-1),1);
imaxs = find(iii==llen);
imax = imaxs(1);
%nbtrack = 150;
rla_path_opts.ik_list = iii(1:imax);
else
    %the first track is the direct path
rla_path_opts.ik_list = 1:llen;
end

% start inverse problem
%fac = 0.25*(2^(ijk-1));
%fac = 0.25;
fac = 1;
src_init=[];
%for cf<11
% src_init = geometries.ellipse(0.15,0.15,500);
% src_init.xs = src_init.xs + 0.85;
src_init = geometries.ellipse(1,1,500);
%for cf=11
% src_init = geometries.ellipse(0.6,0.6,500);

tic; [inv_data_all,src_info_out] = rla.rla_inverse_solver(u_meas,bc,...
                          optim_opts,opts,src_init,lam_init,rla_path_opts,fac); toc;
src_init = src_info_out;

nlast = size(inv_data_all);
%allxs{ijk} = inv_data_all{nlast}.src_info_all{1}.xs;
%allys{ijk} = inv_data_all{nlast}.src_info_all{1}.ys;
allxs{ijk} = src_info_out.xs;
allys{ijk} = src_info_out.ys;
ik_lists{ijk} = rla_path_opts.ik_list;

xtmp = allxs{ijk};
ytmp = allys{ijk};

S.src_info
xt = S.src_info.xs;
yt = S.src_info.ys;

% [u_meas2] = re_eval_fields(src_info_out,nk,kh,inc_type,opts,bc);
% [evec] =compute_mismatch(u_meas,u_meas2);
% edict{ijk} = evec;

polytrue = polyshape(xt,yt);
polythis = polyshape(xtmp,ytmp);

pdiff1 = subtract(polytrue,polythis);
pdiff2 = subtract(polythis,polytrue);

err = area(pdiff1)+area(pdiff2);

errs = [errs,err];
reso = [reso,inv_data_all{end}.res_opt];
%resa = [resa,inv_data_all{end}.res_all];
% save(fname_result,'inv_data_all','allxs','allys','errs','ik_lists','reso','resa');
save(fname_result,'allxs','allys','errs','ik_lists','reso','resa');

end                      
fprintf('Saving results\n')
save(fname_result,'-v7.3')
fprintf('Saved results\n')
