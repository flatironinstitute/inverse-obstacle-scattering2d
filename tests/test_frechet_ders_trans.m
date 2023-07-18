clear

run ../inverse-obstacle-scattering2d/startup.m

% define geometry type
n = 300;
nc = 2; 
coefs = zeros(1,2*nc+1);
coefs(1) = 1.2;
coefs(4) = 0.3;
[src_info] = geometries.starn(coefs,nc,n);


%ratio type for contrast
q_plus_1 = 1.1/1.7;

%domain jump
rho_d = 1.9;

%transmission parameters
a = [complex(1.0); complex(1.0)];
b = [complex(1.0); complex(rho_d)];

% define k0 (starting frequncy, dk, spacing in frequency and 
% number of frequencies (nk)
% Set of frequencies (k_{i})
k0 = 1.7;
kmax = 1.7;
dkinv = 1;
dk = 1.0/dkinv;
kh = k0:dk:kmax;
nk = length(kh);
kh = 1.7;

% setting incident waves and receptors
% inc_type = 0, 1 inc direction with 6 receptors
% inc_type = 1, 10 inc direction with 200 receptors
% inc_type = 2, 50 inc directions with 200 receptors
% inc_type = 3, 100 inc directions with 200 receptors
% inc_type = 4  2*k incident directions with 8*k receptors
% inc_type = 5, 10*k incident directions with 10*k receptors
inc_type = 0;


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

%the first id for the cavity the second for charlie's
src0 = [0.0;0.0];

opts.test_analytic = true;
opts.src_in = src0;
opts.verbose=true;


%Calculating f(q)
n  = 300;

src_info = geometries.starn(coefs,nc,n);
% [src_info] = geometries.ellipse(1,1,n);

L = src_info.L;

u_meas = cell(nk,1);

nppw = 30;

for ik=1:nk
    if(inc_type == 0)
        n_dir = 1;
        n_tgt = 5;
    end
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

    %Calculating F(q)
    fprintf('F(q)\n')
    n = ceil(nppw*L*abs(kh(ik))/2/pi);
    n = max(n,408);
    src_info = geometries.starn(coefs,nc,n);
%     [src_info] = geometries.ellipse(1,1,n);

    %options of curve update
    opts_update_geom = [];
    eps_curv =1e-1;
    opts_update_geom.eps_curv = eps_curv;
    ncoeff_boundary = floor(2*abs(kh(ik)));
    n_curv_min = 20;
    opts_update_geom.n_curv = max(n_curv_min,ncoeff_boundary);
    rlam = 2*pi/real(kh(ik));    
    opts_update_geom.nppw = nppw;
    opts_update_geom.rlam = rlam;

    %re-parameterizing
    coefs_h=zeros(1,2*ncoeff_boundary+1);
    [src_out,ier,a,xu,yu] = rla.update_geom(src_info,ncoeff_boundary,coefs_h,opts_update_geom);
    src_info = src_out;
    
    %checking 
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
   src_orig = src_info;   
   
   %getting erivative
   frechet_mats = rla.get_frechet_ders(kh(ik),mats,src_info,sensor_info,fields,bc,opts);
   Minv = frechet_mats.bdry;  
   
   for ivar = 1 : 2*ncoeff_boundary+1
       fprintf('ivar=%d\n',ivar)
       coefs_h = zeros(1,2*ncoeff_boundary+1);       
        for jj = 1 : 6
            coefs_h(ivar) = .1.^jj;
       
           %Calculating F(q+dq)
           fprintf('F(q+dq)\n')   
           fprintf('coefs_h=%d\n',coefs_h(ivar))
           
           [src_out,ier,a,xu,yu] = rla.update_geom(src_info,ncoeff_boundary,coefs_h,opts_update_geom);
           [mats_out] = rla.get_fw_mats(kh(ik),src_out,bc,sensor_info,opts);
           fields_dq = rla.compute_fields(kh(ik),src_out,mats_out,sensor_info,bc,opts);
                      
           field_approx = fields_dq.uscat_tgt - (fields.uscat_tgt+Minv*(coefs_h'));
           
           fprintf('Error=%d\n',max(abs(field_approx(:))))                      

        end
        
        pause
        
   end
   
   
end

    
    
