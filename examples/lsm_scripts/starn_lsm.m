% This script generates the data for a star shaped domain, for a fixed
% number of sensors and incident directions where data is available for all
% sensors at each incident direction


n  = 300;
addpath('../');


% max number of wiggles
nc = 3;

% parameters 'a'
coefs = zeros(2*nc+1,1);
coefs(1) = 1;
coefs(nc+1) = 0.3;


src_info = geometries.starn(coefs,nc,n);
L = src_info.L;





% Test obstacle Frechet derivative for Dirichlet problem
bc = [];
bc.type = 'Dirichlet';
bc.invtype = 'o';


src0 = [0.01;-0.12];
opts = [];
opts.test_analytic = true;
opts.src_in = src0;
opts.verbose=true;


m  = ceil(60*1.5);
m = 100;
% set target locations
%receptors (r_{\ell})
r_tgt = 10;
n_tgt = m;
t_tgt = 0:2*pi/n_tgt:2*pi-2*pi/n_tgt;

% Incident directions (d_{j})
n_dir = 100;
%n_dir = m;
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



kh = 1;
nppw = 20;


n = ceil(nppw*L*abs(kh)/2/pi);
n = max(n,300);
src_info = geometries.starn(coefs,nc,n);

[mats,erra] = rla.get_fw_mats(kh,src_info,bc,sensor_info,opts);
fields = rla.compute_fields(kh,src_info,mats,sensor_info,bc,opts);

u_meas = [];
u_meas.kh = kh;
u_meas.uscat_tgt = fields.uscat_tgt;
u_meas.tgt = sensor_info.tgt;
u_meas.t_dir = sensor_info.t_dir;
u_meas.err_est = erra;



alpha = 1e-3;

[Ig,xgrid0,ygrid0] = lsm.lsm_tensor(n_tgt,n_dir,u_meas,alpha);
figure; surf(xgrid0,ygrid0,Ig); shading interp; view(2); hold on;  
plot3(src_info.xs,src_info.ys,100*ones(size(src_info.xs)),'k')

figure; contour(xgrid0,ygrid0,Ig); hold on; plot(src_info.xs,src_info.ys,'k')
