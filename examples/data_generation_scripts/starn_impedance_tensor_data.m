% This script generates the data for a star shaped domain, for a fixed
% number of sensors and incident directions where data is available for all
% sensors at each incident direction


n  = 300;
addpath('../');
addpath('../../');


% max number of wiggles
nc = 3;

% parameters 'a'
coefs = zeros(2*nc+1,1);
coefs(1) = 1;
coefs(nc+1) = 0.0;


src_info = geometries.starn(coefs,nc,n);
L = src_info.L;

t = 0:2*pi/n:2*pi*(1.0-1.0/n);
lam = @(t) 2 + sin(2*t) + 0.3*cos(6*t) + 0.1*sin(10*t);
%lam = @(t) 1.4;
%lam = @(t) 1;
src_info.lambda = lam(t)';
nk = 17;


% Test obstacle Frechet derivative for Dirichlet problem
bc = [];
bc.type = 'Impedance';
bc.invtype = 'o';

fname = ['../data/star0_ik1_nk' int2str(nk) '_tensor_data_' bc.type '.mat'];
save(fname,'src_info','lam');


src0 = [0.01;-0.12];
opts = [];
opts.test_analytic = true;
opts.src_in = src0;
opts.verbose=true;


m  = ceil(60*1.5);
% set target locations
%receptors (r_{\ell})
r_tgt = 10;
n_tgt = m;
t_tgt = 0:2*pi/n_tgt:2*pi-2*pi/n_tgt;

% Incident directions (d_{j})
n_dir = m;
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



src_info = geometries.starn(coefs,nc,n);

% Set of frequencies (k_{i})
dk = 0.25;
kh = 1:dk:(1+(nk-1)*dk);
%kh = 7.5*1.5;

u_meas = cell(nk,1);



nppw = 20;

for ik=1:nk
   n = ceil(nppw*L*abs(kh(ik))/2/pi);
   n = max(n,300);
   src_info = geometries.starn(coefs,nc,n);
   nh = 1;
   hcoefs = zeros(2*nh+1,1);
   [src_info,varargout] = rla.update_geom(src_info,nh,hcoefs);
   
   
   t = 0:2*pi/n:2*pi*(1.0-1.0/n);
   src_info.lambda = lam(t)';
   
   
   [mats,erra] = rla.get_fw_mats(kh(ik),src_info,bc,sensor_info,opts);
   fields = rla.compute_fields(kh(ik),src_info,mats,sensor_info,bc,opts);
   
   u_meas0 = [];
   u_meas0.kh = kh(ik);
   u_meas0.uscat_tgt = fields.uscat_tgt;
   u_meas0.tgt = sensor_info.tgt;
   u_meas0.t_dir = sensor_info.t_dir;
   u_meas0.err_est = erra;
   u_meas{ik} = u_meas0;
   disp(erra)
end


save(fname,'u_meas','-append');
uscat_tgt = reshape(fields.uscat_tgt,[m,m]);
uhat = fft2(uscat_tgt);
d = log10(abs(fftshift(uhat)));
imagesc(d);
colorbar();


figure;
clf();
plot(src_info.xs,src_info.ys,'b.');
axis equal;


