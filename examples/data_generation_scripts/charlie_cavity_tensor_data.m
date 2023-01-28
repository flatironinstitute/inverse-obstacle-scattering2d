% This script generates the data for a star shaped domain, for a fixed
% number of sensors and incident directions where data is available for all
% sensors at each incident direction


n  = 500;
addpath('../');



juse = 11;
src_info = geometries.charlie_cavity(juse,n);
L = src_info.L;




nk = 27;


% Test obstacle Frechet derivative for Dirichlet problem
bc = [];
bc.type = 'Dirichlet';
bc.invtype = 'o';


fname = ['../data/ccavity10_ik1_nk' int2str(nk) '_tensor_data_' bc.type 'carlos_test.mat'];
save(fname,'src_info');


src0 = [0.6;-0.012];
opts = [];
opts.test_analytic = true;
opts.src_in = src0;
opts.verbose=true;


m  = ceil(60*1.5);
m = 60;
% set target locations
%receptors (r_{\ell})
r_tgt = 10;
n_tgt = 100;
t_tgt = 0:2*pi/n_tgt:2*pi-2*pi/n_tgt;

% Incident directions (d_{j})
n_dir = 120;
n_dir = 16;
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


n = 500;
juse = 10;
src_info = geometries.charlie_cavity(juse,n);

kh = 10;

[mats,erra] = rla.get_fw_mats(kh,src_info,bc,sensor_info,opts);
plot(src_info.xs,src_info.ys,'k.');
fprintf('error in soln=%d\n',erra);



% Set of frequencies (k_{i})
dk = 0.25;

% kh = 7.5*1.5;
% kh = 10;
kh = 1:dk:(1+(nk-1)*dk);



u_meas = cell(nk,1);



nppw = 20;

for ik=1:nk
   n = ceil(nppw*L*abs(kh(ik))/2/pi);
   n = max(n,1000);
   src_info = geometries.charlie_cavity(juse,n);
   
   [mats,erra] = rla.get_fw_mats(kh(ik),src_info,bc,sensor_info,opts);
   fields = rla.compute_fields(kh(ik),src_info,mats,sensor_info,bc,opts);
   
   u_meas0 = [];
   u_meas0.kh = kh(ik);
   u_meas0.uscat_tgt = fields.uscat_tgt;
   u_meas0.tgt = sensor_info.tgt;
   u_meas0.t_dir = sensor_info.t_dir;
   u_meas0.err_est = erra;
   u_meas{ik} = u_meas0;
   fprintf('error in soln=%d\n',erra);
   
end


save(fname,'u_meas','-append');

figure
clf
uscat_tgt = reshape(fields.uscat_tgt,[n_dir,n_tgt]);
uhat = fft2(uscat_tgt);
d = abs(fftshift(uhat));
imagesc(abs(uscat_tgt));

colorbar();


figure
clf
imagesc(d)
colorbar();

figure;
clf();
plot(src_info.xs,src_info.ys,'b.');
axis equal;


