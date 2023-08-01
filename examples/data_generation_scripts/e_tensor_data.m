% This script generates the data for a letter e like domain, for a fixed
% number of sensors and incident directions where data is available for all
% sensors at each incident direction



addpath('../');


% choose width for e domain
% iwid = 1, all legs equal
% iwid = 2, middle leg half of other two
% iwid = 3, middle leg quarater of other two
iwid = 3;


srcin = load('+geometries/letter_e1.dat');
srcin = srcin.';
L0 = load('+geometries/letter_e1_length.dat');


n = 400;
src_info = geometries.get_e(n,iwid);

L = src_info.L;
plot(src_info.xs,src_info.ys,'k.');
axis equal





nk = 27;


% Test obstacle Frechet derivative for Dirichlet problem
bc = [];
bc.type = 'Dirichlet';
bc.invtype = 'o';

fname = ['../data/letter_e' int2str(iwid) '_10_ik1_nk' int2str(nk) '_tensor_data_' bc.type 'carlos_test.mat'];
save(fname,'src_info');


src0 = [0.6;0.12];
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



src_info = geometries.get_h(n);

% Set of frequencies (k_{i})
dk = 0.25;

% kh = 7.5*1.5;
% kh = 10;
kh = 1:dk:(1+(nk-1)*dk);



u_meas = cell(nk,1);



nppw = 20;

for ik=1:nk
   n = ceil(nppw*L*abs(kh(ik))/2/pi);
   n = max(n,384);
   src_info = geometries.get_e(n,iwid);
   
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


