addpath('../src');
nterms = 70;
n = 800;
src_info = geometries.smooth_plane(nterms,n);
figure(1)
clf
plot(src_info.xs,src_info.ys,'b.'); axis equal


src0 = [0.01;-0.12];
opts = [];
opts.test_analytic = true;
opts.src_in = src0;

% set target locations
%receptors
r_tgt = 10;
n_tgt = 100;
t_tgt = 0:2*pi/n_tgt:2*pi-2*pi/n_tgt;
x_t   = r_tgt * cos(t_tgt);
y_t   = r_tgt * sin(t_tgt);    
tgt   = [ x_t; y_t];


sensor_info = [];
sensor_info.tgt = tgt;


kh = 2.2;

bc = [];
bc.type = 'Dirichlet';


[mats,erra] = rla.get_fw_mats(kh,src_info,bc,sensor_info,opts);

if(erra>1e-11) 
    ipass = 0;
    fprintf('failed Dirichlet test in fw_mats test: %d\n',erra);
end
