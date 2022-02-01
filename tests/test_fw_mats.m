
n  = 300;
addpath('../');
a = 1.1;  % thickness: cannot be >1/3 otherwise not smooth to emach
b = 1.3;  % controls approx opening angle in radians (keep small for resonant)
coefs = [1.0 0 0 0.4 0 0 0];
nc = 3;
src_info = geometries.starn(coefs,nc,n);

plot(src_info.xs,src_info.ys,'b.');
drawnow;

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
fprintf('Error in dirichlet problem: %d\n',erra);


bc = [];
bc.type = 'Neumann';
[mats,erra] = rla.get_fw_mats(kh,src_info,bc,sensor_info,opts);
fprintf('Error in Neumann problem: %d\n',erra);


src_info.lambda = ones(n,1);
bc = [];
bc.type = 'Impedance';
[mats,erra] = rla.get_fw_mats(kh,src_info,bc,sensor_info,opts);
fprintf('Error in Impedance problem: %d\n',erra);
