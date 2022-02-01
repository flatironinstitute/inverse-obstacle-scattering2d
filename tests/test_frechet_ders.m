
n  = 300;
addpath('../');
a = 1.1;  
b = 1.3;  
src_info = geometries.ellipse(a,b,n);

plot(src_info.xs,src_info.ys,'b.');
drawnow;

src0 = [0.01;-0.12];
opts = [];
opts.test_analytic = true;
opts.src_in = src0;

% set target locations
%receptors
r_tgt = 10;
n_tgt = 88;
t_tgt = 0:2*pi/n_tgt:2*pi-2*pi/n_tgt;

% Incident directions
n_dir = 100;
t_dir = 0:2*pi/n_dir:2*pi-2*pi/n_dir;

[t_tgt_grid,t_dir_grid] = meshgrid(t_tgt,t_dir);
t_tgt_grid = t_tgt_grid(:);
t_dir_grid = t_dir_grid(:);
xtgt = r_tgt*cos(t_tgt_grid);
ytgt = r_tgt*sin(t_tgt_grid);
tgt   = [ xtgt'; ytgt'];

% use a random subset of sensor and target location
% so that they are no longer tensor product
ifac = 1;
iind_use = randperm(n_dir*n_tgt,ceil(ifac*n_dir*n_tgt));
iind_use = 1:(n_dir*n_tgt);
sensor_info = [];
sensor_info.tgt = tgt(:,iind_use);
sensor_info.t_dir = t_dir_grid(iind_use);


[t_dir_uni,~,idir] = unique(sensor_info.t_dir);
[tgt_uni,~,itgt] = unique(sensor_info.tgt','rows');
nt_uni = length(tgt_uni(:,1));
induse = itgt + (idir-1)*nt_uni;


kh = 2.2;

% Test obstacle Frechet derivative for Dirichlet problem
bc = [];
bc.type = 'Dirichlet';
bc.invtype = 'o';

[mats,~] = rla.get_fw_mats(kh,src_info,bc,sensor_info,opts);
fields = rla.compute_fields(kh,src_info,mats,sensor_info,bc,opts);

nh = 4;
opts.ncoeff_boundary = nh;
opts.ncoeff_impedance = nh;
frechet_mats = rla.get_frechet_ders(kh,mats,src_info,sensor_info,fields,bc,opts);

rng(1)
hcoefs = 0.1*rand(1,2*nh+1);
uder = frechet_mats.bdry*hcoefs(:);


errs = zeros(2,1);
for ig=1:2
    dh = 10^(-ig);
    hcoefs_use = dh*hcoefs;
    [src_out,~] = rla.update_geom(src_info,nh,hcoefs_use);
    [mats1,~] = rla.get_fw_mats(kh,src_out,bc,sensor_info,opts);
    fields1 = rla.compute_fields(kh,src_out,mats1,sensor_info,bc,opts);

    
    hcoefs_use2 = -dh*hcoefs;
    [src_out2,~] = rla.update_geom(src_info,nh,hcoefs_use2);
    [mats2,~] = rla.get_fw_mats(kh,src_out2,bc,sensor_info,opts);
    fields2 = rla.compute_fields(kh,src_out2,mats2,sensor_info,bc,opts);
    
    uder_est = (fields1.uscat_tgt(:) - fields2.uscat_tgt(:))/2/dh;
    errs(ig) = norm(uder-uder_est);
end
fprintf('Error in Frechet derivative of Dirichlet problem: %d \n',errs(2));
errs

% Test obstacle Frechet derivative for Neumann problem
bc = [];
bc.type = 'Neumann';
bc.invtype = 'o';

[mats,~] = rla.get_fw_mats(kh,src_info,bc,sensor_info,opts);
fields = rla.compute_fields(kh,src_info,mats,sensor_info,bc,opts);


frechet_mats = rla.get_frechet_ders(kh,mats,src_info,sensor_info,fields,bc,opts);
uder = frechet_mats.bdry*hcoefs(:);


errs = zeros(2,1);
for ig=1:2
    dh = 10^(-ig);
    hcoefs_use = dh*hcoefs;
    [src_out,~] = rla.update_geom(src_info,nh,hcoefs_use);
    [mats1,~] = rla.get_fw_mats(kh,src_out,bc,sensor_info,opts);
    fields1 = rla.compute_fields(kh,src_out,mats1,sensor_info,bc,opts);

    
    hcoefs_use2 = -dh*hcoefs;
    [src_out2,~] = rla.update_geom(src_info,nh,hcoefs_use2);
    [mats2,~] = rla.get_fw_mats(kh,src_out2,bc,sensor_info,opts);
    fields2 = rla.compute_fields(kh,src_out2,mats2,sensor_info,bc,opts);
    
    uder_est = (fields1.uscat_tgt(:) - fields2.uscat_tgt(:))/2/dh;
    errs(ig) = norm(uder-uder_est);
end
fprintf('Error in Frechet derivative of Neumann problem: %d \n',errs(2));


% Test obstacle frechet derivative of impedance problem
t = 0:2*pi/n:2*pi*(1.0-1.0/n);
src_info.lambda = sin(2*t)';
bc = [];
bc.type = 'Impedance';
bc.invtype = 'oi';

[mats,~] = rla.get_fw_mats(kh,src_info,bc,sensor_info,opts);
fields = rla.compute_fields(kh,src_info,mats,sensor_info,bc,opts);


frechet_mats = rla.get_frechet_ders(kh,mats,src_info,sensor_info,fields,bc,opts);
uder = frechet_mats.bdry*hcoefs(:);


errs = zeros(2,1);
for ig=1:2
    dh = 10^(-ig);
    hcoefs_use = dh*hcoefs;
    [src_out,~] = rla.update_geom(src_info,nh,hcoefs_use);
    [mats1,~] = rla.get_fw_mats(kh,src_out,bc,sensor_info,opts);
    fields1 = rla.compute_fields(kh,src_out,mats1,sensor_info,bc,opts);

    
    hcoefs_use2 = -dh*hcoefs;
    [src_out2,~] = rla.update_geom(src_info,nh,hcoefs_use2);
    [mats2,~] = rla.get_fw_mats(kh,src_out2,bc,sensor_info,opts);
    fields2 = rla.compute_fields(kh,src_out2,mats2,sensor_info,bc,opts);
    
    uder_est = (fields1.uscat_tgt(:) - fields2.uscat_tgt(:))/2/dh;
    errs(ig) = norm(uder-uder_est);
end
fprintf('Error in obsctacle Frechet derivative of Impedance problem: %d \n',errs(2));



% Test impedance frechet derivative of impedance problem

uder = frechet_mats.impedance*hcoefs(:);
errs = zeros(2,1);
for ig=1:2
    dh = 10^(-ig);
    hcoefs_use = dh*hcoefs;
    hcoefs_use = hcoefs_use(:);
    src_out = src_info;
    h_upd = (cos(t'*(0:nh))*hcoefs_use(1:(nh+1)) + sin(t'*(1:nh))*hcoefs_use((nh+2):end)).';
    src_out.lambda = src_info.lambda + h_upd';
    [mats1,~] = rla.get_fw_mats(kh,src_out,bc,sensor_info,opts);
    fields1 = rla.compute_fields(kh,src_out,mats1,sensor_info,bc,opts);

    
    
    src_out2 = src_info;
    src_out2.lambda = src_info.lambda - h_upd';
    [mats2,~] = rla.get_fw_mats(kh,src_out2,bc,sensor_info,opts);
    fields2 = rla.compute_fields(kh,src_out2,mats2,sensor_info,bc,opts);
    
    uder_est = (fields1.uscat_tgt(:) - fields2.uscat_tgt(:))/2/dh;
    errs(ig) = norm(uder-uder_est);
end


fprintf('Error in impedance Frechet derivative of Impedance problem: %d \n',errs(2));

