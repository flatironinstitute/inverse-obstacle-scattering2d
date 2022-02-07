fname = '../data/star3_ik1_nk30_tensor_data_Dirichlet.mat';
load(fname);
addpath('../');
addpath('../../src');

optim_opts = [];
opts = [];
opts.verbose=true;
bc = [];
bc.type = 'Dirichlet';
bc.invtype = 'o';
optim_opts.optim_type = 'gn';
%optim_opts.eps_curv = 1e-3;
optim_opts.filter_type = 'gauss-conv';
%optim_opts.eps_res = 1e-10;
%optim_opts.eps_upd = 1e-10;
opts.store_src_info = true;
[inv_data_all,src_info_out] = rla.rla_inverse_solver(u_meas,bc,...
                          optim_opts,opts);
rla.post_process(inv_data_all,fname);