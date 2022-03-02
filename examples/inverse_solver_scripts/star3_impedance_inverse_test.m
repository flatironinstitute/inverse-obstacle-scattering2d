fname = '../data/star3_ik1_nk17_tensor_data_Impedance.mat';
load(fname);
addpath('../');
addpath('../../src');
maxNumCompThreads(4)
warning('off')

optim_opts = [];
opts = [];
opts.verbose=true;
opts.ncoeff_impedance_mult = 2;
bc = [];
bc.type = 'Impedance';
bc.invtype = 'i';
optim_opts.optim_type = 'gn';
optim_opts.eps_curv = 1e-3;
optim_opts.filter_type = 'gauss-conv';
optim_opts.eps_res = 1e-7;
optim_opts.eps_upd = 1e-7;
optim_opts.maxit = 200;
opts.store_src_info = true;

if(strcmpi(bc.invtype,'i'))
    A = load(fname);
    nh = 1;
    hcoefs = zeros(2*nh+1,1);
    [src_init,varargout] = rla.update_geom(A.src_info,nh,hcoefs);
    lam = [];
elseif(strcmpi(bc.invtype,'o'))
    src_init = [];
else
    src_init = [];
    lam = [];
end
    


[inv_data_all,src_info_out] = rla.rla_inverse_solver(u_meas,bc,...
                          optim_opts,opts,src_init,lam);
rla.post_process(inv_data_all,fname);

%save('../data/sol_star3_ik1_nk57_tensor_data_imepdance8.mat','inv_data_all','fname','-v7.3');
