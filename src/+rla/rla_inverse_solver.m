function [inv_data_all,src_info_out] = rla_inverse_solver(u_meas,bc,...
                          optim_opts,opts,src_init,lam_init,rla_path_opts)
% 
%  This subroutine endeavors to find the optimal \Gamma, \lambda such that
%  for a given path in frequency space the objective function
%
%  \sum_{j}| F(k_{ell},\bx_{j},d_{j}, \Gamma, 
%     \lambda - u_meas(ik_list(\ell)).uscat_tgt(j)|^2
%     is minimized for each \ell,
%
%  where the path in k_{\ell}, ik_list(\ell) can be specified by the user
%  and if unspecified, the path taken ik_list(\ell) = 1:nfreq
%  where nfreq is the length of the u_meas struct.
%
%
%  Here F(k,\bx_{j}, d_{j}, \Gamma, \lambda) denotes the scattered
%  field at frequency k, at target location \bx_{j}, due to incident 
%  field given by plane wave with direction $d_{j}$, where
%  the boundary of the curve is given by $\Gamma$ and the impedance
%  function given by $\lambda$ (note that the impedance function)
%  is optional, and u_meas denotes the corresponding measured data.
%
% For each individual optimization problem, we support 5 optimization
% methods
%   sd - Steepest descent with step size set to the Cauchy point
%   gn - Gauss Newton
%   min(gn,sd) or min(sd,gn) - Use the better of steepest descent or Gauss-newton
%   sd-gn - steepest descent for the first few iterations followed by
%     Gauss-newton
%   sd-min(gn,sd) or sd-min(sd,gn) - Steepest descent for first few
%   iterations followed by best of Gauss-newton and steepest descent.
%
% Input arguments
%  u_meas - measurement data struct
%      u_meas.tgt(2,nmeas) - xy cooordinates of sensors
%         u_meas.tgt(1:2,i) = xy coordinates corresponding the ith
%            measurement
%      u_meas.t_dir(nmeas) - incident directions
%         u_meas.t_dir(i) - is the incident direction corresponding to
%            the ith measurement
%      u_meas.uscat_tgt(nmeas) - scattered field corresponding to the 
%          measurements
%   bc - boundary condition struct;
%     bc.type = type of boundary condition;
%        'd' or 'Dirichlet' for dirichlet
%        'n' or 'Neumann' for Neumann
%        'i' or 'Impedance' for impedance
%     bc.invtype = type of inverse problem;
%        'o' or 'obstacle' for obstacle only
%        'oi' or 'obsctacle and impedance' for both obstacle and impedance;
% 
% Optional input arguments
%    optim_opts - optimization options struct, default values in brackets
%       optim_opts.eps_curv: constraint for high frequency content of
%             curvature of geometry (Inf)
%       optim_opts.n_curv: fourier coefficient number for defining 
%           the tail of curvature for imposing above constraint
%           (max(30,opts.ncoeff_boundary))
%       optim_opts.optim_type: optimization type ('gn')
%            'gn' - gauss newton
%            'sd' - steepest descent
%            'min(sd,gn)','min(gn,sd)' - best of gauss newton or steepest
%            descent
%            'sd-gn' - steepest descent followed b gauss-newton
%            'sd-min(sd,gn)' or 'sd-min(gn,sd)' - steepest descent
%               followed by best of steepest descent or gauss-newton
%       optim_opts.filter_type: curve update filteration type, invoked if
%          curve is self intersecting or curvature does not satisfy
%          contraint ('gauss-conv')
%             'gauss-conv' - convolve update with gaussian
%             'step_length' - decrease step size by factor of 2
%       optim_opts.maxit_filter: maximum number of iterations of applying
%          the curve update filter (10)
%       optim_opts.maxit: maximum number of outer optimization iterations
%       (100)
%       optim_opts.sd_iter: number of initial steepest descent iterations
%           if applicable (50)
%       optim_opts.eps_upd: exit criterion for declaring local optimum
%          based on l2 norm of update (1e-5)
%       optim_opts.eps_res: exit criterion for declaring local optimum
%          based on l2 norm of relative residual (1e-5)
%
%   opts - options struct, default values in brackets
%      opts.ncoeff_boundary_mult: multiplier for 
%         number of terms in boundary update given by
%         ceil(ncoeff_boundary_mult*abs(kh)) (2)
%      opts.ncoeff_impedance_mult: multiplier for number of terms in impedance update
%         ceil(ncoeff_impedance_mult*abs(kh)) (0.5)
%      opts.ncoeff_impedance_max: maximum number of terms in impedance 
%         update (Inf)
%      opts.nppw: points per wavelength for discretizing updated curve (10)
%      opts.verbose: flag for displaying verbose messages during run
%      (false)
%      opts.store_fields: store fields corresponding to each iterate
%      (false)
%      opts.store_src_info: store src_info corresponding to each iterate
%      (false)
%
%  src_init: initial shape to start the first optimization problem
%     if unspecified, a unit circle centered at the origin is used
%
%  rla_path_opts.ik_list: path to be taken in frequency space specified through
%    a list of indices between 1:nfreq, can be of longer length (1:1:nfreq)
%
% Output arguments:
%  inv_data_all (npath): inverse_sol_data struct for each optimization
%  problem along the path, here npath is the length of
%  rla_path_opts.ik_list 
%   Each cell element is itself a struct of type inverse_sol_data 
%    which contains the following entries:
%    inverse_sol_data.res_all(1:iter_count): residue at each iterate
%    inverse_sol_data.iter_count: number of iterations before exiting
%       optimization
%    inverse_sol_data.src_info_all: source info struct at each iter (if
%      requested)
%    inverse_sol_data.src_info_opt: same as src_info_out, optimum src_info
%      source struct
%    inverse_sol_data.fields_all: fields at each iterate (if requested)
%    inverse_sol_data.deltas: deltas struct (update info) for each iterate
%    inverse_sol_data.fields_opt: fields corresponding to src_info_opt
%    inverse_sol_data.optim_opts: final optimization parameters used in the
%      run
%    inverse_sol_data.exit_criterion: exit criterion for optimization
%    problem
%      exit_criterion = -1; at the last iterate curve was self intersecting
%        or did not satisfy curvature constraint if requested
%      exit_criterion = 1; if res<eps_res
%      exit_criterion = 2; if |update| < eps_upd
%      exit_criterion = 3; out of iterations, i.e. iter_count = maxit;
%
%  src_info_out: optimal shape at the end of optimization path

                      
                      
   if(nargin<7)
       rla_path_opts = [];
   end
   if(nargin<6) 
       lam_init = [];
   end
   if(nargin<5)
       src_init = [];
   end
   if(nargin<4)
       opts = [];
   end
   if(nargin<3)
       optim_opts = [];
   end
   
   if(isfield(rla_path_opts,'ik_list'))
       ik_list = rla_path_opts.ik_list;
   else
       ik_list = 1:1:length(u_meas);
   end
   npath = length(ik_list);
   
   nppw = 10;
   if(isfield(opts,'nppw'))
       nppw = opts.nppw;
   end
   
   if(isempty(src_init))
       kh0 = u_meas{ik_list(1)}.kh;
       n = max(100,ceil(abs(kh0)*nppw));
       src_init = geometries.ellipse(1,1,n);     
   end
   
   
   inv_data_all = cell(1,npath);
   src_info = src_init;
   n = length(src_info.xs);
   if(isempty(lam_init))
       src_info.lambda = ones(n,1);
   else
       % fix this initialization
       t = 0:2*pi/n:2*pi*(1-1/n);
       src_info.lambda = lam_init(t(:));
       
   end
   
   for i=1:npath
       kh = u_meas{ik_list(i)}.kh;
       if(isfield(opts,'use_lscaled_modes'))
          if(opts.use_lscaled_modes)
             opts.ncoeff_boundary_mult = 2*src_info.L/2/pi;
          end
       end
       [inv_data_all{i},src_out] = rla.inverse_solver(kh,src_info,bc, ...
          u_meas{ik_list(i)},optim_opts,opts);
       src_info = src_out;
   end
   src_info_out = src_info;
end
