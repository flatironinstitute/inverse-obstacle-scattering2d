function [inverse_sol_data,src_info_out] = inverse_solver(kh,src_info,bc, ...
     u_meas,optim_opts,opts)
%
%  This subroutine endeavors to find the optimal \Gamma, \lambda such that
%  for a given frequency k the objective function
%
%  \sum_{j}| F(k,\bx_{j},d_{j}, \Gamma, \lambda - u_meas.uscat_tgt(j)|^2
%     is minimized
%
%  Here F(k,\bx_{j}, d_{j}, \Gamma, \lambda) denotes the scattered
%  field at frequency k, at target location \bx_{j}, due to incident 
%  field given by plane wave with direction $d_{j}$, where
%  the boundary of the curve is given by $\Gamma$ and the impedance
%  function given by $\lambda$ (note that the impedance function)
%  is optional, and u_meas denotes the corresponding measured data.
%
%  We support 5 optimization methods for updating the boundary/impedance
%   sd - Steepest descent with step size set to the Cauchy point
%   gn - Gauss Newton
%   min(gn,sd) or min(sd,gn) - Use the better of steepest descent or Gauss-newton
%   sd-gn - steepest descent for the first few iterations followed by
%     Gauss-newton
%   sd-min(gn,sd) or sd-min(sd,gn) - Steepest descent for first few
%   iterations followed by best of Gauss-newton and steepest descent.
%
%   If during any optimization step,
%   the updated curve is self-intersecting or if the curvature
%   is requested to have a nearly band limited tail then,
%   the above updates are modified by either scaling them by 2^i 
%   for i<=maxit_filtering, so that the conditions are met, or
%   by damping the high frequency components of the update using 
%   a Gaussian until the conditions are met.
%
%   If the conditions aren't met with the maximum filtering iterations, 
%   then the input curve is returned with an error code
%   
% Input:
%   kh - Helmholtz wave number
%   src_info - source info struct;
%      src_info.xs = x coordinates;
%      src_info.ys = y coordinates;
%      src_info.dxs = dxdt;
%      src_info.dys = dydt;
%      src_info.ds = sqrt(dxdt^2 + dydt^2);
%      src_info.h = h in trapezoidal parametrization;
%      src_info.lambda - imepdance value at discretization nodes 
%           (optional if solving impedance boundary value problem);
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
%
% Optional input arguments
%    optim_opts - optimization options struct, default values in brackets
%       optim_opts.eps_curv: constraint for high frequency content of
%             curvature of geometry (Inf)
%       optim_opts.ncurv: fourier coefficient number for defining 
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
%      opts.ncoeff_impedance_max: maximum number of coeffs to be used
%         in representing impedance. (Inf)
%      opts.nppw: points per wavelength for discretizing updated curve (10)
%      opts.verbose: flag for displaying verbose messages during run
%      (false)
%      opts.store_fields: store fields corresponding to each iterate
%      (false)
%      opts.store_src_info: store src_info corresponding to each iterate
%      (false)
%
% Output
% inverse_sol_data: struct containing info through the optimization process
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
% src_info_out: updated source struct upon exiting optimization loop
 
    inverse_sol_data = [];
    
    if(nargin<5)
      opts = [];
    end

    if(nargin<4)
      optim_opts = [];
    end
    
    verbose = false;
    if(isfield(opts,'verbose'))
        verbose = opts.verbose;
    end
    
    
    opts_use = opts;
    if(isfield(opts_use,'ncoeff_boundary_mult'))
        opts_use.ncoeff_boundary = floor(kh*opts_use.ncoeff_boundary_mult);
    end
    
    if(isfield(opts_use,'ncoeff_impedance_mult'))
        opts_use.ncoeff_impedance = floor(kh*opts_use.ncoeff_impedance_mult);
        if (isfield(opts_use,'ncoeff_impedance_max'))
            opts_use.ncoeff_impedance = ...
                min(opts_use.ncoeff_impedance, ...
                opts_use.ncoeff_impedance_max);
        end
    end
    
    maxit = 100;
    if(isfield(optim_opts,'maxit'))
        maxit = optim_opts.maxit;
    end
    
   
    optim_type = 'gn';
    if(isfield(optim_opts,'optim_type'))
        optim_type = optim_opts.optim_type;
    end
    
    sd_iter = 0;
    
    if(strcmpi(optim_type,'sd-gn') || strcmpi(optim_type,'sd-min(sd,gn)') || strcmpi(optim_type,'sd-min(gn,sd)'))
        sd_iter = 50;
        if(isfield(optim_opts,'sd_iter'))
           sd_iter = optim_opts.sd_iter;
        end
    end
    
    eps_res = 1e-5;
    if(isfield(optim_opts,'eps_res'))
        eps_res = optim_opts.eps_res;
    end
    
    eps_upd = 1e-5;
    if(isfield(optim_opts,'eps_upd'))
        eps_upd = optim_opts.eps_upd;
    end
    
    src_info_all = cell(maxit,1);
    res_all = zeros(maxit,1);
    exit_criterion = 0;
    fields_all = cell(maxit,1);
    deltas = cell(maxit,1);
    ier = zeros(maxit,1);
    store_fields = false;
    if(isfield(opts,'store_fields'))
        store_fields = opts.store_fields;
    end
    
    store_src_info = false;
    if(isfield(opts,'store_src_info'))
        store_src_info = opts.store_src_info;
    end
    
%
%  Initialize matrices and fields for iteration 0
%
    src_use = src_info;
    mats = rla.get_fw_mats(kh,src_use,bc,u_meas,opts);
    fields = rla.compute_fields(kh,src_use,mats,u_meas,bc,opts);
    
    if(verbose)
      fprintf('-------------------------------------\n')
      fprintf('Starting inverse interation for frequency kh: %d\n',kh);
    end
    
    optim_opts_use0 = optim_opts;
    optim_opts_use0.optim_type = optim_type;
    optim_opts_use0.eps_res = eps_res;
    optim_opts_use0.eps_upd = eps_upd;
    optim_opts_use0.sd_iter = sd_iter;
    optim_opts_use0.maxit = maxit;
    ict = 1;
    for iter=1:maxit
        optim_opts_use = optim_opts_use0;
        if(strcmpi(optim_type,'sd-gn') || strcmpi(optim_type,'sd-min(sd,gn)') || strcmpi(optim_type,'sd-min(gn,sd)'))
            if(iter<=sd_iter)
                optim_opts_use.optim_type = 'sd';
            else
                optim_opts_use.optim_type=optim_type(4:end);
            end
        end
        
        [deltas{iter},src_info_all{iter},mats_out,fields_all{iter},res_all(iter),ier(iter)] = ...
           rla.update_inverse_iterate(kh,src_use,mats,fields,u_meas,bc,optim_opts,opts_use);
        if(verbose)
            fprintf('iter number: %d \t optim_type: %s \t residue: %d \t ier: %d\n',iter,optim_opts_use.optim_type,res_all(iter),ier(iter));
        end
        
        if(ier(iter) ~=0) 
            exit_criterion = -1;
            break;
        end
        
        
        src_use = src_info_all{iter};
        mats = mats_out;
        fields = fields_all{iter};
        
        if(res_all(iter) <= eps_res)
            exit_criterion = 1;
            break;
        end
        
        res_upd_norm = 0;
        if(isfield(deltas{iter},'delta_bdry'))
            res_upd_norm = res_upd_norm + norm(deltas{iter}.delta_bdry(:));
        end
        
        if(isfield(deltas{iter},'delta_impedance'))
            res_upd_norm = res_upd_norm + norm(deltas{iter}.delta_impedance(:));
        end
        
        if(res_upd_norm <= eps_upd)
            exit_criterion = 2;
            break;
        end
        ict = ict + 1;
    end
    if(ict == (maxit +1)) 
        exit_criterion = 3;
    end
    if(verbose)
        fprintf('Exit criterion: %d \t iteration count: %d\n',exit_criterion,iter);
    end
    src_info_out = src_use;
    inverse_sol_data.res_all = res_all(1:iter);
    inverse_sol_data.iter_count = iter;
    inverse_sol_data.res_opt = res_all(iter);
    if(store_src_info)
        inverse_sol_data.src_info_all = src_info_all(1:iter);
    end
    inverse_sol_data.src_info_opt = src_use;
    if(store_fields)
        inverse_sol_data.fields_all = fields_all(1:iter); 
    end
    inverse_sol_data.deltas = deltas(1:iter);
    inverse_sol_data.fields_opt = fields_all{iter};
    inverse_sol_data.exit_criterion = exit_criterion;
    inverse_sol_data.optim_opts = optim_opts;
    inverse_sol_data.kh = kh;
    if(verbose)
      fprintf('Completing inverse interation for frequency kh: %d\n',kh);
      fprintf('-------------------------------------\n');
    end
    
end
