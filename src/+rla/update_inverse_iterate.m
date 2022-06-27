function [deltas,src_out,mats_out,fields_out,res,ier] = ...
   update_inverse_iterate(kh,src_info,mats,fields,u_meas,bc,optim_opts,opts)
%
%  This subroutine updates the boundary and/or the impedance function
%  of the curve, based on one optimization step of 
%  minimizing the objective function of the scattered field
%  at a collection of target locations, and incident directions
%  for a single frequency kh.
%
%  Let F(k,\bx_{j}, d_{j}, \Gamma, \lambda) denote the scattered
%  field at frequency k, at target location \bx_{j}, due to incident 
%  field given by plane wave with direction $d_{j}$, where
%  the boundary of the curve is given by $\Gamma$ and the impedance
%  function given by $\lambda$ (note that the impedance function)
%  is optional, and let u_meas denote the corresponding measured data.
%
%  We support 3 optimization methods for updating the boundary/impedance
%   sd - Steepest descent with step size set to the Cauchy point
%   gn - Gauss Newton
%   min(gn,sd) or min(sd,gn) - Use the better of steepest descent or Gauss-newton
%
%  The code sequentially updates the boundary of the obstacle first if requested,
%  followed by the impedance if requested.
%
%  Suppose J_b is the Frechet derivative of F with respect to the boundary
%  where the update is assumed to be a scalar function \delta_b in the
%  normal direction, and the scalar function is discretized in a fourier
%  basis.
%  Then if optimization method = sd, then
%    \delta_{b} = t J_{b}^{T} (u_meas - F); 
%      with t = (u_meas - F)^{T} J_{b} J_{b}^{T} (u_meas - F)/|J_{b}
%      J_{b}^{T} (u_meas-F)|^2;
%  
%   If the optimization method = gn, then
%     \delta_{b} = J_{b}^{+} (u_meas -F); 
%
%   where J_{b}^{+} is the pseudo inverse of J_{b}
%
%   If the updated curve is self-intersecting or if the curvature
%   is requested to have a nearly band limited tail then,
%   the above updates are modified by either scaling them by 2^i 
%   for i<=maxit_filtering, so that the conditions are met, or
%   by damping the high frequency components of the update using 
%   a Gaussian until the conditions are met.
%
%   If the conditions aren't met with the maximum filtering iterations, 
%   then the input curve is returned with an error code
%
%   If the optimization scheme is min(sd,gn) or min(gn,sd), then
%   both the gn and the sd updates are computed, and the updated curve
%   which has a smaller residue |u_meas - F(.,.,\Gamma_new,.)|
%   is returned
%
%   For the impedance, a similar procedure is applied after updating the
%   boundary. However, there are no filtering involved in the case
%   of the impedance update.
%-------
%   NOTE: IMPEDANCE currently only runs in gauss-newton mode. This needs
%    to be fixed
%------
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
%   mats - matrix structure
%    mats.Fw_mat - Matrix corresponding to discretizing the boundary
%    integral equation
%    mats.inv_Fw_mat = inverse of mats.Fw_mat
%    mats.Fw_dir_mat = matrix to obtain dirichlet data from the given
%       representation
%    mats.Fw_neu_mat = matrix to obtain Neumann data from given
%       representation (note that this is slightly different when
%       solving Dirichlet problem, see documentation below)
%    mats.sol_to_receptor - matrix mapping solution on boundary to 
%       potential at receptors
%    mats.bdrydata_to_receptor - matrix mapping boundary data to 
%       potential at receptors
%  fields - fields struct
%    fields.uinc(n,ndir) - incident field on the boundary due to
%      the incident fields given by the unique (ndir) directions 
%      in the sensor_info.tdir array
%    fields.dudninc(n,ndir) - the corresponding normal derivative of 
%       the incident field on the boundary
%    fields.uscat(n,ndir) - scattered field on the boundary due to
%      the incident fields given by the unique (ndir) directions 
%      in the sensor_info.tdir array
%    fields.dudnscat(n,ndir) - the corresponding normal derivative of 
%       scattered field on the boundary
%    fields.uscat_tgt(nmeas) - scattered field at the sensor locations
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
%       optim_opts.filter_type: curve update filteration type, invoked if
%          curve is self intersecting or curvature does not satisfy
%          contraint ('gauss-conv')
%             'gauss-conv' - convolve update with gaussian
%             'step_length' - decrease step size by factor of 2
%       optim_opts.maxit_filter: maximum number of iterations of applying
%          the curve update filter (10)
%
%   opts - options struct, default values in brackets
%      opts.ncoeff_boundary: number of terms in boundary update
%         ceil(2*abs(kh))
%      opts.ncoeff_impedance: number of terms in impedance update
%         ceil(0.5*abs(kh))
%      opts.nppw: points per wavelength for discretizing updated curve (10)
%      opts.verbose: flag for displaying verbose messages during run
%      (false)
%      opts.constphasefactor : flag, allow for a constant phase difference 
%                        between measurements and model. for each update,
%                        the current optimal phase is used and the
%                        derivatives are updated accordingly (false)
%                        
%
% Output argument
%   deltas: update struct
%      deltas.nmodes_bdry: number of modes used to update the boundary
%      deltas.nmodes_imp: number of modes used to update the impedance
%      deltas.delta_bdry: boundary update coefficients
%      deltas.delta_imp: impedance update coefficients
%      deltas.iter_filter_bdry: number of iterations of filter used to
%         obtain the update
%      deltas.iter_type: type of optimization used in the update step
%      deltas.phase: optimal constant phase if used (if not 
%                           used, set to 1) 
%   src_out: ouput source info struct
%   mats_out: output mats struct
%   fields_out: output fields struct
%   res: relative residue
%   ier: error code
%      ier = 0 implies successful execution
%

    deltas = [];
    if(nargin<8)
        opts = [];
    end

    if(nargin<7)
        optim_opts = [];
    end
   
    opts_use = opts;
    if(~isfield(opts,'ncoeff_boundary'))
        ncoeff_boundary = floor(2*abs(kh));
        opts_use.ncoeff_boundary = floor(2*abs(kh));
    else
        ncoeff_boundary = opts.ncoeff_boundary;
    end
    if(~isfield(opts,'ncoeff_impedance'))
        opts_use.ncoeff_impedance = floor(0.5*abs(kh));
        ncoeff_impedance = floor(0.5*abs(kh));
    else
        ncoeff_impedance = opts.ncoeff_impedance;
    end
    
    lambdareal = true;
    if(isfield(opts,'lambdareal'))
        lambdareal = opts.lambdareal;
    end
    opts_use.lambdareal = lambdareal;
    
    constphasefactor = false;
    if(isfield(opts,'constphasefactor'))
        constphasefactor = opts.constphasefactor;
    end
    opts_use.constphasefactor = constphasefactor;
    
    nppw = 10;
    rlam = Inf;
    if(isfield(opts,'nppw'))
        nppw = opts.nppw;
        rlam = 2*pi/real(kh);
    end
    
    eps_curv = Inf;
    if(isfield(optim_opts,'eps_curv'))
        eps_curv = optim_opts.eps_curv;
        
    end
    
    optim_type = 'gn';
    if(isfield(optim_opts,'optim_type'))
        optim_type = optim_opts.optim_type;
    end
    
    filter_type = 'gauss-conv';
    if(isfield(optim_opts,'filter_type'))
        filter_type = optim_opts.filter_type;
    end
    
    maxit_filter = 10;
    if(isfield(optim_opts,'maxit_filter'))
        maxit_filter = optim_opts.maxit_filter;
    end
    
    verbose = false;
    if(isfield(opts,'verbose'))
        verbose = opts.verbose;
    end
    
    
    % update the geometry part
    
    bc_use = [];
    bc_use.type = bc.type;
    bc_use.invtype = 'o';
    
    frechet_mats = rla.get_frechet_ders(kh,mats,src_info,u_meas,fields,bc_use,opts_use);

    if (constphasefactor)
        [phase,Minv] = optimal_phase_and_jacobian( ...
            u_meas.uscat_tgt(:),fields.uscat_tgt(:),frechet_mats.bdry);
    else
        phase = 1;
        Minv = frechet_mats.bdry;
    end
    
    rhs = u_meas.uscat_tgt(:) - phase*fields.uscat_tgt(:);
    
    opts_update_geom = [];
    opts_update_geom.eps_curv = eps_curv;
    
    n_curv_min = 1;
    if(isfield(optim_opts,'n_curv_min'))
       n_curv_min = optim_opts.n_curv_min;
    end
    
    if(isfield(optim_opts,'n_curv'))
       opts_update_geom.n_curv = optim_opts.n_curv;
    else
        opts_update_geom.n_curv = max(n_curv_min,ncoeff_boundary);
    end
    
    
 
    opts_update_geom.nppw = nppw;
    opts_update_geom.rlam = rlam;
    
    res_in = norm(rhs(:))/norm(u_meas.uscat_tgt(:));
    src_out = src_info;
    mats_out = mats;
    fields_out = fields;
    res = res_in;
    ier = 0;
    
    if(strcmpi(bc.invtype,'o') || strcmpi(bc.invtype,'oi') || strcmpi(bc.invtype,'io'))
        deltas.nmodes_bdry = opts_use.ncoeff_boundary;
        
        Minv = [real(Minv); imag(Minv)];
        
        if(strcmpi(optim_type,'gn') || strcmpi(optim_type,'min(sd,gn)') || strcmpi(optim_type,'min(gn,sd)'))
            delta_bdry_gn = Minv \ [real(rhs(:)); imag(rhs(:))];
            delta_bdry_gn0 = delta_bdry_gn;
            ier_gn = 10;
            iter_filter_bdry_gn = -1;
            for iter=1:maxit_filter
                [src_out_gn,ier_gn] = rla.update_geom(src_info,ncoeff_boundary,delta_bdry_gn,opts_update_geom);
                [mats_out_gn] = rla.get_fw_mats(kh,src_out_gn,bc,u_meas,opts);
                fields_out_gn = rla.compute_fields(kh,src_out_gn,mats_out_gn,u_meas,bc,opts);
                if (constphasefactor)
                    phase = optimal_phase_and_jacobian( ... 
                        u_meas.uscat_tgt(:), fields_out_gn.uscat_tgt(:));
                else
                    phase = 1;
                end
                rhs_gn = u_meas.uscat_tgt(:) - phase*fields_out_gn.uscat_tgt(:);
                res_gn = norm(rhs_gn(:))/norm(u_meas.uscat_tgt(:));
                
                if(ier_gn == 0 && res_gn <= res_in) 
                    iter_filter_bdry_gn = iter-1;
                    break
                else
                    if(strcmpi(filter_type,'gauss-conv'))
                        sigma_bd = 1.0/10^(iter-1);
                        hg = 2/(2*ncoeff_boundary);
                        N_var = ncoeff_boundary;
                        tg = -1:hg:1;
                        gauss_val = exp(-tg.*tg/sigma_bd);
                        gauss_new = zeros(1,length(gauss_val));
                        gauss_new(1:N_var+1) = gauss_val(N_var+1:end);
                        gauss_new(N_var+2:end) = gauss_val(N_var:-1:1);
                        delta_bdry_gn = (delta_bdry_gn0'.*gauss_new)';
                    elseif(strcmpi(filter_type,'step_length'))
                        delta_bdry_gn = delta_bdry_gn0/(2^(iter));
                    end
                end
            end
            
            if(iter_filter_bdry_gn == -1)
                iter_filter_bdry_gn = maxit_filter;
            end
            
            if(verbose)
                fprintf('Post gn filter: Filter type: %s \t filter iteration count: %d \t ier: %d\n',filter_type,iter_filter_bdry_gn,ier_gn);
            end
            if(ier_gn ~=0)
                mats_out_gn = mats;
                fields_out_gn = fields;
                src_out_gn = src_out;
                res_gn = res_in;    
            end
            
            if(res_gn >= res_in)
                mats_out_gn = mats;
                fields_out_gn = fields;
                src_out_gn = src_out;
                res_gn = res_in;
                ier_gn = -5;
            end
        end
        if(strcmpi(optim_type,'sd')  || strcmpi(optim_type,'min(sd,gn)') || strcmpi(optim_type,'min(gn,sd)'))
            delta_bdry_sd = Minv'*[real(rhs(:)); imag(rhs(:))];
            t = delta_bdry_sd'*delta_bdry_sd/norm(Minv*delta_bdry_sd)^2;
            delta_bdry_sd = t*delta_bdry_sd;
            delta_bdry_sd0 = delta_bdry_sd;
            %ier_sd = 10;
            iter_filter_bdry_sd = -1;
            ier_sd = 10;
            for iter=1:maxit_filter
                [src_out_sd,ier_sd] = rla.update_geom(src_info,ncoeff_boundary,delta_bdry_sd,opts_update_geom);
                [mats_out_sd] = rla.get_fw_mats(kh,src_out_sd,bc,u_meas,opts);
                fields_out_sd = rla.compute_fields(kh,src_out_sd,mats_out_sd,u_meas,bc,opts);
                if (constphasefactor)
                    phase = optimal_phase_and_jacobian( ... 
                        u_meas.uscat_tgt(:), fields_out_sd.uscat_tgt(:));
                else
                    phase = 1;
                end
                rhs_sd = u_meas.uscat_tgt(:) - phase*fields_out_sd.uscat_tgt(:);
                res_sd = norm(rhs_sd(:))/norm(u_meas.uscat_tgt(:));
                if(ier_sd == 0 && res_sd<=res_in) 
                    iter_filter_bdry_sd = iter-1;
                    break
                else
                    if(strcmpi(filter_type,'gauss-conv'))
                        sigma_bd = 1.0/10^(iter-1);
                        hg = 2/(2*ncoeff_boundary);
                        N_var = ncoeff_boundary;
                        tg = -1:hg:1;
                        gauss_val = exp(-tg.*tg/sigma_bd);
                        gauss_new = zeros(1,length(gauss_val));
                        gauss_new(1:N_var+1) = gauss_val(N_var+1:end);
                        gauss_new(N_var+2:end) = gauss_val(N_var:-1:1);
                        delta_bdry_sd = (delta_bdry_sd0'.*gauss_new)';
                    elseif(strcmpi(filter_type,'step_length'))
                        delta_bdry_sd = delta_bdry_sd0/(2^(iter));
                    end
                end
            end
            if(iter_filter_bdry_sd == -1)
                iter_filter_bdry_sd = maxit_filter;
            end
               
            
            if(verbose)
                fprintf('Post sd filter: Filter type: %s \t filter iteration count: %d \t ier: %d\n',filter_type,iter_filter_bdry_sd,ier_sd);
            end
            
            if(ier_sd ~=0)
                mats_out_sd = mats;
                fields_out_sd = fields;
                src_out_sd  = src_out;
                res_sd = res_in;
            end
            
            if(res_sd >= res_in)
                mats_out_sd = mats;
                fields_out_sd = fields;
                src_out_sd = src_out;
                res_sd = res_in;
                ier_sd = -5;
            end
        end
        
        
        if(strcmpi(optim_type,'gn'))
            deltas.iter_type = 'gn';
        elseif(strcmpi(optim_type,'sd'))
            deltas.iter_type = 'sd';
        elseif(strcmpi(optim_type,'min(sd,gn)') || strcmpi(optim_type,'min(gn,sd)'))
            if(res_gn < res_sd)
                deltas.iter_type = 'gn';
            else
                deltas.iter_type = 'sd';
            end
            
            fprintf('optimization method used = %s\n',deltas.iter_type);
        end
        
        if(strcmpi(deltas.iter_type,'gn'))
            deltas.delta_bdry = delta_bdry_gn;
            deltas.iter_filter_bdry = iter_filter_bdry_gn;
            mats_out = mats_out_gn;
            fields_out = fields_out_gn;
            
            
            ier = ier_gn;
            src_out = src_out_gn;
            res = res_gn;
        elseif(strcmpi(deltas.iter_type,'sd'))
            deltas.delta_bdry = delta_bdry_sd;
            deltas.iter_filter_bdry = iter_filter_bdry_sd;
            
            mats_out = mats_out_sd;
            fields_out = fields_out_sd;
            ier = ier_sd;
            src_out = src_out_sd;
            res = res_sd;
        end      

        if (constphasefactor)
            phase = optimal_phase_and_jacobian( ... 
            u_meas.uscat_tgt(:), fields_out.uscat_tgt(:));
            deltas.phase = phase;
        else
            deltas.phase = 1;
        end
        
        
    end
    
    % Now update impedance holding boundary fixed
    
    
    if((strcmpi(bc.invtype,'i') || strcmpi(bc.invtype,'oi') || strcmpi(bc.invtype,'io')) && ncoeff_impedance>=0)
        fprintf('Inside update iterate, kh = %d, ncoeff_impedance=%d \n',kh,ncoeff_impedance);
        res_in = res;
        bc_use = [];
        bc_use.type = bc.type;
        bc_use.invtype = 'i';
        
        frechet_mats = rla.get_frechet_ders(kh,mats_out,src_out,u_meas, ...
          fields_out,bc_use,opts_use);
    
        if (constphasefactor)
            [phase,Minv] = optimal_phase_and_jacobian( ...
                u_meas.uscat_tgt(:),fields.uscat_tgt(:),...
                frechet_mats.impedance);
        else
            phase = 1;
            Minv = frechet_mats.impedance;
        end
            
        rhs = u_meas.uscat_tgt(:) - phase*fields_out.uscat_tgt(:);
        if (lambdareal)
            Minv = [real(Minv); imag(Minv)];
            rhs = [real(rhs(:)); imag(rhs(:))];
        end
        if(strcmpi(optim_type,'gn') || strcmpi(optim_type,'min(sd,gn)') || strcmpi(optim_type,'min(gn,sd)'))
            delta_imp_gn = Minv\rhs;
            delta_imp_gn0 = delta_imp_gn;
            ier_gn = 10;
            iter_filter_imp_gn = -1;
            
            for iter=1:maxit_filter
                
                src_out_gn = src_out;
                nh = ncoeff_impedance;
                hcoefs_use = delta_imp_gn(:);
                n = length(src_out.xs);
                t = 0:2*pi/n:2*pi*(1.0-1.0/n); 
                if (nh == 0)
                    h_upd = cos(t*0)*hcoefs_use(1);
                else
                    h_upd = (cos(t'*(0:nh))*hcoefs_use(1:(nh+1)) + ...
                        sin(t'*(1:nh))*hcoefs_use((nh+2):end)).';
                end
                src_out_gn.lambda = src_out_gn.lambda + h_upd.';  
        
                [mats_out_gn] = rla.get_fw_mats(kh,src_out_gn,bc,u_meas,opts);
                fields_out_gn = rla.compute_fields(kh,src_out_gn,mats_out_gn,u_meas,bc,opts);
                if (constphasefactor)
                    phase = optimal_phase_and_jacobian( ... 
                        u_meas.uscat_tgt(:), fields_out_gn.uscat_tgt(:));
                else
                    phase = 1;
                end
                
                rhs_gn = u_meas.uscat_tgt(:) - phase*fields_out_gn.uscat_tgt(:);
                res_imp_gn = norm(rhs_gn(:))/norm(u_meas.uscat_tgt(:));
                
                if(res_imp_gn <= res_in) 
                    iter_filter_imp_gn = iter-1;
                    break
                else
                    if(strcmpi(filter_type,'gauss-conv'))
                        sigma_bd = 1.0/10^(iter-1);
                        hg = 2/(2*ncoeff_impedance);
                        N_var = ncoeff_impedance;
                        tg = -1:hg:1;
                        gauss_val = exp(-tg.*tg/sigma_bd);
                        gauss_new = zeros(1,length(gauss_val));
                        gauss_new(1:N_var+1) = gauss_val(N_var+1:end);
                        gauss_new(N_var+2:end) = gauss_val(N_var:-1:1);
                        delta_imp_gn = (delta_imp_gn0.'.*gauss_new).';
                    elseif(strcmpi(filter_type,'step_length'))
                        delta_imp_gn = delta_imp_gn0/(2^(iter));
                    end
                end
            end
            
            if(iter_filter_imp_gn == -1)
                iter_filter_imp_gn = maxit_filter;
            end
            
            if(verbose)
                fprintf('Post gn filter: Filter type: %s \t filter iteration count: %d \t ier: %d\n',filter_type,iter_filter_imp_gn,ier_gn);
            end
            
            
            if(res_imp_gn >= res_in)
                mats_out_gn = mats_out;
                fields_out_gn = fields_out;
                res_imp_gn = res_in;
                %ier_imp_gn = -5;
            end
        end
        if(strcmpi(optim_type,'sd')  || strcmpi(optim_type,'min(sd,gn)') || strcmpi(optim_type,'min(gn,sd)'))
            delta_imp_sd = Minv'*rhs;
            
            t = delta_imp_sd'*delta_imp_sd/norm(Minv*delta_imp_sd)^2;
            delta_imp_sd = t*delta_imp_sd;
            delta_imp_sd0 = delta_imp_sd;
            %ier_sd = 10;
            iter_filter_imp_sd = -1;
            ier_sd = 10;
            
            for iter=1:maxit_filter
                
                src_out_sd = src_out;
                nh = ncoeff_impedance;
                hcoefs_use = delta_imp_sd(:);
                n = length(src_out.xs);
                t = 0:2*pi/n:2*pi*(1.0-1.0/n); 
                if (nh == 0)
                    h_upd = cos(t*0)*hcoefs_use(1);
                else
                    h_upd = (cos(t'*(0:nh))*hcoefs_use(1:(nh+1)) + sin(t'*(1:nh))*hcoefs_use((nh+2):end)).';
                end
                src_out_sd.lambda = src_out_sd.lambda + h_upd.';  
        
                [mats_out_sd] = rla.get_fw_mats(kh,src_out_sd,bc,u_meas,opts);
                fields_out_sd = rla.compute_fields(kh,src_out_sd,mats_out_sd,u_meas,bc,opts);
                if (constphasefactor)
                    phase = optimal_phase_and_jacobian( ... 
                        u_meas.uscat_tgt(:), fields_out_sd.uscat_tgt(:));
                else
                    phase = 1;
                end
                
                rhs_sd = u_meas.uscat_tgt(:) - phase*fields_out_sd.uscat_tgt(:);
                res_imp_sd = norm(rhs_sd(:))/norm(u_meas.uscat_tgt(:));
                if(res_imp_sd<=res_in) 
                    iter_filter_imp_sd = iter-1;
                    break
                else
                    if(strcmpi(filter_type,'gauss-conv'))
                        sigma_bd = 1.0/10^(iter-1);
                        hg = 2/(2*ncoeff_impedance);
                        N_var = ncoeff_impedance;
                        tg = -1:hg:1;
                        gauss_val = exp(-tg.*tg/sigma_bd);
                        gauss_new = zeros(1,length(gauss_val));
                        gauss_new(1:N_var+1) = gauss_val(N_var+1:end);
                        gauss_new(N_var+2:end) = gauss_val(N_var:-1:1);
                        delta_imp_sd = (delta_imp_sd0.'.*gauss_new).';
                    elseif(strcmpi(filter_type,'step_length'))
                        delta_imp_sd = delta_imp_sd0/(2^(iter));
                    end
                end
            end
            if(iter_filter_imp_sd == -1)
                iter_filter_imp_sd = maxit_filter;
            end
               
            
            if(verbose)
                fprintf('Post sd filter: Filter type: %s \t filter iteration count: %d \t ier: %d\n',filter_type,iter_filter_imp_sd,ier_sd);
            end
            
            
            
            if(res_imp_sd >= res_in)
                mats_out_sd = mats_out;
                fields_out_sd = fields_out;
                res_imp_sd = res_in;
                %ier_imp_sd = -5;
                
            end
        end
        
        
        if(strcmpi(optim_type,'gn'))
            deltas.iter_imp_type = 'gn';
        elseif(strcmpi(optim_type,'sd'))
            deltas.iter_imp_type = 'sd';
        elseif(strcmpi(optim_type,'min(sd,gn)') || strcmpi(optim_type,'min(gn,sd)'))
            if(res_imp_gn < res_imp_sd)
                deltas.iter_imp_type = 'gn';
            else
                deltas.iter_imp_type = 'sd';
            end
            
            fprintf('optimization method used = %s\n',deltas.iter_imp_type);
        end
        
        if(strcmpi(deltas.iter_imp_type,'gn'))
            deltas.delta_impedance = delta_imp_gn;
            deltas.iter_filter_impedance = iter_filter_imp_gn;
            mats_out = mats_out_gn;
            fields_out = fields_out_gn;
            src_out = src_out_gn;
            res = res_imp_gn;
            
        elseif(strcmpi(deltas.iter_imp_type,'sd'))
            deltas.delta_impedance = delta_imp_sd;
            deltas.iter_filter_impedance = iter_filter_imp_sd;
            
            mats_out = mats_out_sd;
            fields_out = fields_out_sd;
            src_out = src_out_sd;
            res = res_imp_sd;
        end      
        if (constphasefactor)
            deltas.phase = optimal_phase_and_jacobian( ... 
                u_meas.uscat_tgt(:), fields_out.uscat_tgt(:));
            
        else
            deltas.phase = 1;
        end
        
    end 
    
end
