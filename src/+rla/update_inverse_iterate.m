function [deltas,src_out,mats_out,fields_out,res,ier_obs,ier_imp] = ...
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
%      src_info.lambda - imepedance value at discretization nodes 
%           (only if solving impedance boundary value problem with
%            opts.impedance_type = 'fourier')
%      src_info.lamcfs - imepedance coefficients in constant + kappa model
%           (only if solving impedance boundary value problem with
%            opts.impedance_type = 'constkappa') or one of the
%            antoine-barucq type models ('antbar2','antbar3')
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
%      opts.impedance_type: character array or string indicating
%         which representation of impedance to use 
%             "fourier" (default) - Fourier series with ncoeff_impedance
%             terms used 
%             "constkappa" - a function of the form c1 + c2*kappa is used, 
%             where kappa is the signed curvature     
%             "antbar" - a function of the form c1 + c2*kappa is used, 
%             where kappa is the signed curvature     
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
    if(~isfield(opts,'impedance_type'))
        impedance_type = 'fourier';
        opts_use.impedance_type = impedance_type;
    else
        impedance_type = opts.impedance_type;
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

    if (strcmpi(impedance_type,'antbar2') || strcmpi(impedance_type,'antbar3'))
        if ~lambdareal
            warning('for antbar models, lambda is real. overriding...')
            lambdareal = true;
        end
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

    if (strcmpi(optim_type,'min(gn,sd)') || ...
            strcmpi(optim_type,'min(sd,gn)'))
        optim_list = {'gn','sd'};
    else
        optim_list = {optim_type};
    end

    optim_list_imp = optim_list;
    optim_type_imp = optim_type;
    if(isfield(optim_opts,'optim_type_imp'))
        optim_type_imp = optim_opts.optim_type_imp;
        if (strcmpi(optim_type_imp,'min(gn,sd)') || ...
                strcmpi(optim_type_imp,'min(sd,gn)'))
            optim_list_imp = {'gn','sd'};
        else
            optim_list_imp = {optim_type_imp};
        end
    end

    stepfac = 2;
    if isfield(optim_opts,'stepfac')
        stepfac = optim_opts.stepfac;
    end

    
    filter_type = 'gauss-conv';
    if(isfield(optim_opts,'filter_type'))
        filter_type = optim_opts.filter_type;
    end

    if (strcmpi(filter_type,'min(step_length,gauss-conv)') || ...
            strcmpi(filter_type,'min(gauss-conv,step_length)'))
        filter_list = {'gauss-conv','step_length'};
    else
        filter_list = {filter_type};
    end
    
    % default is to use the same as obstacle
    filter_type_imp = filter_type;
    filter_list_imp = filter_list;
    if (strcmpi(impedance_type,'antbar2') || strcmpi(impedance_type,'antbar3') || ...
            strcmpi(impedance_type,'constkappa'))
        % for constkappa type models default is just step length
        filter_type_imp = 'step_length';
        filter_list_imp = {filter_type_imp};
    end

    if isfield(optim_opts,'filter_type_imp')
        filter_type_imp = optim_opts.filter_type_imp;
        if (strcmpi(filter_type_imp,'min(step_length,gauss-conv)') || ...
                strcmpi(filter_type_imp,'min(gauss-conv,step_length)'))
            filter_list_imp = {'gauss-conv','step_length'};
        else
            filter_list_imp = {filter_type_imp};
        end
    end        

    if (strcmpi(impedance_type,'antbar2') || strcmpi(impedance_type,'antbar3') || ...
            strcmpi(impedance_type,'constkappa'))
        if (length(filter_list_imp) > 1 || ~strcmpi(filter_list_imp{1},'step_length'))
            warning('for constkappa type models, only step_length filtering makes sense. overriding');
            filter_list_imp = {'step_length'};
        end
    end
    

    maxit_filter = 10;
    if(isfield(optim_opts,'maxit_filter'))
        maxit_filter = optim_opts.maxit_filter;
    end
    
    verbose = false;
    if(isfield(opts,'verbose'))
        verbose = opts.verbose;
    end

    %
    opts_update_geom = [];
    opts_update_geom.eps_curv = eps_curv;
    
    n_curv_min = 1;
    if(isfield(optim_opts,'n_curv_min'))
       n_curv_min = optim_opts.n_curv_min;
    end

    % ninner = number of impedance solves per boundary solve
    ninner = 1;
    if (isfield(optim_opts,'ninner'))
        ninner = optim_opts.ninner;
    end
    
    if(isfield(optim_opts,'n_curv'))
       opts_update_geom.n_curv = optim_opts.n_curv;
    else
        opts_update_geom.n_curv = max(n_curv_min,ncoeff_boundary);
    end
    
    
    opts_update_geom.nppw = nppw;
    opts_update_geom.rlam = rlam;

    opts_update_geom.impedance_type = impedance_type;
    opts_update_geom.kh = kh;
    
    % update the geometry part

    bc_use = [];
    bc_use.type = bc.type;
    bc_use.invtype = 'o';
    
    frechet_mats = rla.get_frechet_ders(kh,mats,src_info,u_meas,...
        fields,bc_use,opts_use);

    if (constphasefactor)
        [phase,Minv,Minvbar] = optimal_const_and_jacobian( ...
            u_meas.uscat_tgt(:),fields.uscat_tgt(:),frechet_mats.bdry);
    else
        phase = 1;
        Minv = frechet_mats.bdry;
        Minvbar = zeros(size(Minv));
    end
    
    rhs = u_meas.uscat_tgt(:) - phase*fields.uscat_tgt(:);
    
    res = norm(rhs(:))/norm(u_meas.uscat_tgt(:));
    res0 = res;
    src_out = src_info;
    mats_out = mats;
    fields_out = fields;
    ier_obs = 0;
    ier_imp = 0;
    Minv = [real(Minv+Minvbar); imag(Minv+Minvbar)];
        
        
    if(strcmpi(bc.invtype,'o') || strcmpi(bc.invtype,'oi') || strcmpi(bc.invtype,'io'))
        deltas.nmodes_bdry = opts_use.ncoeff_boundary;
        
        ier_obs = -5;
        
        % loop over optimization methods
        for ioptim = 1:length(optim_list)
            optim_type0 = optim_list{ioptim};
            if (strcmpi(optim_type0,'gn'))
                delta_bdry = Minv \ [real(rhs(:)); imag(rhs(:))];
                delta_bdry0 = delta_bdry;
            elseif (strcmpi(optim_type0,'sd'))
                delta_bdry = Minv'*[real(rhs(:)); imag(rhs(:))];
                t = delta_bdry'*delta_bdry/norm(Minv*delta_bdry)^2;
                delta_bdry = t*delta_bdry;
                delta_bdry0 = delta_bdry;
            else
                error('unknown optimization type %s',optim_type0);
            end

            % loop over filtration strategies
            for ifilt = 1:length(filter_list)
                filter_type0 = filter_list{ifilt};
                iter_filter_bdry_tmp = -1;

                delta_bdry = delta_bdry0;

                for iter=1:maxit_filter
                    [src_out_tmp,ier_tmp] = rla.update_geom(src_info,ncoeff_boundary,delta_bdry,opts_update_geom);
                    [mats_out_tmp] = rla.get_fw_mats(kh,src_out_tmp,bc,u_meas,opts);
                    fields_out_tmp = rla.compute_fields(kh,src_out_tmp,mats_out_tmp,u_meas,bc,opts);
                    if (constphasefactor)
                        phase = optimal_const_and_jacobian( ... 
                            u_meas.uscat_tgt(:), fields_out_tmp.uscat_tgt(:));
                    else
                        phase = 1;
                    end
                    rhs_tmp = u_meas.uscat_tgt(:) - phase*fields_out_tmp.uscat_tgt(:);
                    res_tmp = norm(rhs_tmp(:))/norm(u_meas.uscat_tgt(:));
                
                    if(ier_tmp == 0 && res_tmp < res0) 
                        iter_filter_bdry_tmp = iter-1;
                        break
                    else
                        if(strcmpi(filter_type0,'gauss-conv'))
                            sigma_bd = 1.0/10^(iter-1);
                            hg = 2/(2*ncoeff_boundary);
                            N_var = ncoeff_boundary;
                            tg = -1:hg:1;
                            gauss_val = exp(-tg.*tg/sigma_bd);
                            gauss_new = zeros(1,length(gauss_val));
                            gauss_new(1:N_var+1) = gauss_val(N_var+1:end);
                            gauss_new(N_var+2:end) = gauss_val(N_var:-1:1);
                            delta_bdry = (delta_bdry0.'.*gauss_new).';
                        elseif(strcmpi(filter_type0,'step_length'))
                            delta_bdry = delta_bdry0/(stepfac^(iter));
                        else
                            error('unknown filter type %s',filter_type0);
                        end
                    end
                end
                if (iter_filter_bdry_tmp == -1)
                    iter_filter_bdry_tmp = maxit_filter;
                end
                if(verbose)
                    fprintf(['Post %s filter: Filter type: %s \t ', ...
                         'filter iteration count: %d \t ier: %d\n'], ...
                         optim_type0,filter_type0,iter_filter_bdry_tmp,ier_tmp);
                end

                if res_tmp < res
                    mats_out = mats_out_tmp;
                    src_out = src_out_tmp;
                    fields_out = fields_out_tmp;
                    res = res_tmp;
                    ier_obs = 0;
                    deltas.delta_bdry = delta_bdry;
                    deltas.iter_filter_bdry = iter_filter_bdry_tmp;
                    deltas.filter_type = filter_type0;
                    deltas.iter_type = optim_type0;
                    deltas.phase = phase;
                end

                if (iter_filter_bdry_tmp == 0)
                    break
                end

            end
        end

        if verbose
            if ier_obs == 0
                msg = ['Post update geometry. Success. Used optim type: %s ' ...
                    '\t filter type: %s\n'];
                fprintf(msg,deltas.iter_type,deltas.filter_type)
            else
                msg = 'Post update geometry. Failed to update\n';
                fprintf(msg);
            end
        end
    end
    
    % Now update impedance holding boundary fixed
    
    if((strcmpi(bc.invtype,'i') || strcmpi(bc.invtype,'oi') || strcmpi(bc.invtype,'io')) && ...
            (ncoeff_impedance>=0 || strcmpi(impedance_type,'constkappa')))
        res0 = res;
        for jinner = 1:ninner
            ier_imp = -5;

            if (strcmpi(impedance_type,'fourier'))
                fprintf('Inside update iterate, kh = %d, ncoeff_impedance=%d \n',kh,ncoeff_impedance);
            else
                fprintf('Inside update iterate, kh = %d, fitting constant plus curvature impedance\n',kh);
                disp(src_out.lamcfs)
            end
    
            bc_use = [];
            bc_use.type = bc.type;
            bc_use.invtype = 'i';
            
            frechet_mats = rla.get_frechet_ders(kh,mats_out,src_out,u_meas, ...
              fields_out,bc_use,opts_use);
    
            if (constphasefactor)
                [phase,Minv,Minvbar] = optimal_const_and_jacobian( ...
                    u_meas.uscat_tgt(:),fields.uscat_tgt(:),...
                    frechet_mats.impedance);
            else
                phase = 1;
                Minv = frechet_mats.impedance;
                Minvbar = zeros(size(Minv));
            end
                
            rhs = u_meas.uscat_tgt(:) - phase*fields_out.uscat_tgt(:);
            if (lambdareal)
                Minv = [real(Minv+Minvbar); imag(Minv+Minvbar)];
                rhs = [real(rhs(:)); imag(rhs(:))];
            else
                Minv = [real(Minv+Minvbar), imag(Minvbar-Minv); ...
                    imag(Minv+Minvbar), real(Minv-Minvbar)];
                rhs = [real(rhs(:)); imag(rhs(:))];
            end

            for ioptim = 1:length(optim_list_imp)
                optim_type0 = optim_list_imp{ioptim};
                if(strcmpi(optim_type0,'gn')) 
    
                    if (strcmpi(impedance_type,'constkappa') && cond(frechet_mats.impedance) > 10^3)
                        % if ill conditioned, revert to steepest descent
                        delta_imp_tmp = Minv'*rhs;
                        t = delta_imp_tmp'*delta_imp_tmp/norm(Minv*delta_imp_tmp)^2;
                        delta_imp_tmp = t*delta_imp_tmp;
                        fprintf(['bad basis for impedance, switching ', ...
                            'to sd...\n'])
                    else
                        delta_imp_tmp = Minv\rhs;
                    end
                elseif(strcmpi(optim_type0,'sd')) 
                    delta_imp_tmp = Minv'*rhs;
                    t = delta_imp_tmp'*delta_imp_tmp/norm(Minv*delta_imp_tmp)^2;
                    delta_imp_tmp = t*delta_imp_tmp;
                else
                    error('optim_type %s not recognized',optim_type0)
                end
                if (~lambdareal)
                    delta_imp_tmp = delta_imp_tmp(1:end/2) + ...
                        1i*delta_imp_tmp(end/2+1:end);   
                end
                delta_imp_tmp0 = delta_imp_tmp;
                ier_imp_tmp = 0;
                iter_filter_imp_tmp = -1;
                
                for ifilt = 1:length(filter_list_imp)
                    filter_type0 = filter_list_imp{ifilt};
                    delta_imp_tmp = delta_imp_tmp0;
                    for iter=1:maxit_filter
                        
                        src_out_tmp = src_out;
                        if strcmpi(impedance_type,'fourier')
                            nh = ncoeff_impedance;
                            hcoefs_use = delta_imp_tmp(:);
                            n = length(src_out.xs);
                            t = 0:2*pi/n:2*pi*(1.0-1.0/n); 
                            if (nh == 0)
                                h_upd = cos(t*0)*hcoefs_use(1);
                            else
                                h_upd = (cos(t'*(0:nh))*hcoefs_use(1:(nh+1)) + ...
                                    sin(t'*(1:nh))*hcoefs_use((nh+2):end)).';
                            end
                            src_out_tmp.lambda = src_out_tmp.lambda + h_upd.';  
                        elseif ( strcmpi(impedance_type,'constkappa') || ...
                            strcmpi(impedance_type,'antbar2') || ...
                            strcmpi(impedance_type,'antbar3') )
    
                            n = length(src_out.xs);
                            src_out_tmp.lamcfs(1:length(delta_imp_tmp)) = ...
                                src_out_tmp.lamcfs(1:length(delta_imp_tmp)) ...
                                + delta_imp_tmp(:);
                            if (isfield(opts,'lambdaprox'))
                                src_out_tmp.lamcfs = opts.lambdaprox(src_out_tmp.lamcfs);
                            end
                            ckcoefs = constkappa_models_convert(src_out_tmp.lamcfs,...
                                impedance_type,kh);
                            src_out_tmp.lambda(:) = [ones(n,1) src_out_tmp.H(:)]*ckcoefs;
                        else
                            error('unknown impedance_type %s',impedance_type);
                        end
    
                
                        [mats_out_tmp] = rla.get_fw_mats(kh,src_out_tmp,bc,u_meas,opts);
                        fields_out_tmp = rla.compute_fields(kh,src_out_tmp,mats_out_tmp,u_meas,bc,opts);
                        if (constphasefactor)
                            phase = optimal_const_and_jacobian( ... 
                                u_meas.uscat_tgt(:), fields_out_tmp.uscat_tgt(:));
                        else
                            phase = 1;
                        end
                        
                        rhs_tmp = u_meas.uscat_tgt(:) - phase*fields_out_tmp.uscat_tgt(:);
                        res_imp_tmp = norm(rhs_tmp(:))/norm(u_meas.uscat_tgt(:));
                        
                        if(res_imp_tmp < res0)
                            iter_filter_imp_tmp = iter-1;
                            break
                        else
                            if(strcmpi(filter_type0,'gauss-conv'))
                                sigma_bd = 1.0/10^(iter-1);
                                hg = 2/(2*ncoeff_impedance);
                                N_var = ncoeff_impedance;
                                tg = -1:hg:1;
                                gauss_val = exp(-tg.*tg/sigma_bd);
                                gauss_new = zeros(1,length(gauss_val));
                                gauss_new(1:N_var+1) = gauss_val(N_var+1:end);
                                gauss_new(N_var+2:end) = gauss_val(N_var:-1:1);
                                delta_imp_tmp = (delta_imp_tmp0.'.*gauss_new).';
                            elseif(strcmpi(filter_type0,'step_length'))
                                delta_imp_tmp = delta_imp_tmp0/(stepfac^(iter));
                            end
                        end
                    end
                    
    
                    if(iter_filter_imp_tmp == -1)
                         iter_filter_imp_tmp = maxit_filter;
                    end
                
                    if(res_imp_tmp < res)
                        deltas.delta_impedance = delta_imp_tmp;
                        deltas.iter_filter_impedance = iter_filter_imp_tmp;
                        deltas.iter_imp_type = optim_type0;
                        deltas.phase = phase;
                        mats_out = mats_out_tmp;
                        src_out = src_out_tmp;
                        fields_out= fields_out_tmp;
                        res = res_imp_tmp;
                        ier_imp = 0;
                    end
                    if(verbose)
                        fprintf('Post %s filter: Filter type: %s \t filter iteration count: %d \t ier: %d\n',optim_type0,filter_type0,iter_filter_imp_tmp,ier_imp_tmp);
                    end

                    if (iter_filter_imp_tmp == 0)
                        break
                    end
                end
            end
        end % end of inner iterations loop        
    end % end of impedance section 

