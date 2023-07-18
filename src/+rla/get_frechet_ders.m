function [frechet_mats] = get_frechet_ders(kh,mats,src_info,sensor_info,fields,bc,opts)
%
% get_frechet_ders returns the relevant Frechet derivative for solving
% Dirichlet, Neumann/ Impedance boundary value problem
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
%   bc - boundary condition struct;
%     bc.type = type of boundary condition;
%        'd' or 'Dirichlet' for dirichlet
%        'n' or 'Neumann' for Neumann
%        'i' or 'Impedance' for impedance
%        't' ot 'Transmission'  for transmission
%     the transmission problem requires additional parameters:
%        these are all length 2 arrays, first entry is interior value,
%        second is exterior value
%        bc.transk = transmission wave numbers
%        bc.transa = coefficients for jump in potential
%        bc.transb = coefficients for jump in normal
%     bc.invtype = type of inverse problem;
%        'o' or 'obstacle' for obstacle only
%        'oi' or 'obsctacle and impedance' for both obstacle and impedance;
%
%   opts - options struct
%      opts.ifflam - whether to use FLAM to compress the matrices (false)
%      opts.verbose - flag for displaying detailed output (false)
%      opts.ncoeff_boundary - number of coeffs to use for boundary update (floor(kh*2))
%      opts.ncoeff_impedance - number of coeffs to use for impedance update
%      (floor(kh*0.5))
%      
%      
%      
% Output
%    frechet_mats - Frechet derivative matrix
%        frechet_mats.bdry: Frechet derivative for boundary update
%        frechet_mats.impedance: Frechet derivative for impedance updates

   if(nargin < 7)
       opts = [];
   end
   ifflam = false;
   if(isfield(opts,'ifflam'))
       ifflam = opts.ifflam;
   end

   verbose = false;
   if(isfield(opts,'verbose'))
       verbose = opts.verbose;
   end
   nc_b= floor(kh*2);
   if(isfield(opts,'ncoeff_boundary'))
       nc_b = opts.ncoeff_boundary;
   end
   
   nc_i= floor(kh*0.5);
   if(isfield(opts,'ncoeff_impedance'))
       nc_i = opts.ncoeff_impedance;
   end
   
   u = fields.uinc + fields.uscat;
   dudn = fields.dudninc + fields.dudnscat;
   [~,n_dir] = size(dudn);
   du = src_info.Der*u;
   m = length(fields.uscat_tgt(:));
   
   n = length(src_info.xs);
   t = 0:(2*pi/n):2*pi*(1.0-1.0/n);
   
   [t_dir_uni,~,idir] = unique(sensor_info.t_dir);
   [tgt_uni,~,itgt] = unique(sensor_info.tgt','rows');
   nt_uni = length(tgt_uni(:,1));
   induse = itgt + (idir-1)*nt_uni;

   if(strcmpi(bc.invtype,'o') || strcmpi(bc.invtype,'oi') || strcmpi(bc.invtype,'io'))
     DFw_boundary = complex(zeros(m,2*nc_b+1));
     
     %%%%%%%%%%%%%Here goes the code to generate J with a larger amoun of
     %%%%%%%%%%%%%modes at specific
% % %      if ((kh == 1.0)||(kh == 5.0)||(kh == 10.0)||(kh == 15.0)||(kh == 20.0)||(kh == 25.0)||(kh == 30.0))
% % %          
% % %          Ntot = floor(10*kh);
% % %          for ivar = 1 : ( 2*Ntot + 1 )             
% % %             if(ivar<=Ntot+1)
% % %                  h_upd = cos((ivar-1)*t);
% % %             else
% % %                  h_upd = sin((ivar-Ntot-1)*t);
% % %             end
% % %             h_upd_nu = repmat(h_upd',1,n_dir);
% % %             bd_data_delta = -h_upd_nu.*dudn;
% % % 
% % %             if(~ifflam)
% % %                  DFw_col = mats.bdrydata_to_receptor*bd_data_delta;
% % %             end
% % %             DFw_mat(:,ivar) = DFw_col(induse); 
% % %             
% % %          end
% % % %          size(DFw_mat)
% % %          
% % %          for ivar = 1: Ntot
% % %             vec_cond(ivar) = cond([DFw_mat(:,1) DFw_mat(:,2:ivar+1) DFw_mat(:,Ntot+2:Ntot+1+ivar)]);
% % %          end
% % %          
% % %          switch kh
% % %             case 1.0    
% % %                  save('cond_number_J_k1.mat','vec_cond','-v7.3')            
% % %             case 5.0
% % %                  save('cond_number_J_k5.mat','vec_cond','-v7.3')            
% % %             case 10.0
% % %                  save('cond_number_J_k10.mat','vec_cond','-v7.3')            
% % %             case 15.0
% % %                  save('cond_number_J_k15.mat','vec_cond','-v7.3')            
% % %             case 20.0
% % %                  save('cond_number_J_k20.mat','vec_cond','-v7.3')            
% % %             case 25.0
% % %                  save('cond_number_J_k25.mat','vec_cond','-v7.3')            
% % %             case 30.0
% % %                  save('cond_number_J_k30.mat','vec_cond','-v7.3')            
% % %              otherwise
% % % 
% % %          end
% % %      end
     %%%%%%%%%%%%%%%
     
     for ivar=1:(2*nc_b + 1)
         if(ivar<=nc_b+1)
             h_upd = cos((ivar-1)*t);
         else
             h_upd = sin((ivar-nc_b-1)*t);
         end
         h_upd_nu = repmat(h_upd',1,n_dir);

         
         if(strcmpi(bc.type,'d') || strcmpi(bc.type,'Dirichlet'))
             bd_data_delta = -h_upd_nu.*dudn;
         end
         
         if(strcmpi(bc.type,'n') || strcmpi(bc.type,'Neumann'))
             bd_data_delta = kh^2*h_upd_nu.*u + src_info.Der*(h_upd_nu.*du);
         end
         
         if(strcmpi(bc.type,'i') || strcmpi(bc.type,'Impedance'))
             bd_data_delta = kh^2*h_upd_nu.*u + src_info.Der*(h_upd_nu.*du);
             lambda_rep = repmat(src_info.lambda,1,n_dir);
             bd_data_delta = bd_data_delta - 1i*kh*lambda_rep.*h_upd_nu.*(dudn + repmat(src_info.H',1,n_dir).*u);
         end
         
         if(strcmpi(bc.type,'t') || strcmpi(bc.type,'Transmission'))    
             
             zks = bc.transk; a = bc.transa; b = bc.transb;
             q = (a(1)/b(1) + a(2)/b(2))*0.5;
             kh = zks(2); khi = zks(1);             
             
            %derivative u        
            der_rhs_u  = h_upd_nu.*(dudn - fields.dudnin);
            
            %derivative dudn        
            part1      = h_upd_nu.*(khi^2*fields.uin -kh^2*b(2)*u);                
            part2      = src_info.Der*(h_upd_nu.*(src_info.Der*(fields.uin-b(2)*u)));            
            der_rhs_du = part1 + part2;

            bd_data_delta = zeros(2*size(der_rhs_u,1),n_dir,'like',1.0+1i);
            bd_data_delta(2:2:end) = - der_rhs_du;
            bd_data_delta(1:2:end) = - der_rhs_u/q;
            
         end        
         
         if(~ifflam)
            DFw_col = mats.bdrydata_to_receptor*bd_data_delta;
         end
         DFw_boundary(:,ivar) = DFw_col(induse); 
     end
     frechet_mats.bdry = DFw_boundary;
   end
   

   if(strcmpi(bc.invtype,'i') || strcmpi(bc.invtype,'oi') || strcmpi(bc.invtype,'io'))
     DFw_impedance = complex(zeros(m,2*nc_i+1));
     for ivar=1:(2*nc_i + 1)
         if(ivar<=nc_i+1)
             h_upd = cos((ivar-1)*t);
         else
             h_upd = sin((ivar-nc_i-1)*t);
         end  
         
         if(strcmpi(bc.type,'i') || strcmpi(bc.type,'Impedance'))
             delta_lambda_rep = repmat(h_upd.',1,n_dir);
             bd_data_delta = -1i*kh*delta_lambda_rep.*u;
         end
         
         
         if(~ifflam)
            DFw_col = mats.bdrydata_to_receptor*bd_data_delta;
         end
         DFw_impedance(:,ivar) = DFw_col(induse); 
     end
     frechet_mats.impedance = DFw_impedance;
   end
end