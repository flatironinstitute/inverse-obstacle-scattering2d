function fields = compute_fields(kh,src_info,mats,sensor_info,bc,opts)
%
% This subroutine computes the scattered field and it's normal derivative
% on the boundary of the obstacle and the scattered field
% at a collection of target locations due to the incident directions
% prescribed by the user, and stores both the incident and the scattered
% fields along with their normal derivatives on the boundary, and also
% the scattered field at the sensor locations
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
%   sensor_info - sensor information struct
%      sensor_info.tgt(2,nmeas) - xy cooordinates of sensors
%         sensor_info.tgt(1:2,i) = xy coordinates corresponding the ith
%            measurement
%      sensor_info.t_dir(nmeas) - incident directions
%         sensor_info.t_dir(i) - is the incident direction corresponding to
%            the ith measurement
%   bc - boundary condition struct;
%     bc.type = type of boundary condition;
%        'd' or 'Dirichlet' for dirichlet
%        'n' or 'Neumann' for Neumann
%        'i' or 'Impedance' for impedance
%        't' or 'Transmission' for transmission
%     the transmission problem requires additional parameters:
%        these are all length 2 arrays, first entry is interior value,
%        second is exterior value
%        bc.transk = transmission wave numbers
%        bc.transa = coefficients for jump in potential
%        bc.transb = coefficients for jump in normal
  
%     bc.invtype = type of inverse problem;
%        'o' or 'obstacle' for obstacle only
%        'oi' or 'obsctacle and impedance' for both obstacle and impedance;
% Optional input arguments
%   opts - options struct
%     opts.ifflam - whether fast direct solver was used to compress
%        matrices and their inverses in the mats structure (false)
%     NOTE: This option is currently unavailable
%
% Output arguments:
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

   fields = [];
   
   if(nargin < 6)
       opts = [];
   end
   
   ifflam = false;
   if(isfield(opts,'ifflam'))
       ifflam = opts.ifflam;
   end
   
   [t_dir_uni,~,idir] = unique(sensor_info.t_dir);
   [tgt_uni,~,itgt] = unique(sensor_info.tgt','rows');
   nt_uni = length(tgt_uni(:,1));
   induse = itgt + (idir-1)*nt_uni;
   
   
   t_dir_uni = t_dir_uni(:);
   x_dir = cos(t_dir_uni)';
   y_dir = sin(t_dir_uni)';
   n_dir = length(x_dir);
   xs = src_info.xs(:)';
   ys = src_info.ys(:)';
   ds = src_info.ds(:)';
   dxs = src_info.dxs(:)';
   dys = src_info.dys(:)';
   
   fields.uinc  = exp(1i *kh * (bsxfun(@times,xs',x_dir)+bsxfun(@times,ys',y_dir)));
   fields.dudninc = 1i* kh * (bsxfun(@times,dys',x_dir)-bsxfun(@times,dxs',y_dir))./repmat(ds',1,n_dir) .* ...
            exp(1i *kh * (bsxfun(@times,xs',x_dir)+bsxfun(@times,ys',y_dir)));
   
   if(strcmpi(bc.type,'d') || strcmpi(bc.type,'Dirichlet'))
      bd_data = -fields.uinc;
   end

   if(strcmpi(bc.type,'n') || strcmpi(bc.type,'Neumann'))
      bd_data = -fields.dudninc;
   end

   if(strcmpi(bc.type,'i') || strcmpi(bc.type,'Impedance'))
     lambda_rep = repmat(src_info.lambda(:),1,n_dir);
     bd_data = -(fields.dudninc + 1i*kh*lambda_rep.*fields.uinc);
   end

   if(strcmpi(bc.type,'t') || strcmpi(bc.type,'Transmission'))
     a = bc.transa; b = bc.transb; zks = bc.transk;
     q = 0.5*(a(1)/a(2) + a(2)/b(2));
     n1 = numel(xs);
     bd_data = zeros(2*n1,n_dir,'like',1.0+1i);
     bd_data(1:2:end,:) = -a(2)*fields.uinc/q;
     bd_data(2:2:end,:) = -b(2)*fields.dudninc;
   end
   
   if(~ifflam)
     uscat_tgt_all = mats.bdrydata_to_receptor*bd_data;
     uscat_tgt_all = uscat_tgt_all(:);
     fields.uscat_tgt = uscat_tgt_all(induse);
     sigma = mats.inv_Fw_mat*bd_data;
     fields.uscat = mats.Fw_dir_mat*sigma;

     if(strcmpi(bc.type,'d') || strcmpi(bc.type,'Dirichlet'))
         bd_data2 = fields.dudninc - 1i * kh *fields.uinc;
         fields.dudnscat = mats.Fw_neu_mat*bd_data2 - fields.dudninc;
     end
     
     if(strcmpi(bc.type,'n') || strcmpi(bc.type,'Neumann') || strcmpi(bc.type,'i') || strcmpi(bc.type,'Impedance'))
        fields.dudnscat = mats.Fw_neu_mat*sigma; 
     end
   end
   
end
