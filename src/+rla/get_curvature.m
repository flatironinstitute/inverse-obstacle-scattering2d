function [H] = get_curvature(src_info)
%  This subroutine evaluates the curvature for given source struct
%
%
   n_bd = length(src_info.xs);
   if(~isfield(src_info,'Der_param'))
       if(isfield(src_info,'paramL'))
          rsc = 2*pi/src_info.paramL;
       else
          rsc = 1.0;
       end
       Der = specdiffmat_ds(n_bd,ones(size(src_info.ds)))*rsc;
   else
       Der = src_info.Der_param;
   end
   H = ( src_info.dxs .* (Der*src_info.dys')' - src_info.dys .* (Der*src_info.dxs')' )  ./ ( src_info.ds.^3 );
       
end