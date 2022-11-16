function [Ig,xgrid0,ygrid0] = lsm_tensor(n_tgt,n_dir,u_meas,alpha)
% This function evaluates the level set using the 
% linear sampling method. The code currently assumes
% that the obstacle is a Dirichlet obstacle, that the data
% is generated in a tensor product of sensors and incident
% directions, and that u_meas.uscat_tgt stores the measured
% data at the sensors ordered in (n_dir,n_tgt) ordering
% and that the sensors are located at a radius of 10
%
% The object is assumed to be contained in [-1.5,1.5], 
% and the target grid is between [-3,3];
%
%
%  Input arugments:
%    n_tgt: number of equispaced targets on the radius 10 disk
%    n_dir: number of equispaced incident directions, full aperture
%    u_meas: measurement struct
%        u_meas.uscat_tgt - scattered field at target locations, due
%         to incident waves from all directions, ordered as
%         uscat_tgt(n_dir,n_tgt)
%        u_meas.kh - wavenumber
%    alpha: regularization parameter (Recommended value 1e-3)
%
%  Output arguments:
%    Ig: level set function
%    [xgrid0,ygrid0]: mesh grid version of points in volume grid where
%       level set function is evaluated

    Nd = n_tgt;
    hd = 2*pi/Nd;
    t_tgt = 0:hd:2*pi-hd;

    % sensor locations
    radius = 10;
    tgt = radius * [cos(t_tgt);sin(t_tgt)];


    kh = u_meas.kh;
    uscat = u_meas.uscat_tgt;
    uscat = reshape(uscat,[n_dir,n_tgt]);
    uscat = uscat.';
    
    tgrid0 = -3;
    tgrid1 = 3;
    nppw = 40;
    Ngrid = max(ceil(nppw*6*abs(kh)/2/pi),100);
    % setting up grid for inverse problem
    hgrid = (tgrid1 -tgrid0)/(Ngrid-1);
    tgrid = tgrid0:hgrid:tgrid1;
    [xgrid0,ygrid0] = meshgrid(tgrid);
    xgrid = xgrid0(:);
    ygrid = ygrid0(:);

    % setting up rhs for LSM
    b = exp(1i*pi/4)*sqrt(pi*kh/2)*besselh(0,1,kh*sqrt( ...
                                  bsxfun(@minus,xgrid,tgt(1,:)).^2 + ...
                                  bsxfun(@minus,ygrid,tgt(2,:)).^2));
    b = transpose(b);

    % setting-up operator
    % alpha -> regularization parameter (you can play around, 10^-6 does the
    % trick here for the kite
    %
    A = uscat*hd;
    F = A' * A + alpha * eye(n_dir);
    fprintf('Condition number of LSM operator:%d\n',cond(A))
    fprintf('Condition number of LSM problem after regularization:%d\n',cond(F))

    % Solving for LSM
    % For each point you will have a vector. g is a matrix Nd times Ngrid
    g = F \ (A' * b);

    % Calculating the indicator function at each point in the grid. 
    % Just like the paper above we have
    % Ig(z) = log(\|g(z)\|_L^2). 
    % I didn't add the hd here, it is unecessary.
    % Other indicator functions can be used. A search in the literature must be
    % done to find a more suitable one.
    Ig = reshape(log(sqrt(sum(conj(g).*g,1))),Ngrid,Ngrid);


end