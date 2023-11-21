function ipass = test_frechet_ders()
    ipass = 1;
    n  = 300;
    a = 1.1;  
    b = 1.8;  
    src_info = geometries.ellipse(a,b,n);
    narm = 3; rad = 1.1; amp = 0.2;
    starcoefs = zeros(2*narm+1,1);
    starcoefs(1) = rad;
    if narm > 0
        starcoefs(narm+1) = amp;
    end

    n  = max(300,100*narm);
    src_info = geometries.starn(starcoefs,narm,n);
    

    plot(src_info.xs,src_info.ys,'b.');
    drawnow;

    src0 = [0.01;-0.12];
    opts = [];
    opts.test_analytic = true;
    opts.src_in = src0;

    % set target locations
    %receptors
    r_tgt = 10;
    n_tgt = 88;
    t_tgt = 0:2*pi/n_tgt:2*pi-2*pi/n_tgt;

    % Incident directions
    n_dir = 100;
    t_dir = 0:2*pi/n_dir:2*pi-2*pi/n_dir;

    [t_tgt_grid,t_dir_grid] = meshgrid(t_tgt,t_dir);
    t_tgt_grid = t_tgt_grid(:);
    t_dir_grid = t_dir_grid(:);
    xtgt = r_tgt*cos(t_tgt_grid);
    ytgt = r_tgt*sin(t_tgt_grid);
    tgt   = [ xtgt'; ytgt'];

    % use a random subset of sensor and target location
    % so that they are no longer tensor product
    ifac = 1;
    iind_use = randperm(n_dir*n_tgt,ceil(ifac*n_dir*n_tgt));
    iind_use = 1:(n_dir*n_tgt);
    sensor_info = [];
    sensor_info.tgt = tgt(:,iind_use);
    sensor_info.t_dir = t_dir_grid(iind_use);


    [t_dir_uni,~,idir] = unique(sensor_info.t_dir);
    [tgt_uni,~,itgt] = unique(sensor_info.tgt','rows');
    nt_uni = length(tgt_uni(:,1));
    induse = itgt + (idir-1)*nt_uni;


    kh = 2.2;

    % Test obstacle Frechet derivative for Dirichlet problem
    bc = [];
    bc.type = 'Dirichlet';
    bc.invtype = 'o';

    [mats,~] = rla.get_fw_mats(kh,src_info,bc,sensor_info,opts);
    fields = rla.compute_fields(kh,src_info,mats,sensor_info,bc,opts);

    nh = 4;
    opts.ncoeff_boundary = nh;
    opts.ncoeff_impedance = nh;
    frechet_mats = rla.get_frechet_ders(kh,mats,src_info,sensor_info,fields,bc,opts);

    rng(1)
    hcoefs = 0.1*rand(1,2*nh+1);
    uder = frechet_mats.bdry*hcoefs(:);

    errs = zeros(2,1);
    for ig=1:2
        dh = 10^(-ig);
        hcoefs_use = dh*hcoefs;
        [src_out,~] = rla.update_geom(src_info,nh,hcoefs_use);
        [mats1,~] = rla.get_fw_mats(kh,src_out,bc,sensor_info,opts);
        fields1 = rla.compute_fields(kh,src_out,mats1,sensor_info,bc,opts);


        hcoefs_use2 = -dh*hcoefs;
        [src_out2,~] = rla.update_geom(src_info,nh,hcoefs_use2);
        [mats2,~] = rla.get_fw_mats(kh,src_out2,bc,sensor_info,opts);
        fields2 = rla.compute_fields(kh,src_out2,mats2,sensor_info,bc,opts);

        uder_est = (fields1.uscat_tgt(:) - fields2.uscat_tgt(:))/2/dh;
        errs(ig) = norm(uder-uder_est);
    end
    if(errs(2)>1e-4) 
        ipass = 0;
        fprintf('failed Dirichlet test in frechet_ders test: %d\n',errs(2));
    end
    
    

    % Test obstacle Frechet derivative for Neumann problem
    bc = [];
    bc.type = 'Neumann';
    bc.invtype = 'o';

    [mats,~] = rla.get_fw_mats(kh,src_info,bc,sensor_info,opts);
    fields = rla.compute_fields(kh,src_info,mats,sensor_info,bc,opts);


    frechet_mats = rla.get_frechet_ders(kh,mats,src_info,sensor_info,fields,bc,opts);
    uder = frechet_mats.bdry*hcoefs(:);


    errs = zeros(2,1);
    for ig=1:2
        dh = 10^(-ig);
        hcoefs_use = dh*hcoefs;
        [src_out,~] = rla.update_geom(src_info,nh,hcoefs_use);
        [mats1,~] = rla.get_fw_mats(kh,src_out,bc,sensor_info,opts);
        fields1 = rla.compute_fields(kh,src_out,mats1,sensor_info,bc,opts);


        hcoefs_use2 = -dh*hcoefs;
        [src_out2,~] = rla.update_geom(src_info,nh,hcoefs_use2);
        [mats2,~] = rla.get_fw_mats(kh,src_out2,bc,sensor_info,opts);
        fields2 = rla.compute_fields(kh,src_out2,mats2,sensor_info,bc,opts);

        uder_est = (fields1.uscat_tgt(:) - fields2.uscat_tgt(:))/2/dh;
        errs(ig) = norm(uder-uder_est);
    end
    if(errs(2)>1e-4) 
        ipass = 0;
        fprintf('failed Neumann test in frechet_ders test: %d\n',errs(2));
    end
    


    % Test obstacle frechet derivative of impedance problem
    t = 0:2*pi/n:2*pi*(1.0-1.0/n);
    src_info.lambda = sin(2*t).';
    src_info.lambda = (sin(2*t) + 0.1*1i*(1+cos(2*t))).';
    bc = [];
    bc.type = 'Impedance';
    bc.invtype = 'oi';

    [mats,~] = rla.get_fw_mats(kh,src_info,bc,sensor_info,opts);
    fields = rla.compute_fields(kh,src_info,mats,sensor_info,bc,opts);


    frechet_mats = rla.get_frechet_ders(kh,mats,src_info,sensor_info,fields,bc,opts);
    uder = frechet_mats.bdry*hcoefs(:);


    errs = zeros(2,1);
    for ig=1:2
        dh = 10^(-ig);
        hcoefs_use = dh*hcoefs;
        [src_out,~] = rla.update_geom(src_info,nh,hcoefs_use);
        [mats1,~] = rla.get_fw_mats(kh,src_out,bc,sensor_info,opts);
        fields1 = rla.compute_fields(kh,src_out,mats1,sensor_info,bc,opts);


        hcoefs_use2 = -dh*hcoefs;
        [src_out2,~] = rla.update_geom(src_info,nh,hcoefs_use2);
        [mats2,~] = rla.get_fw_mats(kh,src_out2,bc,sensor_info,opts);
        fields2 = rla.compute_fields(kh,src_out2,mats2,sensor_info,bc,opts);

        uder_est = (fields1.uscat_tgt(:) - fields2.uscat_tgt(:))/2/dh;
        errs(ig) = norm(uder-uder_est);
    end
    if(errs(2)>1e-4) 
        ipass = 0;
        fprintf('failed Impedance test obstacle in frechet_ders test: %d\n',errs(2));
    end
    
    % Test impedance frechet derivative of impedance problem

    hcoefs =  (0.1*rand(1,2*nh+1) + 0.1*1i*rand(1,2*nh+1))/sqrt(2);
    
    uder = frechet_mats.impedance*hcoefs(:);
    errs = zeros(2,1);
    for ig=1:2
        dh = 10^(-ig);
        hcoefs_use = dh*hcoefs;
        hcoefs_use = hcoefs_use(:);
        src_out = src_info;
        h_upd = (cos(t.'*(0:nh))*hcoefs_use(1:(nh+1)) + sin(t.'*(1:nh))*hcoefs_use((nh+2):end)).';
        src_out.lambda = src_info.lambda + h_upd.';
        [mats1,~] = rla.get_fw_mats(kh,src_out,bc,sensor_info,opts);
        fields1 = rla.compute_fields(kh,src_out,mats1,sensor_info,bc,opts);

        src_out2 = src_info;
        src_out2.lambda = src_info.lambda - h_upd.';
        [mats2,~] = rla.get_fw_mats(kh,src_out2,bc,sensor_info,opts);
        fields2 = rla.compute_fields(kh,src_out2,mats2,sensor_info,bc,opts);

        uder_est = (fields1.uscat_tgt(:) - fields2.uscat_tgt(:))/2/dh;
        errs(ig) = norm(uder-uder_est);
    end
    if(errs(2)>1e-4)
        ipass = 0;
        fprintf('failed Impedance test impedance in frechet_ders test: %d\n',errs(2));
    end
    
    % Test obstacle frechet derivative of impedance problem
    % impedance_type = constkappa
    t = 0:2*pi/n:2*pi*(1.0-1.0/n);
    rng(123)
    nh = 4;
    hcoefs = 0.1*rand(1,2*nh+1);
    
    src_info.lamcfs = [1.1+0.5*1i; 0.4-0.3*1i];
    src_info.lambda = src_info.lamcfs(1) + src_info.lamcfs(2)*src_info.H;
    bc = [];
    bc.type = 'Impedance';
    bc.invtype = 'oi';
    opts.impedance_type = 'constkappa';
    opts.kh = kh;

    [mats,~] = rla.get_fw_mats(kh,src_info,bc,sensor_info,opts);
    fields = rla.compute_fields(kh,src_info,mats,sensor_info,bc,opts);

    frechet_mats = rla.get_frechet_ders(kh,mats,src_info,sensor_info,fields,bc,opts);
    uder = frechet_mats.bdry*hcoefs(:);

    errs = zeros(2,1);
    for ig=1:2
        dh = 10^(-ig);
        hcoefs_use = dh*hcoefs;
        [src_out,~] = rla.update_geom(src_info,nh,hcoefs_use,opts);
        [mats1,~] = rla.get_fw_mats(kh,src_out,bc,sensor_info,opts);
        fields1 = rla.compute_fields(kh,src_out,mats1,sensor_info,bc,opts);


        hcoefs_use2 = -dh*hcoefs;
        [src_out2,~] = rla.update_geom(src_info,nh,hcoefs_use2,opts);
        [mats2,~] = rla.get_fw_mats(kh,src_out2,bc,sensor_info,opts);
        fields2 = rla.compute_fields(kh,src_out2,mats2,sensor_info,bc,opts);

        uder_est = (fields1.uscat_tgt(:) - fields2.uscat_tgt(:))/2/dh;
        errs(ig) = norm(uder-uder_est);
    end
    if(errs(2)>1e-4) 
        ipass = 0;
        fprintf('failed Impedance test obstacle in frechet_ders test (constkappa): %d\n',errs(2));
    end
    
    % Test impedance frechet derivative of impedance problem
    % impedance_type = const_kappa

    hcoefs =  (0.1*rand(2,1) + 0.1*1i*rand(2,1));
    
    uder = frechet_mats.impedance*hcoefs(:);
    errs = zeros(2,1);
    for ig=1:2
        dh = 10^(-ig);
        hcoefs_use = dh*hcoefs;
        hcoefs_use = hcoefs_use(:);
        src_out = src_info;
        src_out.lamcfs = src_out.lamcfs(:) + hcoefs_use(:);
        src_out.lambda = src_out.lamcfs(1) + src_out.lamcfs(2)*src_out.H;
        [mats1,~] = rla.get_fw_mats(kh,src_out,bc,sensor_info,opts);
        fields1 = rla.compute_fields(kh,src_out,mats1,sensor_info,bc,opts);

        src_out2 = src_info;
        src_out2.lamcfs = src_out2.lamcfs(:) - hcoefs_use(:);
        src_out2.lambda = src_out2.lamcfs(1) + src_out2.lamcfs(2)*src_out2.H;
        [mats2,~] = rla.get_fw_mats(kh,src_out2,bc,sensor_info,opts);
        fields2 = rla.compute_fields(kh,src_out2,mats2,sensor_info,bc,opts);

        uder_est = (fields1.uscat_tgt(:) - fields2.uscat_tgt(:))/2/dh;
        errs(ig) = norm(uder-uder_est);
    end
    if(errs(2)>1e-4)
        ipass = 0;
        fprintf('failed Impedance test impedance in frechet_ders test (constkappa): %d\n',errs(2));
    end
    
    % Test obstacle frechet derivative of impedance problem
    % impedance_type = antbar2
    t = 0:2*pi/n:2*pi*(1.0-1.0/n);
    rng(12)
    nh = 4;
    hcoefs = 0.1*rand(1,2*nh+1);

    opts.impedance_type = 'antbar2';    
    opts.kh = kh;

    src_info.lamcfs = [1.3;0.7];
    ckcfs = constkappa_models_convert(src_info.lamcfs,opts.impedance_type,kh);
    src_info.lambda = ckcfs(1) + ckcfs(2)*src_info.H;
    bc = [];
    bc.type = 'Impedance';
    bc.invtype = 'oi';


    [mats,~] = rla.get_fw_mats(kh,src_info,bc,sensor_info,opts);
    fields = rla.compute_fields(kh,src_info,mats,sensor_info,bc,opts);

    frechet_mats = rla.get_frechet_ders(kh,mats,src_info,sensor_info,fields,bc,opts);
    uder = frechet_mats.bdry*hcoefs(:);

    errs = zeros(2,1);
    for ig=1:2
        dh = 10^(-ig);
        hcoefs_use = dh*hcoefs;
        [src_out,~] = rla.update_geom(src_info,nh,hcoefs_use,opts);
        [mats1,~] = rla.get_fw_mats(kh,src_out,bc,sensor_info,opts);
        fields1 = rla.compute_fields(kh,src_out,mats1,sensor_info,bc,opts);


        hcoefs_use2 = -dh*hcoefs;
        [src_out2,~] = rla.update_geom(src_info,nh,hcoefs_use2,opts);
        [mats2,~] = rla.get_fw_mats(kh,src_out2,bc,sensor_info,opts);
        fields2 = rla.compute_fields(kh,src_out2,mats2,sensor_info,bc,opts);

        uder_est = (fields1.uscat_tgt(:) - fields2.uscat_tgt(:))/2/dh;
        errs(ig) = norm(uder-uder_est);
    end
    if(errs(2)>2*1e-4) 
        ipass = 0;
        fprintf('failed Impedance test obstacle in frechet_ders test (antbar2): %d\n',errs(2));
    end
    
    % Test impedance frechet derivative of impedance problem
    % impedance_type = antbar2

    hcoefs =  0.1*rand(2,1);
    
    uder = frechet_mats.impedance*hcoefs(:);
    errs = zeros(2,1);
    for ig=1:2
        dh = 10^(-ig);
        hcoefs_use = dh*hcoefs;
        hcoefs_use = hcoefs_use(:);
        src_out = src_info;
        src_out.lamcfs = src_out.lamcfs(:) + hcoefs_use(:);
        ckcfs = constkappa_models_convert(src_out.lamcfs,opts.impedance_type,kh);
        src_out.lambda = ckcfs(1) + ckcfs(2)*src_out.H;
        [mats1,~] = rla.get_fw_mats(kh,src_out,bc,sensor_info,opts);
        fields1 = rla.compute_fields(kh,src_out,mats1,sensor_info,bc,opts);

        src_out2 = src_info;
        src_out2.lamcfs = src_out2.lamcfs(:) - hcoefs_use(:);
        ckcfs = constkappa_models_convert(src_out2.lamcfs,opts.impedance_type,kh);
        src_out2.lambda = ckcfs(1) + ckcfs(2)*src_out2.H;
        [mats2,~] = rla.get_fw_mats(kh,src_out2,bc,sensor_info,opts);
        fields2 = rla.compute_fields(kh,src_out2,mats2,sensor_info,bc,opts);

        uder_est = (fields1.uscat_tgt(:) - fields2.uscat_tgt(:))/2/dh;
        errs(ig) = norm(uder-uder_est);
    end
    if(errs(2)>1e-4)
        ipass = 0;
        fprintf('failed Impedance test impedance in frechet_ders test (antbar2): %d\n',errs(2));
    end
    
    % Test obstacle frechet derivative of impedance problem
    % impedance_type = antbar3
    t = 0:2*pi/n:2*pi*(1.0-1.0/n);
    rng(12)
    nh = 4;
    hcoefs = 0.1*rand(1,2*nh+1);

    opts.impedance_type = 'antbar3';    
    opts.kh = kh;
    src_info.lamcfs = [1.3;0.7;0.4];
    ckcfs = constkappa_models_convert(src_info.lamcfs,opts.impedance_type,kh);
    src_info.lambda = ckcfs(1) + ckcfs(2)*src_info.H;
    bc = [];
    bc.type = 'Impedance';
    bc.invtype = 'oi';


    [mats,~] = rla.get_fw_mats(kh,src_info,bc,sensor_info,opts);
    fields = rla.compute_fields(kh,src_info,mats,sensor_info,bc,opts);

    frechet_mats = rla.get_frechet_ders(kh,mats,src_info,sensor_info,fields,bc,opts);
    uder = frechet_mats.bdry*hcoefs(:);

    errs = zeros(2,1);
    for ig=1:2
        dh = 10^(-ig);
        hcoefs_use = dh*hcoefs;
        [src_out,~] = rla.update_geom(src_info,nh,hcoefs_use,opts);
        [mats1,~] = rla.get_fw_mats(kh,src_out,bc,sensor_info,opts);
        fields1 = rla.compute_fields(kh,src_out,mats1,sensor_info,bc,opts);


        hcoefs_use2 = -dh*hcoefs;
        [src_out2,~] = rla.update_geom(src_info,nh,hcoefs_use2,opts);
        [mats2,~] = rla.get_fw_mats(kh,src_out2,bc,sensor_info,opts);
        fields2 = rla.compute_fields(kh,src_out2,mats2,sensor_info,bc,opts);

        uder_est = (fields1.uscat_tgt(:) - fields2.uscat_tgt(:))/2/dh;
        errs(ig) = norm(uder-uder_est);
    end
    if(errs(2)>2*1e-4) 
        ipass = 0;
        fprintf('failed Impedance test obstacle in frechet_ders test (antbar3): %d\n',errs(2));
    end
    
    % Test impedance frechet derivative of impedance problem
    % impedance_type = antbar3

    hcoefs =  0.1*rand(3,1);
    
    uder = frechet_mats.impedance*hcoefs(:);
    errs = zeros(2,1);
    for ig=1:2
        dh = 10^(-ig);
        hcoefs_use = dh*hcoefs;
        hcoefs_use = hcoefs_use(:);
        src_out = src_info;
        src_out.lamcfs = src_out.lamcfs(:) + hcoefs_use(:);
        ckcfs = constkappa_models_convert(src_out.lamcfs,opts.impedance_type,kh);
        src_out.lambda = ckcfs(1) + ckcfs(2)*src_out.H;
        [mats1,~] = rla.get_fw_mats(kh,src_out,bc,sensor_info,opts);
        fields1 = rla.compute_fields(kh,src_out,mats1,sensor_info,bc,opts);

        src_out2 = src_info;
        src_out2.lamcfs = src_out2.lamcfs(:) - hcoefs_use(:);
        ckcfs = constkappa_models_convert(src_out2.lamcfs,opts.impedance_type,kh);
        src_out2.lambda = ckcfs(1) + ckcfs(2)*src_out2.H;
        [mats2,~] = rla.get_fw_mats(kh,src_out2,bc,sensor_info,opts);
        fields2 = rla.compute_fields(kh,src_out2,mats2,sensor_info,bc,opts);

        uder_est = (fields1.uscat_tgt(:) - fields2.uscat_tgt(:))/2/dh;
        errs(ig) = norm(uder-uder_est);
    end
    if(errs(2)>1e-4)
        ipass = 0;
        fprintf('failed Impedance test impedance in frechet_ders test (antbar3): %d\n',errs(2));
    end
    
    % Test impedance frechet derivative of impedance problem
    % impedance_type = antbarphys
    
    opts.impedance_type = 'antbarphys';

    hcoefs =  0.1*rand(3,1);
    src_info.lamcfs = [1.3;0.7;0.4];

    ckcfs = constkappa_models_convert(src_info.lamcfs,opts.impedance_type,kh);
    src_info.lambda = ckcfs(1) + ckcfs(2)*src_info.H;
    bc = [];
    bc.type = 'Impedance';
    bc.invtype = 'oi';
    
    [mats,~] = rla.get_fw_mats(kh,src_info,bc,sensor_info,opts);
    fields = rla.compute_fields(kh,src_info,mats,sensor_info,bc,opts);
    frechet_mats = rla.get_frechet_ders(kh,mats,src_info,sensor_info,fields,bc,opts);

    uder = frechet_mats.impedance*hcoefs(:);
    errs = zeros(2,1);

    for ig=1:2
        dh = 10^(-ig);
        hcoefs_use = dh*hcoefs;
        hcoefs_use = hcoefs_use(:);
        src_out = src_info;
        src_out.lamcfs = src_out.lamcfs(:) + hcoefs_use(:);
        ckcfs = constkappa_models_convert(src_out.lamcfs,opts.impedance_type,kh);
        src_out.lambda = ckcfs(1) + ckcfs(2)*src_out.H;
        [mats1,~] = rla.get_fw_mats(kh,src_out,bc,sensor_info,opts);
        fields1 = rla.compute_fields(kh,src_out,mats1,sensor_info,bc,opts);

        src_out2 = src_info;
        src_out2.lamcfs = src_out2.lamcfs(:) - hcoefs_use(:);
        ckcfs = constkappa_models_convert(src_out2.lamcfs,opts.impedance_type,kh);
        src_out2.lambda = ckcfs(1) + ckcfs(2)*src_out2.H;
        [mats2,~] = rla.get_fw_mats(kh,src_out2,bc,sensor_info,opts);
        fields2 = rla.compute_fields(kh,src_out2,mats2,sensor_info,bc,opts);

        uder_est = (fields1.uscat_tgt(:) - fields2.uscat_tgt(:))/2/dh;
        errs(ig) = norm(uder-uder_est);
    end
    if(errs(2)>1e-4)
        errs
        ipass = 0;
        fprintf('failed Impedance test impedance in frechet_ders test (antbarphys): %d\n',errs(2));
    end
    
end

