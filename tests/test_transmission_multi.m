
    addpath ../src/
    ipass = 1;
    n  = 300;
    a = 1.1;  
    b = 1.3;  
    coefs = [1.0 0 0 0.4 0 0 0];
    
    nc = 3;
    src_info = geometries.starn(coefs,nc,n);
    
    
    centers = [3,4,0; 3,0,4];
    scales = [1.0;1.0;1.0];
    src_info = geometries.get_multi_comp(src_info,centers,scales);
    

    plot(src_info.xs,src_info.ys,'b.');
    drawnow;

    src0 = [0.01+centers(1,1);-0.12+centers(2,1)];
    src_out = [10.1;1.1];
    opts = [];
    opts.test_analytic = true;
    opts.src_in = src0;
    opts.src_out = src_out;

    % set target locations
    %receptors
    r_tgt = 10;
    n_tgt = 100;
    t_tgt = 0:2*pi/n_tgt:2*pi-2*pi/n_tgt;
    x_t   = r_tgt * cos(t_tgt);
    y_t   = r_tgt * sin(t_tgt);    
    tgt   = [ x_t; y_t];


    sensor_info = [];
    sensor_info.tgt = tgt;


    kh = 2.2;

    bc = [];
    bc.type = 'Transmission';
    bc.transk = [2;1]; bc.transa = [1.3;1]; bc.transb = [1.5,1.1];
    [mats,erra] = rla.get_fw_mats(kh,src_info,bc,sensor_info,opts);
    if(erra>1e-11) 
        ipass = 0;
        fprintf('failed Transmission test in fw_mats test: %d\n',erra);
    end
