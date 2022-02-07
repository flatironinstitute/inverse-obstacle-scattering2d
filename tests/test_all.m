ntest  = 3;
ipass = zeros(ntest,1);
ipass(1) = test_update_geom();
ipass(2) = test_fw_mats();
ipass(3) = test_frechet_ders();
fprintf('Succesfully completed %d out of %d tests in the ios2d testing suite\n',sum(ipass),ntest);