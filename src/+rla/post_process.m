function [] = post_process(inv_data,fname)
   inv_tmp = cell2mat(inv_data);
   iter_count = vertcat(inv_tmp.iter_count);
   exit_criterion = vertcat(inv_tmp.exit_criterion);
   res_all = vertcat(inv_tmp.res_all);
   res_opt = vertcat(inv_tmp.res_opt);
   kh = vertcat(inv_tmp.kh);
   
   kmin = min(kh);
   kmax = max(kh);
   
   
   maxiter = max(iter_count);
   figure;
   semilogy(kh,res_opt,'k.','MarkerSize',10)
   title('Residue of optimal shape');
   rmin = min(res_opt);
   rmax = max(res_opt);
   ylim([rmin/10,rmax*10]);
   xlim([kmin-0.25,kmax+0.25]);
   xlabel('Frequency');
   
   figure;
   plot(kh,iter_count,'k.','MarkerSize',10)
   miniter = min(iter_count);
   ylim([miniter-1,maxiter+1]);
   title('Iteration count');
   xlim([kmin-0.25,kmax+0.25]);
   xlabel('Frequency');
   
   figure;
   plot(kh,exit_criterion,'k.','MarkerSize',10)
   ylim([-2,5])
   title('Exit criterion');
   xlim([kmin-0.25,kmax+0.25]);
   yticks([-1,1,2,3]);
   yticklabels({'Self-intersection','small res','small upd','max iter'})
   xlabel('Frequency')
   
   
end