function [] = post_process(inv_data,fname)
   inv_tmp = cell2mat(inv_data);
   iter_count = vertcat(inv_tmp.iter_count);
   exit_criterion = vertcat(inv_tmp.exit_criterion);
   res_all = vertcat(inv_tmp.res_all);
   res_opt = vertcat(inv_tmp.res_opt);
   kh = vertcat(inv_tmp.kh);
   iter_vec = zeros(sum(iter_count),1);
   
   maxiter = max(iter_count);
   khvec2 = zeros(sum(iter_count),1);
   dkh = (kh(2)-kh(1))/(maxiter+5);
   istart = 1;
   for i=1:length(iter_count)
       iter_vec(istart:(istart+iter_count(i)-1)) = 1:iter_count(i).';
       khvec2(istart:(istart+iter_count(i)-1)) = kh(i) + (0:(iter_count(i)-1).')*dkh;
       istart = istart + iter_count(i);
   end
   
      
      
   
   kmin = min(kh);
   kmax = max(kh);
   
   
   
   figure;
   semilogy(kh,res_opt,'k.','MarkerSize',10)
   title('Residue of optimal shape');
   rmin = min(res_opt);
   rmax = max(res_opt);
   ylim([rmin/10,rmax*10]);
   xlim([kmin-0.25,kmax+0.25]);
   xlabel('Frequency');
   drawnow;
   
   
   figure;
   semilogy(khvec2,res_all,'k.','MarkerSize',10)
   title('Residue');
   rmin = min(res_all);
   rmax = max(res_all);
   ylim([rmin/10,rmax*10]);
   xlim([kmin-0.25,kmax+0.25]);
   xlabel('Frequency');
   drawnow;
   
   
   figure;
   plot(kh,iter_count,'k.','MarkerSize',10)
   miniter = min(iter_count);
   ylim([miniter-1,maxiter+1]);
   title('Iteration count');
   xlim([kmin-0.25,kmax+0.25]);
   xlabel('Frequency');
   drawnow;
   
   figure;
   plot(kh,exit_criterion,'k.','MarkerSize',10)
   ylim([-2,5])
   title('Exit criterion');
   xlim([kmin-0.25,kmax+0.25]);
   yticks([-1,1,2,3]);
   yticklabels({'Self-intersection','small res','small upd','max iter'})
   xlabel('Frequency')
   drawnow;
   
   n = length(res_all);
   F(n) = struct('cdata',[],'colormap',[]);
   khvec = repelem(kh,iter_count);
   
   
   stmp = vertcat(inv_tmp.src_info_all);
   xsall = []; ysall = [];
   for i = 1:length(stmp)
       xsall = [xsall; stmp{i}.xs(:)];
       ysall = [ysall; stmp{i}.ys(:)];
   end
   
   
   xmin = min(xsall);
   xmax = max(xsall);
   ymin = min(ysall);
   ymax = max(ysall);
   
   if(nargin == 2)
       S = load(fname);
       src_ex = S.src_info;
       xmin0 = min(src_ex.xs);
       xmax0 = max(src_ex.xs);
       ymin0 = min(src_ex.ys);
       ymax0 = max(src_ex.ys);
       xmin = min(xmin0,xmin);
       ymin = min(ymin0,ymin);
       xmax = max(xmax0,xmax);
       ymax = max(ymax0,ymax);
   end
   xdiff = xmax-xmin;
   ydiff = ymax-ymin;
   xmid = 0.5*(xmax+xmin);
   ymid = 0.5*(ymax+ymin);
   xmin = xmid - 0.6*xdiff;
   xmax = xmid + 0.6*xdiff;
   ymin = ymid - 0.6*ydiff;
   ymax = ymid + 0.6*ydiff;
   figure
   figure('Renderer','zbuffer')
   if(nargin == 1)
      plot(stmp{i}.xs,stmp{i}.ys,'k.')
   else
      plot(stmp{i}.xs,stmp{i}.ys,'k.',src_ex.xs,src_ex.ys,'b.')
   end
   title(['kh=' num2str(khvec(1)) ]);
   axis equal
   xlim([xmin,xmax])
   ylim([ymin,ymax])
   drawnow;
   set(gca,'NextPlot','replacechildren');
   for i=1:n
       if(nargin == 1)
          plot(stmp{i}.xs,stmp{i}.ys,'k.')
       else
          plot(stmp{i}.xs,stmp{i}.ys,'k.',src_ex.xs,src_ex.ys,'b.')
       end
       
       title(['kh=' num2str(khvec(i)) '  iteration=' num2str(iter_vec(i))]);
       axis equal
       xlim([xmin,xmax])
       ylim([ymin,ymax])
       drawnow;
       F(i) = getframe(gcf);
   end
   close(gcf)
   implay(F,1);
     
   
end