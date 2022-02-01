        subroutine dumb_congrad(matvec,a,ww,y,n,numit,x,
     1      errs,niter,w,forms)
        implicit real *8 (a-h,o-z)
        save
        dimension a(n,n),x(1),y(1),errs(1),w(1),forms(1),ww(1)
c
c       This subroutine applies a somewhat strange form of
c       Conjugate Gradient algorithm to the system of linear
c       equations. The strangeness of the scheme is the fact that 
c       it it is performed with complete reorthogonilization
c
c                 Ax=y
c
c            
c                    Input parameters:
c
c  matvec - the subroutine applying the matrix a to arbitrary vectors
c        The calling sequence of matvec must be
c
c                 matvec(a,n,x,y)
c
c  a - the matrix of the system
c  y - the right-hand side
c  n - the dimensionality of the system
c  eps - the accuracy to which the problem will be solved 
c       (in the sense of residual)
c  numit - the maximum number of iterations to be performed
c
c                    Output parameters:
c
c  x - the solution
c  errs - the array of errors produced by consecutive iterations
c      (niter of them things). Please note that 
c
c               errs(i)= ||A(X_i) - Y||,
c      and DO NOT decrease monotonically
c
c  niter - the number of iterations actually performed, until one
c      of three conditions has been achieved:
c    1. numit iterations have been performed, or
c    2. the quadratic form (see below) stopped decreasing
c  forms - forms(i) is the value of the quadratic form 
c
c         1/2 (AX,X) - (X,Y)
c
c      obtained diring the i-th iteration. This paramer is supplied
c      entirely for the user's ediification. 
c
c                    Work arrays:
c
c  w - must be at least 2*numit*n+n+30 real *8 elements long
c
c
c        . . . allocate memory
c
c
        ixs=1
        lxs=numit*n+10
c
        iaxs=ixs+lxs
        laxs=numit*n+10
c
        iyr=iaxs+laxs
        lyr=n+2
c
        call dumb_congrad0(matvec,a,ww,y,n,numit,w(ixs),w(iaxs),x,
     1      w(iyr),errs,niter,forms)
c
        return
        end
c 
c 
c 
c 
c 
        subroutine dumb_congrad0(matvec,a,ww,y,n,numit,xs,axs,x,yr,
     1      errs,niter,forms)
        implicit real *8 (a-h,o-z)
        save
        dimension a(n,n),xs(n,1),axs(n,1),x(1),yr(1),y(1),
     1      errs(1),forms(1),ww(1)
c
c        initialize the process
c
        call dumb_conres_copy(y,xs,n)
        call matvec(a,n,y,axs,ww)
c
        call dumb_conres_scap(xs,axs,n,d)
        d=1/sqrt(d)
        do 1200 i=1,n
c
        axs(i,1)=axs(i,1)*d
        xs(i,1)=xs(i,1)*d
 1200 continue
c
        call dumb_conres_scap(xs,y,n,d)
c
        do 1400 i=1,n
c
        x(i)=xs(i,1)*d
        yr(i)=y(i)-axs(i,1)*d
 1400 continue
c
c        conduct the iterations
c
        n_numit=numit-1
        if(n .lt. numit) n_numit=n-1
        
        do 3000 i=1,n_numit
c
c       construct the new direction
c
        niter=i

cccc        call prinf('i=*',i,1)
c
        call dumb_conres_copy(yr,xs(1,i+1),n)
        call matvec(a,n,yr,axs(1,i+1),ww )
c
c        orthogonalize the new direction to all preceding ones
c
        do 1900 ijk=1,2
c
        do 1800 j=1,i

ccc        ii=i-1
ccc        if(ii .lt. 1) ii=1
ccc        do 1800 j=ii,i
c
        call dumb_conres_scap(axs(1,j),xs(1,i+1),n,d)
c
        do 1600 jj=1,n
c
        axs(jj,i+1)=axs(jj,i+1)-axs(jj,j)*d
        xs(jj,i+1)=xs(jj,i+1)-xs(jj,j)*d
 1600 continue
 1800 continue
c
 1900 continue
c
c       . . . normalize the new vector
c
        call dumb_conres_scap(xs(1,i+1),axs(1,i+1),n,d)
        d=1/sqrt(d)

        do 2200 j=1,n
c
        xs(j,i+1)=xs(j,i+1)*d
        axs(j,i+1)=axs(j,i+1)*d
 2200 continue
c
        call matvec(a,n,xs(1,i+1),axs(1,i+1),ww )
c
c       update the solution and the residual
c
        call dumb_conres_scap(xs(1,i+1),yr,n,d)
c
        dd=0
        do 2400 j=1,n
c
        yr(j)=yr(j)-axs(j,i+1)*d
        x(j)=x(j)+xs(j,i+1)*d
c
        dd=dd+yr(j)**2
 2400 continue
c
c       calculate the quadratic form
c
        call dumb_conres_scap(x,yr,n,form)
        call dumb_conres_scap(x,y,n,d)
c
        form=(form-d)/2        
        forms(i)=form
c
        errs(i)=sqrt(dd)
c
        if( (i .gt. 1) .and. (forms(i) .gt. forms(i-1)) ) then
            niter=i
cccc            call prin2('and forms=*',forms,niter)
            return
        endif 
c
 3000 continue
c
 3100 continue
c
        if(2 .ne. 3) return
c
        do 3400 i=1,niter
        do 3200 j=i,niter
c
        call dumb_conres_scap(xs(1,i),axs(1,j),n,d)
c
 3250 format('   i= ',i4,', j= ',i4,', ', '(xs(i),axs(j))= ',e11.5)
c
        write(6,3250) i,j,d
        write(13,3250) i,j,d
c
 3200 continue
 3400 continue
c
        return
        end
c 
c 
c 
c 
c 
        subroutine dumb_conres(matvec,a,ww,y,n,eps,numit,x,errs,niter,w)
        implicit real *8 (a-h,o-z)
        save
        dimension a(n,n),x(1),y(1),errs(1),w(1),ww(1)
c
c       This subroutine applies a somewhat strange form of
c       Minimum Residual algorithm to the system of linear
c       equations. The strangeness of the scheme is the fact that 
c       it it is performed with complete reorthogonilization
c
c                 Ax=y
c
c            
c                    Input parameters:
c
c  matvec - the subroutine applying the matrix a to arbitrary vectors
c        The calling sequence of matvec must be
c
c                 matvec(a,n,x,y)
c
c  a - the matrix of the system
c  y - the right-hand side
c  n - the dimensionality of the system
c  eps - the accuracy to which the problem will be solved 
c       (in the sense of residual)
c  numit - the maximum number of terations to be performed
c
c                    Output parameters:
c
c  x - the solution
c  errs - the array of errors produced by consecutive iterations
c      (niter of them things)
c  niter - the number of iterations actually performed, until one
c      of three conditions has been achieved:
c    1. relative error eps has been achieved, or 
c    2. numit iterations have been performed, or
c    3. the error stopped decreasing
c
c                    Work arrays:
c
c  w - must be at least 2*numit*n+n+30 real *8 elements long
c
c
c        . . . allocate memory
c
        ixs=1
        lxs=numit*n+10
c
        iaxs=ixs+lxs
        laxs=numit*n+10
c
        iyr=iaxs+laxs
        lyr=n+2
c
        call dumb_conres0(matvec,a,ww,y,n,numit,eps,
     1      w(ixs),w(iaxs),x,w(iyr),errs,niter)
c
        return
        end
c 
c 
c 
c 
c 
        subroutine dumb_conres0(matvec,a,ww,y,n,numit,eps,
     1      xs,axs,x,yr,errs,niter)
        implicit real *8 (a-h,o-z)
        save
        dimension a(n,n),xs(n,1),axs(n,1),x(1),yr(1),y(1),
     1      errs(1),ww(1)
c
c        initialize the process
c
        nrec=5
        call dumb_conres_copy(y,xs,n)
        call matvec(a,n,y,axs,ww)
c
        call dumb_conres_scap(axs,axs,n,d)
        d=1/sqrt(d)
        do 1200 i=1,n
c
        axs(i,1)=axs(i,1)*d
        xs(i,1)=xs(i,1)*d
 1200 continue
c
        call dumb_conres_scap(axs,y,n,d)
        do 1400 i=1,n
c
        x(i)=xs(i,1)*d
        yr(i)=y(i)-axs(i,1)*d
 1400 continue
c
c        conduct the iterations
c
        irec=0
        do 3000 i=1,numit-1
c
c       construct the new direction
c
        irec=irec+1
        call dumb_conres_copy(yr,xs(1,i+1),n)
        call matvec(a,n,yr,axs(1,i+1),ww )
c
c        orthogonalize the new direction to all preceding ones
c
        do 1900 ijk=1,2
c
        do 1800 j=1,i

ccc        ii=i-1
ccc        if(ii .lt. 1) ii=1
ccc        do 1800 j=ii,i
c
        call dumb_conres_scap(axs(1,j),axs(1,i+1),n,d)


cccc        call prin2('and d=*',d,1)
c
        do 1600 jj=1,n
c
        axs(jj,i+1)=axs(jj,i+1)-axs(jj,j)*d
        xs(jj,i+1)=xs(jj,i+1)-xs(jj,j)*d
 1600 continue
c
        if((ijk .eq. 1) .and. (irec .eq. nrec) ) then
            irec=0
            call matvec(a,n,xs(1,i+1),axs(1,i+1),ww )
        endif
c
 1800 continue
c
 1900 continue
c
c       . . . normalize the new vector
c
        call dumb_conres_scap(axs(1,i+1),axs(1,i+1),n,d)
        d=1/sqrt(d)

        do 2200 j=1,n
c
        xs(j,i+1)=xs(j,i+1)*d
        axs(j,i+1)=axs(j,i+1)*d
 2200 continue
c
c       update the solution and the residual
c
        call dumb_conres_scap(axs(1,i+1),yr,n,d)
c
        dd=0
        do 2400 j=1,n
c
        yr(j)=yr(j)-axs(j,i+1)*d
        x(j)=x(j)+xs(j,i+1)*d
c
        dd=dd+yr(j)**2
 2400 continue
c
        errs(i)=sqrt(dd)
c
        if(i .eq. 1) goto 3000
        if(errs(i) .ge. errs(i-1)) then
            niter=i
            return
        endif 
c
c
        if(errs(i) .lt. errs(1)*eps) then
            niter=i
            return
        endif 
c
 3000 continue
c
        niter=numit
c
        if(2 .ne. 3) return

        do 3400 i=1,numit
        do 3200 j=1,numit
c
        call dumb_conres_scap(axs(1,i),axs(1,j),n,d)

        xs(j,i)=d

 3200 continue
 3400 continue

        call prin2('and prods=*',xs,numit**2)



        return
        end
c 
c 
c 
c 
c 
        subroutine dumb_conres_scap(x,y,n,prod)
        implicit real *8 (a-h,o-z)
        save
        dimension x(1),y(1)
c 
        prod=0
        do 1200 i=1,n
        prod=prod+x(i)*(y(i))
 1200 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine dumb_conres_copy(x,y,n)
        implicit real *8 (a-h,o-z)
        save
        real *8 x(1),y(1)
c 
        do 2400 i=1,n
        y(i)=x(i)
 2400 continue
        return
        end  
