cc Copyright (C) 2010: Leslie Greengard and Mike O'Neil
cc Contact: greengard@cims.nyu.edu
cc Contact: oneil@cims.nyu.edu
c
c    This file contains matrix generators for the single layer potential,
c    the principal value of the double layer potential and the principal
c    value of the normal derivative of the single layer potential.
c
c
c    The bac suffix refers to the fact that we use Bradley Alpert
c    quadrature and that the matrices are complex.
c
c     The quadrature rule used is taken from:
c
c     B. Alpert,
c     Hybrid Gauss-Trapezoidal Quadrature Rules 
c     SIAM J. Sci. Comput. 20 (1999), pp. 1551-1584. 
c
c
c
        subroutine formmatbac(amat,norder,ns,srcinfo,h,gfun,dpars,
     1      zpars,ipars)
        implicit real *8 (a-h,o-z)
        real *8 srcinfo(6,ns),dpars(*)
        integer ipars(*)
        complex *16 zpars(*)
        real *8 tpts(100),txtra(100),coefs(100),srctmp(5),srcpts(5,100)
        complex *16 amat(ns,ns),ima,u
        integer *4 its(100),its2(100)
        dimension extranodes2(2), extraweights2(2)
        dimension extranodes3(4), extraweights3(4)
        dimension extranodes4(6), extraweights4(6)
        dimension extranodes8(14), extraweights8(14)
        dimension extranodes16(30), extraweights16(30)
        dimension extranodes(30), extraweights(30)
        external gfun
c
c       this routine builds the matrix which applies 
c       the logarithimacally singular kernels defined via 
c       gfun to a vector using alpert quadrature
c
c     input:
c       norder - the order of alpert quadrature to use, 0, 2, 4, 8, 16
c       srcinfo(5,ns) 
c           srcinfo(1:2,*) = xs,ys, srcinfo(3:4,*) = rnx,rny,
c           srcinfo(5,*) = dsdt
c           srcinfo(6,:) - curvature
c       xs,ys - the x,y coordinates of points on the curve
c       rnx,rny - the unit normal derivative at the points xs,ys
c       dsdt - derivative of parameterization
c       h - step size in parameterization variable
c       ns - length of xs and ys
c       gfun - routine that evaluates the kernel, must have a
c           calling sequence of the form
c           gfun(srcinfo,targinfo,dpars,zpars,ipars,uval)
c
c           where dpars - list of real parameters
c                 zpars - list of complex parameters
c                 ipars - list of integer parameters
c                 srcinfo - source info
c                 targinfo - target info
c                 uval - the complex *16 value of the kernel 
c                     at ((xtarg,ytarg),(xsrc,ysrc))
c
c     output:
c       amat - the ns x ns matrix that will apply gfun
c
c
      data extranodes2/
     1 -1.591549430918953D-01,
     1 1.591549430918953D-01/
      data extraweights2/
     1 0.5d0,
     1 0.5d0/

      data extranodes3/
     1 -1.150395811972836D-01,
     1 -9.365464527949632D-01,
     1 1.150395811972836D-01,
     1 9.365464527949632D-01/
      data extraweights3/
     1 3.913373788753340d-01,
     1 1.108662621124666d0,
     1 3.913373788753340d-01,
     1 1.108662621124666d0/

      data extranodes4/
     1 -1.023715124251890D+00,
     1 -2.935370741501914D-01,
     1 -2.379647284118974D-02, 
     1 2.379647284118974D-02, 
     1 2.935370741501914D-01,
     1 1.023715124251890D+00/ 
      data extraweights4/
     1 9.131388579526912D-01,
     1 4.989017152913699D-01,
     1 8.795942675593887D-02,
     1 8.795942675593887D-02,
     1 4.989017152913699D-01,
     1 9.131388579526912D-01/
c
      data extranodes8/
     1 -3.998861349951123D+00,
     1 -2.980147933889640D+00,
     1 -1.945288592909266D+00,
     1 -1.027856640525646D+00,
     1 -3.967966533375878D-01,
     1 -9.086744584657729D-02,
     1 -6.531815708567918D-03,
     1 6.531815708567918D-03,
     1 9.086744584657729D-02,
     1 3.967966533375878D-01,
     1 1.027856640525646D+00,
     1 1.945288592909266D+00,
     1 2.980147933889640D+00,
     1 3.998861349951123D+00/
      data extraweights8/
     1 1.004787656533285D+00,
     1 1.036093649726216D+00,
     1 1.008710414337933D+00,
     1 7.947291148621894D-01,
     1 4.609256358650077D-01,
     1 1.701315866854178D-01,
     1 2.462194198995203D-02,
     1 2.462194198995203D-02,
     1 1.701315866854178D-01,
     1 4.609256358650077D-01,
     1 7.947291148621894D-01,
     1 1.008710414337933D+00,
     1 1.036093649726216D+00,
     1 1.004787656533285D+00/
c
      data extranodes16/
     1 -8.999998754306120d+00,
     1 -7.999888757524622d+00,
     1 -6.997957704791519d+00,
     1 -5.986360113977494d+00,
     1 -4.957203086563112d+00,
     1 -3.928129252248612d+00,
     1 -2.947904939031494d+00,
     1 -2.073471660264395d+00,
     1 -1.348993882467059d+00,
     1 -7.964747731112430d-01,
     1 -4.142832599028031d-01,
     1 -1.805991249601928d-01,
     1 -6.009290785739468d-02,
     1 -1.239382725542637d-02,
     1 -8.371529832014113d-04,
     1 8.371529832014113d-04,
     1 1.239382725542637d-02,
     1 6.009290785739468d-02,
     1 1.805991249601928d-01,
     1 4.142832599028031d-01,
     1 7.964747731112430d-01,
     1 1.348993882467059d+00,
     1 2.073471660264395d+00,
     1 2.947904939031494d+00,
     1 3.928129252248612d+00,
     1 4.957203086563112d+00,
     1 5.986360113977494d+00,
     1 6.997957704791519d+00,
     1 7.999888757524622d+00,
     1 8.999998754306120d+00/
      data extraweights16/
     1 1.000007149422537d+00,
     1 1.000395017352309d+00,
     1 1.004798397441514d+00,
     1 1.020308624984610d+00,
     1 1.035167721053657d+00,
     1 1.014359775369075d+00,
     1 9.362411945698647d-01,
     1 8.051212946181061d-01,
     1 6.401489637096768d-01,
     1 4.652220834914617d-01,
     1 3.029123478511309d-01,
     1 1.704889420286369d-01,
     1 7.740135521653088d-02,
     1 2.423621380426338d-02,
     1 3.190919086626234d-03,
     1 3.190919086626234d-03,
     1 2.423621380426338d-02,
     1 7.740135521653088d-02,
     1 1.704889420286369d-01,
     1 3.029123478511309d-01,
     1 4.652220834914617d-01,
     1 6.401489637096768d-01,
     1 8.051212946181061d-01,
     1 9.362411945698647d-01,
     1 1.014359775369075d+00,
     1 1.035167721053657d+00,
     1 1.020308624984610d+00,
     1 1.004798397441514d+00,
     1 1.000395017352309d+00,
     1 1.000007149422537d+00/
c
       ima=(0,1)
       done=1
       pi=4*atan(done)
c
       if (norder.eq. 0) then
         nskip=1
         nextra=0
         goto 1111
       endif

       if (norder.eq.2) then
         nskip = 1
         nextra = 2
         do i=1,nextra
	      extranodes(i) = extranodes2(i)
	      extraweights(i) = extraweights2(i)
         enddo
         goto 1111
       endif

       if (norder.eq.3) then
         nskip = 2
         nextra = 4
         do i=1,nextra
	       extranodes(i) = extranodes3(i)
	       extraweights(i) = extraweights3(i)
         enddo
        goto 1111
      endif
c
      if (norder.eq.4) then
        nskip = 2
        nextra = 6
	    do i = 1,nextra
	      extranodes(i) = extranodes4(i)
	      extraweights(i) = extraweights4(i)
        enddo
        goto 1111
      endif
c
      if (norder.eq.8) then 
        nskip = 5
        nextra = 14
	    do i = 1,nextra
	      extranodes(i) = extranodes8(i)
	      extraweights(i) = extraweights8(i)
        enddo
        goto 1111
      endif
c
      if (norder.eq.16) then 
        nskip = 10
        nextra = 30
	    do i = 1,nextra
	      extranodes(i) = extranodes16(i)
	      extraweights(i) = extraweights16(i)
        enddo
        goto 1111
      endif      

c
c       quadrature code not available for any other order -> return
c
      ier = 1
      call prinf('wrong order for quadrature, norder=*',norder,1)
      stop
      return
 1111 continue
      ier = 0

c
c       carry out "punctured" trapezoidal rule and fill in matrix
c       entries, skipping entries within nskip of the diagonal
c

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j)
        do j=1,ns
          do i=1,ns
            amat(i,j)=0
          enddo
        enddo
C$OMP END PARALLEL DO       


      imax = 0
      imin = ns

c
      n=ns-2*nskip+1
c
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(i,iii,k,u) 
      do i=1,ns
        iii=i-1+nskip

        do k=0,n-1
          iii=iii+1
          if (iii .gt. ns) iii=iii-ns
          call gfun(srcinfo(1,iii),srcinfo(1,i),dpars,zpars,ipars,u)
          amat(i,iii)=u*srcinfo(5,iii)*h
        enddo
      enddo
C$OMP END PARALLEL DO

        if (norder .eq. 0) return

c
c       now add in corrections and interpolated stuff for alpert
c       first determine all the interpolation coefficients
c
        ninterp=norder+2
c
      do ipt=1,ns
        do i=1,nextra
          txtra(i)=h*(ipt-1)+h*extranodes(i)
        enddo

        do i=1,nextra
c
c       find the closest ninterp points to each of the txtra
c
          n1=txtra(i)/h
          if (txtra(i) .lt. 0) n1=n1-1
          n2=n1+1
          nnn=n1-(ninterp-2)/2

c
          do j=1,ninterp
            its(j)=nnn+j-1
            its2(j)=its(j)+1
            if (its2(j) .le. 0) its2(j)=its2(j)+ns
            if (its2(j) .gt. ns) its2(j)=its2(j)-ns
          enddo

c
c       fill interpolation nodes and function values
c
          do j=1,ninterp
            tpts(j)=its(j)*h
            do l=1,5
              srcpts(l,j)=srcinfo(l,its2(j))
            enddo
          enddo

c
c       now compute the values of xs, ys, dsdt at ttt using barycentric
c       interpolation
c
          ttt=txtra(i)
          call bary1_coefs(ninterp,tpts,ttt,coefs)
          do l=1,5
            srctmp(l)=0
          enddo
c
          do j=1,ninterp
            do l=1,5
              srctmp(l)=srctmp(l)+srcpts(l,j)*coefs(j)
            enddo
          enddo


c
c       evaluate the kernel at the new quadrature point xxx,yyy and
c       add its contribution to the matrix at its interpolation points
c

          call gfun(srctmp,srcinfo(1,ipt),dpars,zpars,ipars,u)
c
          do j=1,ninterp
            jjj=its2(j)
            amat(ipt,jjj)=amat(ipt,jjj)+u*srctmp(5)*h*
     1          extraweights(i)*coefs(j)
          enddo
        enddo
      enddo

c
      return
      end

c
c
c
c
c
c
        subroutine formmatbac_corr(acorr,norder,ns,srcinfo,h,gfun,dpars,
     1      zpars,ipars)
        implicit real *8 (a-h,o-z)
        real *8 srcinfo(6,ns),dpars(*)
        integer ipars(*)
        complex *16 zpars(*)
        real *8 tpts(100),txtra(100),coefs(100),srctmp(5),srcpts(5,100)
        complex *16 acorr(*),ima,u
        integer *4 its(100),its2(100)
        dimension extranodes2(2), extraweights2(2)
        dimension extranodes3(4), extraweights3(4)
        dimension extranodes4(6), extraweights4(6)
        dimension extranodes8(14), extraweights8(14)
        dimension extranodes16(30), extraweights16(30)
        dimension extranodes(30), extraweights(30)
        external gfun
c
c       this routine builds the local correction for  
c       the logarithimacally singular kernels defined via 
c       gfun to a vector using alpert quadrature
c
c     input:
c       norder - the order of alpert quadrature to use, 0, 2, 4, 8, 16
c       srcinfo(6,ns) 
c           srcinfo(1:2,*) = xs,ys, srcinfo(3:4,*) = rnx,rny,
c           srcinfo(5,*) = dsdt
c           srcinfo(6,:) - curvature
c       xs,ys - the x,y coordinates of points on the curve
c       rnx,rny - the unit normal derivative at the points xs,ys
c       dsdt - derivative of parameterization
c       h - step size in parameterization variable
c       ns - length of xs and ys
c       gfun - routine that evaluates the kernel, must have a
c           calling sequence of the form
c           gfun(srcinfo,targinfo,dpars,zpars,ipars,uval)
c
c           where dpars - list of real parameters
c                 zpars - list of complex parameters
c                 ipars - list of integer parameters
c                 srcinfo - source info
c                 targinfo - target info
c                 uval - the complex *16 value of the kernel 
c                     at ((xtarg,ytarg),(xsrc,ysrc))
c
c     output:
c       acorr - correction matrix that will apply acorr 
c           must be of size
c
c       5*ns for order 2,
c       9*ns for order 4,
c       17*ns for order 8,
c       35*ns for order 16
c
c
      data extranodes2/
     1 -1.591549430918953D-01,
     1 1.591549430918953D-01/
      data extraweights2/
     1 0.5d0,
     1 0.5d0/

      data extranodes3/
     1 -1.150395811972836D-01,
     1 -9.365464527949632D-01,
     1 1.150395811972836D-01,
     1 9.365464527949632D-01/
      data extraweights3/
     1 3.913373788753340d-01,
     1 1.108662621124666d0,
     1 3.913373788753340d-01,
     1 1.108662621124666d0/

      data extranodes4/
     1 -1.023715124251890D+00,
     1 -2.935370741501914D-01,
     1 -2.379647284118974D-02, 
     1 2.379647284118974D-02, 
     1 2.935370741501914D-01,
     1 1.023715124251890D+00/ 
      data extraweights4/
     1 9.131388579526912D-01,
     1 4.989017152913699D-01,
     1 8.795942675593887D-02,
     1 8.795942675593887D-02,
     1 4.989017152913699D-01,
     1 9.131388579526912D-01/
c
      data extranodes8/
     1 -3.998861349951123D+00,
     1 -2.980147933889640D+00,
     1 -1.945288592909266D+00,
     1 -1.027856640525646D+00,
     1 -3.967966533375878D-01,
     1 -9.086744584657729D-02,
     1 -6.531815708567918D-03,
     1 6.531815708567918D-03,
     1 9.086744584657729D-02,
     1 3.967966533375878D-01,
     1 1.027856640525646D+00,
     1 1.945288592909266D+00,
     1 2.980147933889640D+00,
     1 3.998861349951123D+00/
      data extraweights8/
     1 1.004787656533285D+00,
     1 1.036093649726216D+00,
     1 1.008710414337933D+00,
     1 7.947291148621894D-01,
     1 4.609256358650077D-01,
     1 1.701315866854178D-01,
     1 2.462194198995203D-02,
     1 2.462194198995203D-02,
     1 1.701315866854178D-01,
     1 4.609256358650077D-01,
     1 7.947291148621894D-01,
     1 1.008710414337933D+00,
     1 1.036093649726216D+00,
     1 1.004787656533285D+00/
c
      data extranodes16/
     1 -8.999998754306120d+00,
     1 -7.999888757524622d+00,
     1 -6.997957704791519d+00,
     1 -5.986360113977494d+00,
     1 -4.957203086563112d+00,
     1 -3.928129252248612d+00,
     1 -2.947904939031494d+00,
     1 -2.073471660264395d+00,
     1 -1.348993882467059d+00,
     1 -7.964747731112430d-01,
     1 -4.142832599028031d-01,
     1 -1.805991249601928d-01,
     1 -6.009290785739468d-02,
     1 -1.239382725542637d-02,
     1 -8.371529832014113d-04,
     1 8.371529832014113d-04,
     1 1.239382725542637d-02,
     1 6.009290785739468d-02,
     1 1.805991249601928d-01,
     1 4.142832599028031d-01,
     1 7.964747731112430d-01,
     1 1.348993882467059d+00,
     1 2.073471660264395d+00,
     1 2.947904939031494d+00,
     1 3.928129252248612d+00,
     1 4.957203086563112d+00,
     1 5.986360113977494d+00,
     1 6.997957704791519d+00,
     1 7.999888757524622d+00,
     1 8.999998754306120d+00/
      data extraweights16/
     1 1.000007149422537d+00,
     1 1.000395017352309d+00,
     1 1.004798397441514d+00,
     1 1.020308624984610d+00,
     1 1.035167721053657d+00,
     1 1.014359775369075d+00,
     1 9.362411945698647d-01,
     1 8.051212946181061d-01,
     1 6.401489637096768d-01,
     1 4.652220834914617d-01,
     1 3.029123478511309d-01,
     1 1.704889420286369d-01,
     1 7.740135521653088d-02,
     1 2.423621380426338d-02,
     1 3.190919086626234d-03,
     1 3.190919086626234d-03,
     1 2.423621380426338d-02,
     1 7.740135521653088d-02,
     1 1.704889420286369d-01,
     1 3.029123478511309d-01,
     1 4.652220834914617d-01,
     1 6.401489637096768d-01,
     1 8.051212946181061d-01,
     1 9.362411945698647d-01,
     1 1.014359775369075d+00,
     1 1.035167721053657d+00,
     1 1.020308624984610d+00,
     1 1.004798397441514d+00,
     1 1.000395017352309d+00,
     1 1.000007149422537d+00/
c
       ima=(0,1)
       done=1
       pi=4*atan(done)

c
       if (norder.eq. 0) then
         nskip=1
         nextra=0
         goto 1111
       endif

       if (norder.eq.2) then
         nskip = 1
         nextra = 2
         ncorr = 2
         do i=1,nextra
	      extranodes(i) = extranodes2(i)
	      extraweights(i) = extraweights2(i)
         enddo
         goto 1111
       endif

c
      if (norder.eq.4) then
        nskip = 2
        nextra = 6
        ncorr = 4
	    do i = 1,nextra
	      extranodes(i) = extranodes4(i)
	      extraweights(i) = extraweights4(i)
        enddo
        goto 1111
      endif
c
      if (norder.eq.8) then 
        nskip = 5
        nextra = 14
        ncorr = 8
	    do i = 1,nextra
	      extranodes(i) = extranodes8(i)
	      extraweights(i) = extraweights8(i)
        enddo
        goto 1111
      endif
c
      if (norder.eq.16) then 
        nskip = 10
        nextra = 30
        ncorr = 17
	    do i = 1,nextra
	      extranodes(i) = extranodes16(i)
	      extraweights(i) = extraweights16(i)
        enddo
        goto 1111
      endif      

c
c       quadrature code not available for any other order -> return
c
      ier = 1
      call prinf('wrong order for quadrature, norder=*',norder,1)
      stop
      return
 1111 continue
      ier = 0

c
c       carry out "punctured" trapezoidal rule and fill in matrix
c       entries, skipping entries within nskip of the diagonal
c

      nn = ns*(2*ncorr+1)
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(j)
        do j=1,nn
          acorr(j) = 0
        enddo
C$OMP END PARALLEL DO       


c
      n=2*ncorr+1


c
      do i=1,ns
        do j=-ncorr,-nskip
          ipt = i+j
          if(ipt.le.0) ipt = ipt + ns
          if(ipt.gt.ns) ipt = ipt - ns

          icorr = (i-1)*n + j+ncorr+1
          call gfun(srcinfo(1,ipt),srcinfo(1,i),dpars,zpars,ipars,u)
          acorr(icorr) = u*srcinfo(5,ipt)*h
        enddo

        do j=nskip,ncorr
          ipt = i+j
          if(ipt.le.0) ipt = ipt + ns
          if(ipt.gt.ns) ipt = ipt - ns

          icorr = (i-1)*n + j+ncorr+1
          call gfun(srcinfo(1,ipt),srcinfo(1,i),dpars,zpars,ipars,u)
          acorr(icorr) = u*srcinfo(5,ipt)*h
        enddo
      enddo


c
c       now add in corrections and interpolated stuff for alpert
c       first determine all the interpolation coefficients
c
      ninterp=norder+2
c
      do ipt=1,ns
        do i=1,nextra
          txtra(i)=h*(ipt-1)+h*extranodes(i)
        enddo

        do i=1,nextra
c
c       find the closest ninterp points to each of the txtra
c
          n1=txtra(i)/h
          if (txtra(i) .lt. 0) n1=n1-1
          n2=n1+1
          nnn=n1-(ninterp-2)/2

c
          do j=1,ninterp
            its(j)=nnn+j-1
            its2(j)=its(j)+1
            if (its2(j) .le. 0) its2(j)=its2(j)+ns
            if (its2(j) .gt. ns) its2(j)=its2(j)-ns
          enddo

c
c       fill interpolation nodes and function values
c
          do j=1,ninterp
            tpts(j)=its(j)*h
            do l=1,5
              srcpts(l,j)=srcinfo(l,its2(j))
            enddo
          enddo

c
c       now compute the values of xs, ys, dsdt at ttt using barycentric
c       interpolation
c
          ttt=txtra(i)
          call bary1_coefs(ninterp,tpts,ttt,coefs)
          do l=1,5
            srctmp(l)=0
          enddo
c
          do j=1,ninterp
            do l=1,5
              srctmp(l)=srctmp(l)+srcpts(l,j)*coefs(j)
            enddo
          enddo


c
c       evaluate the kernel at the new quadrature point xxx,yyy and
c       add its contribution to the matrix at its interpolation points
c
c
          call gfun(srctmp,srcinfo(1,ipt),dpars,zpars,ipars,u)
c
          do j=1,ninterp
            jjj=its2(j)

            icorr = jjj-ipt
            if(icorr.lt.-ncorr) icorr = (icorr + ns)
            if(icorr.gt.ncorr) icorr = icorr-ns

            jcorr = (ipt-1)*n + icorr+ncorr+1
            
            acorr(jcorr) = acorr(jcorr) + u*srctmp(5)*h*
     1          extraweights(i)*coefs(j)
          enddo
        enddo
      enddo

c
      return
      end

c
c
c
c
c

        subroutine bary1_coefs(n,ts,ttt,coefs)
        implicit real *8 (a-h,o-z)
        real *8 ts(1),ttt,whts(1000),coefs(1)
c
c       use barycentric interpolation on ts,xs to evaluate at ttt
c       first evaluate the barycentric weights (real routine)
c
c       input:
c         n - the length of ts,xs
c         ts - nodes with which to interpolate
c         ttt - point at which to interpolate
c
c       output:
c         coefs - coefficients for the interpolation
c
c
        do 1200 i=1,n
        whts(i)=1
 1200   continue
c
        do 1600 i=1,n
        do 1400 j=1,n
        if (i .ne. j) whts(i)=whts(i)/(ts(i)-ts(j))
 1400   continue
 1600   continue
c
c       this uses the '2nd form' of barycentric interpolation, first
c       form the denominator
c
        dd=0
        do 2000 i=1,n
        dd=dd+whts(i)/(ttt-ts(i))
 2000   continue
c
c       and next the interpolation coefficients
c
        do 2400 i=1,n
        coefs(i)=whts(i)/(ttt-ts(i))/dd
 2400   continue
c
        return
        end

