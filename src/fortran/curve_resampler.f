
       subroutine curve_resampler_guru(ier,nx,nh,par1,rl,n,eps,
     1     t,xinfo,h,rltot,w,lw,lsave)
c 
c       This subroutine produces an equispaced (in terms
c       of the arc length) resampling of a user-specified
c       curve. The curve is specified using the routine
c       funcurv_fft whose parameters are supplied by the user above.
c       
c       Explanation: this subroutine
c       assumes that the curve to be resampled is defined by the
c       mapping
c 
c            [0,rl] \to R^2,
c 
c       with the mapping and its derivatives provided via
c       the user-supplied subroutine funcurve (see specifications
c       below). The curve does not need to be closed.
c 
c 
c                    input parameters:
c  nx = number of terms in fourier expansion of x,y,dxdt,dydt
c  nh = number of terms in the fourier expansion of h
c  par1 = fourier coeffs of x+ima*y, dxdt+ima*dydt
c  rl - the length of the t-interval.
c  n - the number of subintervals into which the curve is to be
c       subdivided.
c  eps - the precision with which all operations are to be performed.
c       recommended value: 1.0d-15
c 
c                    output parameters:
c 
c  ier - error return code.
c       ier=0   means normal completion
c       ier=1   means that the discretization might not be
c               sufficiently fine to describe the curve.
c               remedy: increase n. also note that this is
c               an extremely minor problem. what it really
c               says is that the newton process used by
c               the subroutine to find equispaced nodes on
c               the curve had to activate the step-length
c               control procedure at least once.
c       ier=10  means that the adaptive gaussian quadrature
c               used by the subroutine failed to obtain
c               the accuracy requested by the user at least
c               once. this is fairly serious.
c               remedy: increase eps, and/or check your subroutine
c               funcurve
c       ier=100 means that the newton process failed
c               to converge at least once. this is very serious.
c               remedy: increase eps, and/or check your subroutine
c               funcurve
c       ier=16000 means that the length of the user-provided work
c               array, w, is insufficient. this is a fatal error
c  xinfo(4,n) - boundary data at 
c  t - the equispaced discretization of the curve. note that the
c       subroutine will return n+1 values in array t. that is, if the
c       curve is a closed one, then t(n+1) will be equal to t(1).
c  h - the sampling interval along the arc of the resampling returned
c       in array t
c  w - work array; suggested length of 1000 real *8, but see
c       ier=16000 if more space is needed
c  lsave - the number of elements of w that should not be changed
c         between the call to this subroutine and the subsequent calls
c         to anafpt
c
c       . . . determine the total length of the curve
c
      implicit real *8 (a-h,o-z)
      integer nx,nh
      complex *16 par1(*)
      real *8 rl,eps
      integer n
      real *8 h,w(lw),par2(4)
      real *8 xinfo(4,n),t(n+1)
      external funcurv_fft

      par2(1) = nx
      par2(2) = nh
      par2(3) = max(nx,nh)
      par2(4) = rl

      do i=1,n+1
        t(i) = 0
      enddo

      rltot = 0
      h = 0
      lsave = 0

      call anafast(ier,funcurv_fft,par1,par2,rl,n,eps,t,h,rltot,
     1       w,lw,lsave)
      do i=1,n
        call funcurv_fft(t(i),par1,par2,xinfo(1,i),xinfo(2,i),
     1    xinfo(3,i),xinfo(4,i))
        dst = sqrt(xinfo(3,i)**2 + xinfo(4,i)**2)
        xinfo(3,i) = xinfo(3,i)/dst
        xinfo(4,i) = xinfo(4,i)/dst
      enddo

      return
      end
c 
c 
c 
c 
c 
        subroutine funcurv_fft(t,par1,par2,x,y,dxdt,dydt)
c
c  this subroutine generates the curve info for a curve
c  described using it's fourier coeffs
c
c  the curve is assumed to be of the form
c  
c  x(t) + iy(t) = x0(t) + i y0(t)  + h(t)*(dy0dt(t) - i dx0dt(t))
c
c  where x0+iy0 = \sum_{j=nxmin}^{nxmax} c_{j} e^{2\pi j t/L}
c        dx0 + idy0 = \sum_{j=nxmin}^{nxmax} d_{j} e^{2\pi j t/L}
c
c  and h(t) = \sum_{j=nhmin}^{nhmax} h_{j} e^{2\pi j t/L}
c     and dh/dt = \sum_{j=nhmin}^{nhmax} dh_{j} e^{2\pi j t/L}
c
c  nxmin = -(nx-1)/2 if nx is odd and -nx/2 +1 if nx is even
c  nxmax = (nx-1)/2 if nx is odd and nx/2 if nx is even
c
c  similarly for nhmin and nhmax
c
c  par1(1:nx) = coeffs cj ordering is c_{0},c_{1},c_{2}... c_{npos},
c    c_{-nneg},c_{-nneg+1},...c_{-1} 
c  par1(nx+1:2*nx) = coefs d_{j}
c  par1(2*nx + 1:2*nx +nh) = coefs h_{j}
c  par1(2*nx+ nh + 1:2*nx+2*nh) = coefs dh_{j}
c
c  par2(1) = nx
c  par2(2) = nh
c  par2(3) = nmax = max(nx,nh)
c  par2(4) = L
c


c
        implicit real *8 (a-h,o-z)
        complex *16 par1(*)
        complex *16, allocatable :: zexp(:)
        complex *16 ztmp,ima,z,dzdt,d2zdt2,zh,dzhdt
        real *8 par2(*),pi,rl
        data pi/3.141592653589793d0/
        data ima/(0.0d0,1.0d0)/
        
        
        nx = par2(1)
        nh = par2(2)
        nmax = par2(3)
        rl= par2(4)

        nn = nmax/2+1
        
        allocate(zexp(0:nn))

        zexp(0) = 1
        ztmp = exp(pi*ima*2*t/rl)
        do i=1,nn
          zexp(i) = zexp(i-1)*ztmp
        enddo


        x = 0
        y = 0
        dxdt = 0
        dydt = 0
        z = 0
        dzdt = 0 
        d2zdt2 = 0

        zh = 0
        dzhdt = 0

        nxmax = nx/2
        nxmin = (nx-1)/2

        nhmax = nh/2
        nhmin = (nh-1)/2

        do i=0,nxmax
          z = z+par1(i+1)*zexp(i)
        enddo

        do i=1,nxmin
          z = z+conjg(zexp(i))*par1(nx-i+1)
        enddo


        do i=0,nxmax
          dzdt = dzdt+par1(nx+i+1)*zexp(i)
          d2zdt2 = d2zdt2 + par1(nx+i+1)*zexp(i)*2*pi*i/rl*ima
        enddo

        do i=1,nxmin
          dzdt = dzdt+conjg(zexp(i))*par1(2*nx-i+1)
          d2zdt2 = d2zdt2-conjg(zexp(i))*par1(2*nx-i+1)*2*pi*i/rl*ima
        enddo

        

        do i=0,nhmax
          zh = zh+par1(i+1+2*nx)*zexp(i)
        enddo

        do i=1,nhmin
          zh = zh+conjg(zexp(i))*par1(nh-i+1+2*nx)
        enddo

        do i=0,nhmax
          dzhdt = dzhdt+par1(2*nx+nh+i+1)*zexp(i)
        enddo

        do i=1,nhmin
          dzhdt = dzhdt+conjg(zexp(i))*par1(2*nh-i+1+2*nx)
        enddo

        rh = real(zh)
        rdh = real(dzhdt)
        ztmp = z - rh*ima*dzdt/abs(dzdt)
cc        ztmp = z + rh*ima*dzdt

        x = real(ztmp)
        y = imag(ztmp)
        rdot = real(dzdt)*real(d2zdt2) + imag(dzdt)*imag(d2zdt2)

cc        ztmp = dzdt+rh*ima*d2zdt2+rdh*ima*dzdt 
        ztmp = dzdt-rh*ima*d2zdt2/abs(dzdt)-rdh*ima*dzdt/abs(dzdt) 
     1    + rh*ima*dzdt*rdot/abs(dzdt)**3
        dxdt = real(ztmp)
        dydt = imag(ztmp)
        
        return
        end
c 
c 
c 
c 
c 
c
        subroutine anafast(ier,funcurve,par1,par2,rl,n,eps,t,h,rltot,
     1       w,lw,lsave)
c
        implicit real *8 (a-h,o-z)
        real *8 t(1),par1(1),par2(1),w(lw)
        external funcurve
c 
c       This subroutine produces an equispaced (in terms
c       of the arc length) resampling of a user-specified
c       curve. It is an updated version of anaresa that should run
c       in roughly 1/3 the time. The curve is specified in the form
c       of the subroutine funcurve providing some (generally, not
c       equispaced) parametrization of the curve and the length rl of
c       the interval of parametrization. Explanation: this subroutine
c       assumes that the curve to be resampled is defined by the
c       mapping
c 
c            [0,rl] \to R^2,
c 
c       with the mapping and its derivatives provided via
c       the user-supplied subroutine funcurve (see specifications
c       below). The curve does not need to be closed.
c 
c 
c                    input parameters:
c 
c  funcurve - the subroutine providing a parametrization
c        of the curve. the calling sequence of funcurve
c        must be
c 
c                  funcurve(t,par1,par2,x,y,dxdt,dydt),
c 
c        with:
c 
c       t - the parameter by which the user's curve is parametrized
c       x - the x-coordinate of the point on the curve corresponding
c           to the values t of the parameter
c       y - the y-coordinate of the point on the curve corresponding
c           to the values t of the parameter
c       dxdt - the derivative with respect to t of the  x-coordinate
c           of the point on the curve corresponding to the values t
c           of the parameter
c       dydt - the derivative with respect to t of the  y-coordinate
c           of the point on the curve corresponding to the values t
c           of the parameter
c
c  rl - the length of the t-interval.
c  n - the number of subintervals into which the curve is to be
c       subdivided.
c  eps - the precision with which all operations are to be performed.
c       recommended value: 1.0d-15
c 
c                    output parameters:
c 
c  ier - error return code.
c       ier=0   means normal completion
c       ier=1   means that the discretization might not be
c               sufficiently fine to describe the curve.
c               remedy: increase n. also note that this is
c               an extremely minor problem. what it really
c               says is that the newton process used by
c               the subroutine to find equispaced nodes on
c               the curve had to activate the step-length
c               control procedure at least once.
c       ier=10  means that the adaptive gaussian quadrature
c               used by the subroutine failed to obtain
c               the accuracy requested by the user at least
c               once. this is fairly serious.
c               remedy: increase eps, and/or check your subroutine
c               funcurve
c       ier=100 means that the newton process failed
c               to converge at least once. this is very serious.
c               remedy: increase eps, and/or check your subroutine
c               funcurve
c       ier=16000 means that the length of the user-provided work
c               array, w, is insufficient. this is a fatal error
c  t - the equispaced discretization of the curve. note that the
c       subroutine will return n+1 values in array t. that is, if the
c       curve is a closed one, then t(n+1) will be equal to t(1).
c  h - the sampling interval along the arc of the resampling returned
c       in array t
c  w - work array; suggested length of 1000 real *8, but see
c       ier=16000 if more space is needed
c  lsave - the number of elements of w that should not be changed
c         between the call to this subroutine and the subsequent calls
c         to anafpt
c
c       . . . determine the total length of the curve
c 
        ier=0
        jer=0
        ker=0
c
        a=0
        b=rl
        m=6
c
        epsgauss=eps
        epsnewt=sqrt(eps)/1000
c
c       setup work array
c       reserve first couple elements for later
c
        ixg=10
        lxg=m+3
c
        iwg=ixg+lxg
        lwg=m+3
c
        iwr=iwg+lwg
        lwr=(lw-iwr)/2-5
c
        iwv=iwr+lwr
        lwv=lwr
c
c       compute gaussian nodes and weights
c
        itype=1
        call legeexps(itype,m,w(ixg),unil,vnil,w(iwg))
c
        call anafastga(ier,a,b,funcurve,par1,par2,m,w(ixg),w(iwg),
     1       epsgauss,rltot,w(iwr),w(iwv),lwr,lwv,maxrec,numint,nints)
        if (ier.ne.0) return
cccc        call prin2('after first anafastga, rltot=*',rltot,1)
cccc        call prinf('after first anafastga, maxrec=*',maxrec,1)
c
c       compress and finish work array
c
        lwr=nints+3
        iwv2=iwr+lwr
        lwv2=lwr
c
        lsave=iwv2+lwv2
c
        do 100 i=0,lwv2-1
        w(iwv2+i)=w(iwv+i)
 100    continue
c
        w(1)=rl
        w(2)=m
        w(3)=nints
        w(4)=ixg
        w(5)=iwg
        w(6)=iwr
        w(7)=iwv2
c 
c       one subinterval after another, construct the equispaced
c       subdivison of the user-specified curve
c 
        h=rltot/n
        t(1)=0
        ij0=1
c
        dtp0=0
        dtp1=0
        dtp2=0
        dtp3=0
        dtp4=0
        dtpp0 = 0

c
        do 3000 k=1,n
cccc         call prinf('in anafast, k=*',k,1)
c 
c       determine the dl/dt at the point t(k)
c 
        call funcurve(t(k),par1,par2,x,y,dxdt,dydt)
c 
c       get the initial approximation for the next point
c 
        dldt=sqrt(dxdt**2+dydt**2)
        dt0=h/dldt
c
        if (k.le.5) goto 1100
c
c       after the first couple iterations
c       initialize newton by constructing
c       cubic extrapolation of previous dt's
c

   
        ddt0=dtp0
        ddt4=5*dtp0-10*dtp1+10*dtp2-5*dtp3+dtp4
        if(ddt4.le.0) goto 1100

c
        ddd1=abs(dtp0-dtpp0)
        ddd2=abs(ddt4-ddt0)
c
        rdd=ddd2/ddd1
cccc        call prin2_long('rdd = *',rdd,1)
c
        dtpp0=dt0
c
c       protect against cases when cubic extrapolation
c       is not a good choice
c
        if (rdd.gt.4d0) goto 1100
c
        dt0=ddt4
c
 1100   continue
c 
c       use newton to obtain the next point on the curve
c       accurately
c
        dt=dt0
        maxit=20
        do 1400 i=1,maxit
c 
c       determine the actual length of the step
c       - with step-length control
c 
        dtold=0
        derrold=h
        do 1200 j=1,20
c 
c       if the point has jumped outside the interval
c       of definition of the curve, or jumped backwards
c       - bring it back within  reasonable range
c 
        if(t(k)+dt.gt.rl) dt=rl-t(k)
        if(dt.lt.0) dt=0
c 
c       perform step-length control
c
        call anafindint(w(iwr),nints,t(k)+dt,ij0,ijk)
c
        if (ij0.eq.ijk) then
        call anafastg2(t(k),t(k)+dt,funcurve,
     1       par1,par2,w(ixg),w(iwg),m,htrue2)
        goto 2300
        end if
c
        call anafastg2(t(k),w(iwr+ij0-1),funcurve,
     1       par1,par2,w(ixg),w(iwg),m,htrue2)
c
        do 2200 kk=ij0+1,ijk-1
        htrue2=htrue2+w(iwv+kk-1)
 2200   continue
c
        call anafastg2(w(iwr+ijk-2),t(k)+dt,funcurve,
     1       par1,par2,w(ixg),w(iwg),m,htrue3)
        htrue2=htrue2+htrue3
c
 2300   continue
        htrue=htrue2
c 
c       determine the dl/dt at that point and make a
c       newton step
c 
        call funcurve(t(k)+dt,par1,par2,x,y,dxdt,dydt)
        dldt=sqrt(dxdt**2+dydt**2)
c 
        derr=htrue-h
cccc        call prin2('htrue-h=*',derr,1)
c 
        if (abs(derr).lt.abs(derrold)) goto 1300
        dt=(dt+dtold)/2
cccc        call prin2('step-length control activated, dt=*',dt,1)
        ker=1
 1200   continue
 1300   continue
c 
        dt=dt-derr/dldt
c 
c        if the required precision has been achieved
c        -terminate newton
c 
        if (abs(derr).lt.epsnewt) goto 2000
c 
 1400   continue
cccc        call prin2('newton failed, derr=*',derr,1)
        ier=1
 2000   continue
c
c       store previous dt's
c
        dtp4=dtp3
        dtp3=dtp2
        dtp2=dtp1
        dtp1=dtp0
        dtp0=dt
c
        ij0=ijk
        t(k+1)=t(k)+dt
c
 3000   continue

c 
c       if(jer.ne.0) jer=1
        ier=ier*100+jer*10+ker
cccc        call prinf('exiting anafast, ier=*',ier,1)
c
        return
c
        end
c
c
c
c
c
        subroutine anafindint(w,n,t,ii,ij)
c
        implicit real *8 (a-h,o-z)
        real *8 w(1)
c
c       given sorted array w, find interval in which
c       t lives; ii is initialization of this search
c       ij is output, i.e. w(ij-1)< t <=w(ij)
c
        do 100 i=ii,n
        ij=i
        if (t.le.w(i)) return
 100    continue
c
        return
        end
c 
c 
c 
c 
c 
        subroutine anafpt(ier,x,n,h,ts,funcurve,par1,par2,eps,
     1       tout,xout,yout,dxdtout,dydtout,w)
c
        implicit real *8 (a-h,o-z)
        real *8 ts(1),par1(1),par2(1),w(1),tz(10),xz(10),cz(10)
        external funcurve
c 
c       This subroutine finds the location (both on the
c       parametrizing interval [0,rl] and in R^2) of the
c       point specified by the user on the curve. The
c       curve has to have been preprocessed by a preceding
c       call to the subroutine anafast (see), and the point
c       is specified by its distance (in terms of arc-length)
c       from the beginning of the curve. Please note that this
c       subroutine depends on the subroutine anafast for the
c       parameteres h,ts and the work array. This subroutine
c       has no function as a stand-alone device.
c 
c                        input parameters:
c
c  x - the location (in terms of arc-length) on the curve of the
c       point to be found
c  n - the number of nodes in the equispaced discretization of the
c       curve produced by a preceding call to anafast
c  h - the sampling interval along the arc of the resampling returned
c       in array ts by the subroutine anafast
c  ts - the equispaced discretization of the curve produced by the
c       subroutine anafast
c  funcurve - the subroutine providing a parametrization
c        of the curve. the calling sequence of funcurve
c        must be
c  w - work array returned by subroutine anafast
c 
c                  funcurve(t,par1,par2,x,y,dxdt,dydt),
c 
c        with:
c 
c       t - the parameter by which the user's curve is parametrizaed
c       x - the x-coordinate of the point on the curve corresponding
c           to the values t of the parameter
c       y - the y-coordinate of the point on the curve corresponding
c           to the values t of the parameter
c       dxdt - the derivative with respect to t of the  x-coordinate
c           of the point on the curve corresponding to the values t
c           of the parameter
c       dydt - the derivative with respect to t of the  y-coordinate
c           of the point on the curve corresponding to the values t
c           of the parameter
c  eps - the precision to which the point is to be found. Generally, it
c       is a good idea to make this parameter the same as in the
c       preceding call to anafast
c 
c                         output parameters:
c 
c  tout - the value of the curve parameter corresponding to the location
c       x along the arc-length
c  xout, yout - the coordinates in R^2 of the point x on the curve
c  dxdt, dydt - the coordinates of the tangent (of length 1) to the
c       curve at the point (xout,yout)
c

c
c       unwrap work array
c
        rl=w(1)
        mm=w(2)
        nints=w(3)
        ixg=w(4)
        iwg=w(5)
        iwr=w(6)
        iwv=w(7)
c 
c       . . . find the two nodes (in  arc-length) between which the
c       user-supplied point is located
c 
        ier=0
        d=x/h
        m=d
c 
c       if the point is outside the user-supplied curve - bomb
c
        rn=n
        if (d.gt.rn) ier=512
        if (d.gt.rn) return
c
c       if point is end-point of the user-supplied curve - use rl
c
        if (m.eq.n) tk=rl
        if (m.eq.n) goto 3000
c
        tinit=ts(m+1)
        xinit=m*h
c
c       use cubic interpolation to find the approximate location
c       of the curve parameter (as opposed to the arc-length)
c       corresponding to the user-supplied arc-length
c
        if (m.ge.1) tz(1)=ts(m)
        if (m.eq.0) tz(1)=ts(n)-ts(n+1)
        tz(2)=ts(m+1)
        tz(3)=ts(m+2)
        if (m.le.n-2) tz(4)=ts(m+3)
        if (m.eq.n-1) tz(4)=ts(n+1)+ts(2)
c
        xz(1)=(m-1)*h
        xz(2)=m*h
        xz(3)=(m+1)*h
        xz(4)=(m+2)*h
c
        cz(1)=1
        cz(2)=1
        cz(3)=1
        cz(4)=1
c
        do 100 i=1,4
        do 200 j=1,4
        if (i.eq.j) goto 200
        cz(i)=cz(i)*(x-xz(j))/(xz(i)-xz(j))
 200    continue
 100    continue
c
        t0=0
        do 300 i=1,4
        t0=t0+cz(i)*tz(i)
 300    continue
c 
c       use Newton to find the location of the curve parameter
c       (as opposed to the arc-length) corresponding to the
c       user-supplied arc-length
c 
        tk=t0
        tkold=t0
c
        epsgauss=eps
        epsnewt=sqrt(eps)/1000
c
        maxit=20
        fk=h
c
        ij0=1
        call anafindint(w(iwr),nints,tinit,ij0,ija)
c 
        do 1400 i=1,maxit
c 
c       evaluate the function whose root we are seeking
c 
        step=1
        do 1200 j=1,20
c
        call anafindint(w(iwr),nints,tk,ija,ijb)
c
        if (ija.eq.ijb) then
        call anafastg2(tinit,tk,funcurve,par1,par2,
     1       w(ixg),w(iwg),mm,rltot)
        goto 2300
        end if
c
        call anafastg2(tinit,w(iwr+ija-1),funcurve,par1,par2,
     1       w(ixg),w(iwg),mm,rltot)
c
        do 2200 kk=ija+1,ijb-1
        rltot=rltot+w(iwv+kk-1)
 2200   continue
c
        call anafastg2(w(iwr+ijb-2),tk,funcurve,par1,par2,
     1       w(ixg),w(iwg),mm,rltot2)
        rltot=rltot+rltot2
c
 2300   continue
c 
        fkold=fk
        fk=rltot+xinit-x
c 
c       if the target function failed to decrease
c       - activate step-length control
c
        if (abs(fk).le.abs(fkold)) goto 1100
c 
        ier=1
        step=step/2
        tk=tkold-fkold/dldt*step
        goto 1200
c 
 1100   continue
c 
c       determine the dl/dt at that point and make a
c       newton step
c 
        call funcurve(tk,par1,par2,xk,yk,dxdt,dydt)
        dldt=sqrt(dxdt**2+dydt**2)
c 
        tkold=tk
        tk=tk-fk/dldt*step
c 
c        if the required precision has been achieved
c        - terminate newton
c 
        if (abs(fk).lt.epsnewt) goto 2000
c 
        goto 1400
 1200   continue
c 
        ier=64
        return
c 
 1400   continue
 2000   continue
c
 3000   continue
c
        tout=tk
        call funcurve(tout,par1,par2,xout,yout,dxdtout,dydtout)
c 
        d=sqrt(dxdtout**2+dydtout**2)
        dxdtout=dxdtout/d
        dydtout=dydtout/d
c 
        return
c
        end
c 
c 
c 
c 
c 
        subroutine anafastga(ier,a,b,fun,par1,par2,m,t,w,eps,
     1      rint,wright,wvals,lenwr,lenwv,maxrec,numint,nints)
c
        implicit real *8 (a-h,o-z)
        real *8 t(100),w(100),stack(400),vals(200),idep(200),
     1      par1(1),par2(1),wright(1),wvals(1)
c 
c       this subroutine uses the adaptive gauss quadrature
c       to evaluate the length of a user-supplied interval
c       on a curve.
c 
c                       input parameters:
c 
c  a,b - the ends of the interval on which the integral is to
c       be evaluated
c  fun - the user-supplied function describing the curve.to be integrated.
c        the calling sequence of fun must be
c 
c        fun(t,x,y,dxdt,dydt)
c 
c  par1, par2 - dummies
c  m - the order of the quadrature to me used on each subinterval
c  eps - the accuracy (absolute) to which the integral will be
c       evaluated
c 
c                       output parameters:
c 
c  ier - error return code.
c       ier=0  means normal conclusion
c       ier=8  means that at some point, the stack was exceeded. this
c              is a fatal error.
c       ier=16 means that the total number of subintervals in the
c              adaptive subdivision of [a,b] turned out to be greater
c              than 100000. this is a fatal error.
c 
c  rint - the integral as evaluated
c  maxrec - the maximum depth to which the recursion went at its
c       deepest point. can not be greater than 200, since at that
c       point ier is set to 8 and the execution of the subroutine
c       terminated.
c  numint - the total number of intervals in the subdivision. can not
c       be greater than 100000, since at that
c       point ier is set to 16 and the execution of the subroutine
c       terminated.
c 
c       integrate the user-supplied function using the
c       adaptive gaussian quadratures
c 
        nnmax=100000
        maxdepth=200
c
        call anafastg1(ier,stack,a,b,fun,par1,par2,t,w,m,
     1       vals,nnmax,eps,rint,wright,wvals,lenwr,lenwv,nints,idep,
     2       maxdepth,maxrec,numint)
c
        return
c
        end
c 
c 
c 
c 
c 
        subroutine anafastg1(ier,stack,a,b,fun,par1,par2,t,w,m,
     1      vals,nnmax,eps,rint,wright,wvals,lenwr,lenwv,nints,idepth,
     2      maxdepth,maxrec,numint)
c
        implicit real *8 (a-h,o-z)
        real *8 stack(2,1),t(1),w(1),vals(1),par1(1),par2(1),
     1       wright(1),wvals(1),idepth(1)
c 
c       start the recursion
c 
        stack(1,1)=a
        stack(2,1)=b
        idepth(1)=0
        call anafastg2(a,b,fun,par1,par2,t,w,m,vals(1))
c 
c       recursively integrate the thing
c 
        j=1
        rint=0
        ier=0
        maxrec=0
        nints=0
c
        do 3000 i=1,nnmax
cccc        call prinf('i=*',i,1)
c
        numint=i
        if (idepth(j).gt.maxrec) maxrec=idepth(j)
cccc        call prinf('j=*',j,1)
c 
c       subdivide the current subinterval
c
        c=(stack(1,j)+stack(2,j))/2
        call anafastg2(stack(1,j),c,fun,
     1       par1,par2,t,w,m,value2)
c 
        call anafastg2(c,stack(2,j),fun,
     1       par1,par2,t,w,m,value3)
c 
        dd=dabs(value2+value3-vals(j))
cccc        call prin2('in anafastg1, dd=*',dd,1)
c
        ifdone=0
        if (dd.le.eps) ifdone=1
c 
c       if the function on this subinterval has been
c       integrated with sufficient accuracy - add the
c       value to that of the global integral and move up
c       in the stack
c 
        if (ifdone.eq.0) goto 2000
c
        nints=nints+1
c
        if (nints.le.lenwr.and.nints.le.lenwv) goto 1200
        ier=16000
        return
 1200   continue
c
        wright(nints)=stack(2,j)
        wvals(nints)=value2+value3
c 
        rint=rint+value2+value3
        j=j-1
c 
c       if the whole thing has been integrated - return
c 
        if (j.eq.0) return
        goto 3000
c
 2000   continue
c 
c       if the function on this subinterval has not been
c       integrated with sufficient accuracy - move
c       down the stack
c 
        stack(1,j+1)=stack(1,j)
        stack(2,j+1)=(stack(1,j)+stack(2,j))/2
        vals(j+1)=value2
        idepth(j+1)=idepth(j)+1
c 
        stack(1,j)=(stack(1,j)+stack(2,j))/2
        vals(j)=value3
        idepth(j)=idepth(j)+1
c 
        j=j+1
c 
c       if the depth of the recursion has become excessive - bomb
c 
        if (j.le.maxdepth) goto 3000
        ier=8
        return
c
 3000   continue
c
        ier=16
c
        return
c
        end
c 
c 
c 
c 
c 
        subroutine anafastg2(a,b,fun,par1,par2,t,w,m,rint)
c
        implicit real *8 (a-h,o-z)
        real *8 t(1),w(1),par1(1),par2(1)
c 
c       integrate the function fun on the interval [a,b]
c 
        rint=0
        u=(b-a)/2
        v=(b+a)/2
c
        do 1200 i=1,m
c
        tt=u*t(i)+v
c
        call fun(tt,par1,par2,x,y,dxdt,dydt)
c
        rint=rint+w(i)*sqrt(dxdt**2+dydt**2)
c
 1200   continue
c
        rint=rint*u
c
        return
c
        end
