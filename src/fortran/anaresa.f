        subroutine anaresa(ier,funcurve,par1,par2,rl,n,eps,t,h,rltot)
        implicit real *8 (a-h,o-z)
        save
        dimension t(1),par1(1),par2(1)
        external funcurve
c 
c        This subroutine produces an equispaced (in terms
c        of the arc length) resampling of a user-specified
c        curve. The curve is specified in the form of the
c        subroutine funcurve providing some (generally,
c        not equispaced) parametrization of the curve, and
c        the length rl of the interval of parametrization.
c        Explanation: this subroutine assumes that the curve
c        to be resampled is defined by the mapping
c 
c            [0,rl] \to R^2,
c 
c        with the mapping and its derivatives provided via
c        the user-supplied subroutine funcurve (see specifications
c        below). The curve does not need to be closed.
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
c          t - the parameter by which the user's curve is parametrizaed
c          x - the x-coordinate of the point on the curve corresponding
c                to the values t of the parameter
c          y - the y-coordinate of the point on the curve corresponding
c                to the values t of the parameter
c          dxdt - the derivative with respect to t of the  x-coordinate
c                of the point on the curve corresponding to the values t
c                of the parameter
c          dydt - the derivative with respect to t of the  y-coordinate
c                of the point on the curve corresponding to the values t
c                of the parameter
c  rl - the length of the t-interval.
c  n - the number of subintervals into which the curve is to be
c          subdivided.
c  eps - the precision with which all operations are to be performed.
c          recommended value: 1.0d-15
c 
c                    output parameters:
c 
c  ier - error return code.
c                  ier=0 means normal completion
c                  ier=1 means that the discretization might not be
c                        sufficiently fine to describe the curve.
c                    remedy: increase n. also note that this is
c                        an extremely minor problem. what it really
c                        says is that the newton process used by
c                        the subroutine to find equispaced nodes on
c                        the curve had to activate the step-length
c                        control procedure at least once.
c                  ier=10  means that the adaptive gaussian quadrature
c                        used by the subroutine failed to obtain
c                       the accuracy requested by the user at least
c                        once. this is fairly serious.
c                    remedy: increase eps, and/or check your subroutine
c                        funcurve
c                  ier=100 means that the newton process failed
c                        to converge at least once. this is very serious.
c                    remedy: increase eps, and/or check your subroutine
c                        funcurve
c  t - the equispaced discretization of the curve. note that the
c          subroutine will return n+1 values in array t. that is, if the
c          curve is a closed one, then t(n+1) will be equal to t(1).
c  h - the sampling interval along the arc of the resampling returned
c          in array t.
c 
c        . . . determine the total length of the curve
c 
        ier=0
        jer=0
        ker=0
        a=0
        b=rl
        m=12
        epsgauss=eps*1000
        epsnewt=dsqrt(eps)/1000
c 
cccc        call prinf('before anaresga, m=*',m,1)
        call anaresga(jer,a,b,funcurve,par1,par2,m,epsgauss,
     1      rltot,maxrec,numint)
cccc        call prin2('after first anaresga, rltot=*',rltot,1)
c 
c       one subinterval after another, construct the equispaced
c       subdivison of the user-specified curve
c 
        h=rltot/n
        t(1)=0
c 
        do 3000 k=1,n
cccc         call prinf('in anaresa, k=*',k,1)
c 
c       determine the dl/dt at the point t(k)
c 
        call funcurve(t(k),par1,par2,x,y,dxdt,dydt)
c 
c       get the initial approximation for the next point
c 
        dldt=dsqrt(dxdt**2+dydt**2)
        dt0=h/dldt
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
        if(t(k)+dt .gt. rl) dt=rl-t(k)
        if(dt .lt. 0) dt=0
c 
c       perform step-length control
c 
        call anaresga(jer,t(k),t(k)+dt,funcurve,
     1     par1,par2,m,epsgauss,htrue,maxrec,numint)
ccc         call prinf('after anaresga, maxrec=*',maxrec,1)
ccc         call prinf('after anaresga, numint=*',numint,1)
c 
c       determine the dl/dt at that point and make a
c       newton step
c 
        call funcurve(t(k)+dt,par1,par2,x,y,dxdt,dydt)
        dldt=dsqrt(dxdt**2+dydt**2)
c 
        derr=htrue-h
cccc        call prin2('htrue-h=*',derr,1)
c 
        if(dabs(derr) .lt. dabs(derrold) ) goto 1300
        dt=(dt+dtold)/2
cccc          call prin2('step-length control activated, dt=*',dt,1)
        ker=1
 1200 continue
 1300 continue
c 
        dt=dt-derr/dldt
c 
c        if the required precision has been achieved
c        -terminate newton
c 
        if(dabs(derr) .lt. epsnewt) goto 2000
c 
 1400 continue
cccc        call prin2('newton failed, derr=*',derr,1)
        ier=1
 2000 continue
        t(k+1)=t(k)+dt
 3000 continue
c 
c       if(jer .ne. 0) jer=1
        ier=ier*100+jer*10+ker
cccc        call prinf('exiting anaresa, ier=*',ier,1)
        return
        end
c 
c 
c 
c 
c 
         subroutine anapoint(ier,x,n,h,ts,funcurve,par1,par2,
     1       eps,tout,xout,yout,dxdtout,dydtout)
         implicit real *8 (a-h,o-z)
        save
         real *8 ts(1)
         external funcurve
c 
c         This subroutine finds the location (both on the
c         parametrizing interval [0,rl] and in R^2) of the
c         point specified by the user on the curve. The
c         curve has to have been preprocessed by a preceding
c         call to the subroutine anaresa (see), and the point
c         is specified by its distance (in terms of arc-length)
c         from the beginning of the curve. PLEASE NOTE THAT THIS
C         SUBROUTINE DEPENDS ON THE SUBROUTINE ANARESA FOR THE
C         PARAMETERS H,TS;  THIS SUBROUTINE HAAS NO FUNCTION AS
C         A STAND-ALONE DEVICE.
c 
c                        Input parameters:
c 
c  x - the location (in terms of arc-length) on the curve of the
c         point to be found
c  n - the number of nodes in the equispaced discretization of the
c         curve produced by a preceding call to anaresa
c  h - the sampling interval along the arc of the resampling returned
c          in array ts by the subroutine anaresa
c  ts - the equispaced discretization of the curve produced by the
c         subroutine anaresa
c  funcurve - the subroutine providing a parametrization the curve.
c        Should be the same as in the preceding call to anaresa. The
c        calling sequence of funcurve must be
c 
c                  funcurve(t,par1,par2,x,y,dxdt,dydt),
c 
c        with:
c 
c          t - the parameter by which the user's curve is parametrizaed
c          x - the x-coordinate of the point on the curve corresponding
c                to the values t of the parameter
c          y - the y-coordinate of the point on the curve corresponding
c                to the values t of the parameter
c          dxdt - the derivative with respect to t of the  x-coordinate
c                of the point on the curve corresponding to the values t
c                of the parameter
c          dydt - the derivative with respect to t of the  y-coordinate
c                of the point on the curve corresponding to the values t
c                of the parameter
c  eps - the precision to which the point is to be found. Generally, it
c        is a good idea to make this parameter the same as in the
c        preceding call to anaresa
c 
c                         Output parameters:
c 
c  tout - the value of the curve parameter corresponding to the location
c        x along the arc-length
c  xout, yout - the coordinates in R^2 of the point x on the curve
c  dxdt, dydt - the coordinates of the tangent (of length 1) to the
c        curve at the point (xout,yout)
c 
c 
c        . . . find the two nodes (in  arc-length) between which the
c              user-supplied point is located
c 
        ier=0
        d=x/h
        m=d
c 
c       if the point is outside the user-supplied curve - bomb
c 
        if(m .gt. n+1) ier=512
        if(m .gt. n+1) return
c 
c        use linear interpolation to find the approximate location
c        of the curve parameter (as opposed to the arc-length)
c        corresponding to the user-supplied arc-length
c 
        t1=ts(m+1)
        t2=ts(m+2)
        x1=m*h
        x2=x1+h
c 
        alpha=(t2-t1)/(x2-x1)
        beta=t1-alpha*x1
c 
        t0=alpha*x+beta
c 
c       use Newton to find the location of the curve parameter
c       (as opposed to the arc-length) corresponding to the
c       user-supplied arc-length
c 
        tk=t0
        tkold=t0
c 
        mm=12
        epsgauss=eps*1000
        epsnewt=dsqrt(eps)/1000
        maxit=20
        fk=h
c 
        do 1400 i=1,maxit
c 
c       evaluate the function whose root we are seeking
c 
        a=t1
        b=tk
c 
        step=1
        do 1200 j=1,20
c 
        call anaresga(jer,a,b,funcurve,par1,par2,mm,epsgauss,
     1      rltot,maxrec,numint)
c 
        fkold=fk
        fk=rltot+x1-x
c 
c       if the target function failed to decrease
c       - activate step-length control
c 
        if(abs(fk) .le. abs(fkold)) goto 1100
c 
        ier=1
        step=step/2
        tk=tkold-fkold/dldt *step
        goto 1200
c 
 1100 continue
c 
c       determine the dl/dt at that point and make a
c       newton step
c 
        call funcurve(tk,par1,par2,xk,yk,dxdt,dydt)
        dldt=dsqrt(dxdt**2+dydt**2)
c 
        tkold=tk
        tk=tk-fk/dldt *step
c 
c        if the required precision has been achieved
c        - terminate newton
c 
        if(abs(fk) .lt. epsnewt) goto 2000
c 
        goto 1400
 1200 continue
c 
        ier=64
        return
c 
 1400 continue
 2000 continue
c 
        tout=tk
        call funcurve(tout,par1,par2,xout,yout,dxdtout,dydtout)
c 
        d=sqrt(dxdtout**2+dydtout**2)
        dxdtout=dxdtout/d
        dydtout=dydtout/d
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine anaresga(ier,a,b,fun,par1,par2,m,eps,
     1      rint,maxrec,numint)
        implicit real *8 (a-h,o-z)
        save
        dimension t(100),w(100),stack(400),vals(200),
     1      par1(1),par2(1)
        data m7/-2341034/
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
c          ier=0 means normal conclusion
c          ier=8 means that at some point, one subinterval in the
c                subdivision was smaller than (b-a)/2**200. this
c                is a fatal error.
c          ier=16 means that the total number of subintervals in the
c                adaptive subdivision of [a,b] turned out to be greater
c                than 100000.  this is a fatal error.
c 
c  rint - the integral as evaluated
c  maxrec - the maximum depth to which the recursion went at its
c         deepest point. can not be greater than 200, since at that
c         point ier is set to 8 and the execution of the subroutine
c         terminated.
c  numint - the totla number of intervals in the subdivision. can not
c         be greater than 100000,  since at that
c         point ier is set to 16 and the execution of the subroutine
c         terminated.
c 
c 
c          . . . check if the gauss quadarture has been
c                initialized at a preceeding call to this routine
c 
        if(m .eq. m7) goto 1200
        call anagauss(m,t,w)
c 
        m7=m
 1200 continue
c 
c        integrate the user-supplied function using the
c        adaptive gaussian quadratures
c 
        nnmax=100000
        maxdepth=200
c 
        call anaresg1(ier,stack,a,b,fun,par1,par2,t,w,m,
     1      vals,nnmax,eps,rint,maxdepth,maxrec,numint)
        return
        end
c 
c 
c 
c 
c 
        subroutine anaresg1(ier,stack,a,b,fun,
     1      par1,par2,t,w,m,vals,nnmax,eps,
     2      rint,maxdepth,maxrec,numint)
        implicit real *8 (a-h,o-z)
        save
        dimension stack(2,1),t(1),w(1),vals(1),par1(1),par2(1)
c 
c       start the recursion
c 
        stack(1,1)=a
        stack(2,1)=b
        call anaresg2(a,b,fun,par1,par2,t,w,m,vals(1))
c 
c       recursively integrate the thing
c 
        j=1
        rint=0
        ier=0
        maxrec=0
        do 3000 i=1,nnmax
ccc        call prinf('i=*',i,1)
        numint=i
        if(j .gt. maxrec) maxrec=j
ccc        call prinf('j=*',j,1)
c 
c       subdivide the current subinterval
c 
         c=(stack(1,j)+stack(2,j))/2
        call anaresg2(stack(1,j),c,fun,
     1      par1,par2,t,w,m,value2)
c 
        call anaresg2(c,stack(2,j),fun,
     1      par1,par2,t,w,m,value3)
c 
        dd=dabs(value2+value3-vals(j))
cccc         call prin2('in anaresg1, dd=*',dd,1)
        ifdone=0
        if(dd .le. eps) ifdone=1
c 
c       if the function on this subinterval has been
c       integrated with sufficient accuracy - add the
c       value to that of the global integral and move up
c       in the stack
c 
        if(ifdone  .eq. 0) goto 2000
c 
        rint=rint+value2+value3
        j=j-1
c 
c        if the whole thing has been integrated - return
c 
        if(j .eq. 0) return
        goto 3000
 2000 continue
c 
c       if the function on this subinterval has not been
c       integrated with sufficient accuracy - move
c       down the stack
c 
        stack(1,j+1)=stack(1,j)
        stack(2,j+1)=(stack(1,j)+stack(2,j))/2
        vals(j+1)=value2
c 
        stack(1,j)=(stack(1,j)+stack(2,j))/2
        vals(j)=value3
c 
        j=j+1
c 
c       if the depth of the recursion has become excessive - bomb
c 
        if(j .le. maxdepth) goto 3000
        ier=8
        return
 3000 continue
        ier=16
        return
        end
c 
c 
c 
c 
c 
        subroutine anaresg2(a,b,fun,par1,par2,t,w,m,rint)
        implicit real *8 (a-h,o-z)
        save
        dimension t(1),w(1),par1(1),par2(1)
c 
c       integrate the function fun on the interval [a,b]
c 
        rint=0
        u=(b-a)/2
        v=(b+a)/2
        do 1200 i=1,m
        tt=u*t(i)+v
cccc         call prinf('in anaresg2 i=*',i,1)
ccc        call funcurve(tt,x,y,dxdt,dydt)
        call fun(tt,par1,par2,x,y,dxdt,dydt)
cccc         call prin2('in anaresg2 after funcurve, dxdt=*',dxdt,1)
        rint=rint+w(i)*dsqrt(dxdt**2+dydt**2)
cccccccc        rint=rint+fun(tt,par1,par2)*w(i)
 1200 continue
        rint=rint*u
        return
        end
c 
c 
c 
c 
c 
        subroutine anagauss(n,ts,whts)
        implicit real *8 (a-h,o-z)
        save
        dimension ts(1),whts(1)
c 
c        this subroutine constructs the nodes and the
c        weights of the n-point gaussian quadrature on
c        the interval [-1,1]
c 
c                input parameters:
c 
c  n - the number of nodes in the quadrature
c 
c                output parameters:
c 
c  ts - the nodes of the n-point gaussian quadrature
c  w - the weights of the n-point gaussian quadrature
c 
c       . . . construct the array of initial approximations
c             to the roots of the n-th legendre polynomial
c 
        eps=1.0d-14
        ZERO=0
        DONE=1
        pi=datan(done)*4
        h=pi/(2*n)
        do 1200 i=1,n
        t=(2*i-1)*h
        ts(n-i+1)=dcos(t)
1200  CONTINUE
c 
c         use newton to find all roots of the legendre polynomial
c 
        ts(n/2+1)=0
        do 2000 i=1,n/2
c 
        xk=ts(i)
        deltold=1
        do 1400 k=1,10
        call analegpo(xk,n,pol,der)
        delta=-pol/der
cccc         call prin2('delta=*',delta,1)
        xk=xk+delta
        if(deltold .le. eps) goto 1600
        deltold=dabs(delta)
 1400 continue
 1600 continue
        ts(i)=xk
        ts(n-i+1)=-xk
 2000 continue
c 
c       now, use the explicit integral formulae
c       to obtain the weights
c 
        a=-1
        b=1
        do 2200 i=1,n/2+1
        call anaprod(a,ts,n,i,fm)
        call anaprod(b,ts,n,i,fp)
        whts(i)=fp-fm
        whts(n-i+1)=whts(i)
 2200 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine analegpo(x,n,pol,der)
        implicit real *8 (a-h,o-z)
c 
        save
        pkm1=1
        pk=x
c 
        pk=1
        pkp1=x
c 
c        if n=0 or n=1 - exit
c 
        if(n .ge. 2) goto 1200
        pol=1
        der=0
        if(n .eq. 0) return
c 
        pol=x
        der=1
        return
 1200 continue
c 
c       n is greater than 1. conduct recursion
c 
        do 2000 k=1,n-1
        pkm1=pk
        pk=pkp1
        pkp1=( (2*k+1)*x*pk-k*pkm1 )/(k+1)
 2000 continue
c 
c        calculate the derivative
c 
        pol=pkp1
        der=n*(x*pkp1-pk)/(x**2-1)
        return
        end
c 
c 
c 
c 
c 
        subroutine anaprod(x,xs,n,i,f)
        implicit real *8 (a-h,o-z)
        save
        dimension xs(1)
c 
c      evaluate the product
c 
        f=1
        dlarge=1.0d20
        dsmall=f/dlarge
c 
        large=0
        do 2000 j=1,n
        dd=dabs(f)
        if( dd .gt. dsmall) goto 1200
         f=f*10000
         large=large-1
 1200 continue
c 
        if( dd .lt. dlarge) goto 1400
        f=f/10000
        large=large+1
 1400 continue
        if(j .eq. i) goto 2000
        f=f*(x-xs(j))/(xs(i)-xs(j))
 2000 continue
        d10000=10000
        f=f*d10000**large
        f=f**2*(x-xs(i))
        return
        end
