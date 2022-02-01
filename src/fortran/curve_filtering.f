

      subroutine simple_curve_resampler_mem(n,xy,nb,eps,nmax,nlarge,
     1   nout,lsave,lused,ier)
      implicit real *8 (a-h,o-z)
      integer n,nb
      real *8 xy(2,n)
      integer, intent(out) :: nlarge,nout

      real *8, allocatable :: x0(:),y0(:),work(:),ts(:)
      real *8, allocatable :: xver(:),yver(:)

      allocate(x0(n+1),y0(n+1))
      allocate(ts(n+1))
      allocate(xver(n+1),yver(n+1))

      do i=1,n
        x0(i) = xy(1,i)
        y0(i) = xy(2,i)
      enddo

      x0(n+1) = x0(1)
      y0(n+1) = y0(1)

      lused = 0
      lsave = 0
      ier = 0

      epscheck = max(eps,1.0d-12)
 
      do ii=1,nmax
        nlarge = 64*2**(ii)*n
        lenw = 10*nlarge*n + 10000
        allocate(work(lenw))
        curvelen = 0
        derr = 0
        nbl = 0
        iert = 0
        work = 0
        call rsblcurve(iert,x0,y0,n,nb,nlarge,
     1       curvelen,nbl,derr,ts,work,lenw,lsave,lused)
        nout = nlarge


        if(iert.ne.0) goto 1111
        erra = 0
        ra = 0
        do i=1,n
          xver(i) = 0
          yver(i) = 0
          call eval_curve(ier,ts(i),work,xver(i),yver(i),dxt,dyt,curv)
          erra = erra + (xver(i)-xy(1,i))**2
          ra = ra + xy(1,i)**2
          erra = erra + (yver(i)-xy(2,i))**2
          ra = ra + xy(2,i)**2
        enddo
        erra = sqrt(erra/ra)

        if(erra.lt.epscheck) goto 1000 

 1111   continue        
        deallocate(work)
      enddo

      if(iert.ne.0) ier = 4
      if(iert.eq.0) ier = 2
 1000 continue      


      return
      end
c
c
c
c
c
      subroutine simple_curve_resampler_guru(n,xy,nb,
     1   nlarge,lsave,lused,nout,srcinfoout,hout,curvelen,wsave,ts,ier)
      implicit real *8 (a-h,o-z)
      integer n,nout,lsave,lused
      real *8 xy(2,n),srcinfoout(6,nout),curvelen
      real *8 ts(n+1)
      real *8 wsave(lsave)
      real *8, allocatable :: x0(:),y0(:),work(:)

      allocate(x0(n+1),y0(n+1))

      do i=1,n
        x0(i) = xy(1,i)
        y0(i) = xy(2,i)
      enddo

      x0(n+1) = x0(1)
      y0(n+1) = y0(1)
      lenw = lused + 1000
      allocate(work(lenw))
      nbl = 0
      ier = 0

      
      call rsblcurve(ier,x0,y0,n,nb,nlarge,curvelen,nbl,
     1   derr,ts,work,lenw,lsave0,lused)

      
      wsave(1:lsave) = work(1:lsave)

      hout = curvelen/(nout+0.0d0)
      do i=1,nout
        t = (i-1)*hout
        call eval_curve(ier,t,wsave,srcinfoout(1,i),srcinfoout(2,i),
     1    dxt,dyt,srcinfoout(6,i))
        srcinfoout(5,i) = sqrt(dxt**2 + dyt**2)
        srcinfoout(3,i) = dyt/srcinfoout(5,i)
        srcinfoout(4,i) = -dxt/srcinfoout(5,i)
      enddo

      return
      end
c
c
c
c
c
      subroutine eval_curve_multi(n,ts,lsave,wsave,binfo)
      implicit real *8 (a-h,o-z)
      integer n,lsave
      real *8 ts(n),wsave(lsave),binfo(6,n)

      do i=1,n
        call eval_curve(ier,ts(i),wsave,binfo(1,i),binfo(2,i),dxt,dyt,
     1   binfo(6,i))
        binfo(5,i) = sqrt(dxt**2 + dyt**2)
        binfo(3,i) = dyt/binfo(5,i)
        binfo(4,i) = -dxt/binfo(5,i)
      enddo

      return
      end
c
c
c
c 
c
        subroutine rsblcurve(ier,x0,y0,n0,nmin,nlarge,
     1       curvelen,nbl,derr,ts,w,lenw,lsave,lused)
c
        implicit real *8 (a-h,o-z)
        real *8 x0(n0),y0(n0),w(*),ts(n0)
c 
c       this entry constructs a band-limit curve that passes through the
c       points of a user-supplied curve. the user supplies the curve in
c       the form of points (x0,y0) on the curve; the subroutine smooths
c       and resamplse the curve in an equispaced manner. in order to
c       evaluate points along the output curve, the user must call the
c       subroutine rsrespnt. this subroutine has no use as a stand-alone
c       device.
c
c                         input parameters:
c 
c  x0,y0 - user-supplied curve to be resampled
c  n0 - the number of nodes in (x0,y0). note that the input curve will
c       be treated as closed, i.e.,
c       x0(n0+1)=x0(1)
c       y0(n0+1)=y0(1)
c  nmin - desired bandlimit of the applied filter
c  nlarge - the number of nodes into which the curve is to be resampled;
c       must be large
c  lenw - the length (in real *8 words) of the array w provided by the
c       user
c 
c                   output parameters:
c 
c  ier - error return code;
c       ier=0  means normal conclusion
c       ier=4  means that the subroutine attempting to find the minimum
c              distance from some input point to the curve failed to
c              converge; this is a fatal error
c       ier=10 failure during call to anarsbl
c       ier=2000 output curve of rsresa intersects itself or nearly so;
c              this is a fatal error
c              possible remedy: increase nlarge and/or nmin
c       ier=2200 error when solving the linear system in constructing
c              the gaussians (this error may be unnecessary)
c       ier=16000 means that the length of the user-provided array
c              is insufficient; this is a fatal error
c  nbl - a measure of the bandlimit of the curvature of the resampled
c       curve; specifically, if w(1),w(2),...,w(nlarge) are the fourier
c       coefficients, then nbl is the integer such that        
c           \sum_{i=2}^nbl {abs(w(i))} \over 
c           \sum_{i=2}^nlarge {abs(w(i))}      ~    0.99
c  derr - average distance between the curve and the original data
c       (before computation of arc-length parametrization); in the
c       author's experience, this has only occurred when the input
c       data is sampled highly non-uniformly
c  ts - arc-length to the original data on output curve
c  w -  the array to be used by eval_curve (see below); the first
c       lsave elements of w should not be altered between the call to
c       this subroutine and the subsequent calls to eval_curve
c  lsave - the number of elements of w that should not be changed
c       between the call to this entry and the subsequent calls to
c       rsrespnt
c
        ier=0
c
c       allocate memory in the work array for the resampling routine
c
        iwsave=1
        lwsave=4*nlarge+30
c
        if (lwsave.le.lenw) goto 1100
        ier=16000
        return
 1100   continue
c
        call dcffti(nlarge,w(iwsave))

c
        iw=iwsave+lwsave+5
        lw=lenw-iw+1
c
c       resample curve at nlarge points, obtain tangent
c       angles, filter them, close the curve, and then
c       reconstruct the curve at nlarge points
c
        call rsblcurve0(ier,x0,y0,n0,nmin,nlarge,
     1       curvelen,err0,w(iw),lw,lsave1,lused1,w(iwsave))
        if (ier.ne.0) return
cccc    if (ier.ne.0.or.ier.ne.8.or.ier.ne.12) return
c
        lused=iw+lused1
c
        iw3=iw+lsave1+10
        lw3=7*nlarge+70
c
        its=iw3+lw3
        lts=nlarge+4
c 
        izs=its+lts
        lzs=2*nlarge+10
c 
        ider1=izs+lzs
        lder1=2*nlarge+10
c 
        ifis=ider1+lder1
        lfis=nlarge+4
c
        itn=ifis+lfis
        ltn=n0+4
c
        iww=itn+ltn
        lww=2*nlarge+10
c
        iw2=iww+lww
        lw2=lenw-iw2+1
c
c       result of rsblcurve0 is an upsampled closed curve
c       that does not go through the original data; use
c       gaussian perturbations along the curve to fix
c
        call rsblcurve1(ier,x0,y0,n0,nmin,nlarge,
     1       curvelen,nbl,derr,w(itn),w(its),w(izs),w(iw),w(iw3),
     2       w(iw2),lw2,lused2)
        if (ier.ne.0) return
c
        lused2=iw2+lused2
        if (lused2.gt.lused) lused=lused2
c
        ixgs=iw2
        lxgs=100
c
        iwhts=ixgs+lxgs
        lwhts=100
c
        iwr=iwhts+lwhts
        lwr=(lw2-300)/2
c
        iwv=iwr+lwr
        lwv=(lw2-300)/2
c
c       obtain arc-length parametrization
c
        call rsblcurve2(ier,n0,nlarge,curvelen,nbl,w(iw),w(iw3),
     1       m,w(ixgs),w(iwhts),w(iwr),w(iwv),lwr,lwv,nints,
     2       w(its),w(izs),w(ider1),w(ifis),
     3       h,eps,w(itn),ts,w(iww),w(iwsave))
        if (ier.ne.0) return
c
        lused2=iwv+nints+5
        if (lused2.gt.lused) lused=lused2
c
c       garbage collection...
c
        do 100 i=1,lsave1
        w(19+i)=w(iw+i-1)
 100    continue
        iw=20
        lsave=20+lsave1
c
        do 200 i=1,lw3
        w(lsave+9+i)=w(iw3+i-1)
 200    continue
        iw3=lsave+10
        lsave=lsave+10+lw3
c
        do 300 i=1,lts
        w(lsave+9+i)=w(its+i-1)
 300    continue
        its=lsave+10
        lsave=lsave+10+lts
c
        do 400 i=1,nints
        w(lsave+9+i)=w(iwr+i-1)
 400    continue
        iwr=lsave+10
        lsave=lsave+10+nints
c
        do 500 i=1,nints
        w(lsave+9+i)=w(iwv+i-1)
 500    continue
        iwv=lsave+10
        lsave=lsave+10+nints
c
        do 600 i=1,100
        w(lsave+9+i)=w(ixgs+i-1)
 600    continue
        ixgs=lsave+10
        lsave=lsave+10+100
c
        do 700 i=1,100
        w(lsave+9+i)=w(iwhts+i-1)
 700    continue
        iwhts=lsave+10
        lsave=lsave+10+100
c
        w(1)=nlarge+0.1
        w(2)=iw+0.1
        w(3)=iw3+0.1
        w(4)=its+0.1
        w(5)=iwr+0.1
        w(6)=iwv+0.1
        w(7)=nints+0.1
        w(8)=h
        w(9)=eps
        w(10)=m+0.1
        w(11)=ixgs
        w(12)=iwhts
c
        return
c
        end
c
c
c
c
c
        subroutine rsblcurve1(ier,x0,y0,n0,nmin,nlarge,
     1       curvelen,nbl,derr,tn,ts,zs,w,w3,w2,lenw2,lused2)
c
        implicit real *8 (a-h,o-z)
        real *8 x0(n0),y0(n0),tn(n0),ts(nlarge+1),zs(2,nlarge+1),
     1       w(*),w3(*),w2(*)
c
        iw2=1
        lw2=4*n0+10
c
        ixn=iw2+lw2
        lxn=n0+4
c
        iyn=ixn+lxn
        lyn=n0+4
c
        ia=iyn+lyn
        la=n0*n0
c
        iwd=ia+la
        lwd=2*n0+10
c
        irx=iwd+lwd
        lrx=n0+4
c
        iry=irx+lrx
        lry=n0+4
c
        numit=200
        iwk=iry+lry
        lwk=2*numit*n0+n0+60
c
        lused2=iwk+lwk
c
        if(lused2.le.lenw2) goto 1100
        ier=16000
        return
 1100   continue
c
c       find points tn on filtered curve that pass closest
c       to input data
c
        call getnodes(ier,x0,y0,n0,nmin,nlarge,
     1       curvelen,derr,w,w2(ixn),w2(iyn),tn,zs,ts)
        if (ier.ne.0) return
c
cccc        call prin2('and derr = *',derr,1)
c
c       add gaussians centered at tn to the filtered curve
c       so that it passes through the input data
c
        call pertgauss(ier,x0,y0,n0,nmin,nlarge,
     1       curvelen,nbl,derr,w,w3,w2,w2(ixn),w2(iyn),tn,
     2       w2(ia),w2(iwd),w2(irx),w2(iry),numit,w2(iwk))
c
        return
c
        end
c
c
c
c
c
        subroutine rsblcurve2(ier,n0,nlarge,curvelen,nbl,w,w3,
     1       mm,xgs,whts,wright,wvals,lwr,lwv,nints,ts,zs,der1,fis,
     2       h,eps,tn,tsout,ww,wsave)
c
        implicit real *8 (a-h,o-z)
        real *8 ts(nlarge+1),zs(2,nlarge+1),der1(2,nlarge+1),w(1),w3(1),
     1       xgs(1),whts(1),wright(1),wvals(1),tn(n0),tsout(n0),
     2       ww(1),wsave(1)
        real *8 curv
        external funcurv_fast        
c
c       subroutine to obtain arc-length discretization of gaussian
c       perturbed bandlimited curve (computed via funcurv_fast)
c
        rl=curvelen
        eps=1.0d-15
        mm=6
        call anarsbl(jer,funcurv_fast,w,w3,rl,nlarge,eps,ts,h,rltot,
     1       mm,xgs,whts,wright,wvals,lwr,lwv,nints)
c
        if (jer.eq.16000) then
        ier=jer
        return
        end if
c
        if (jer.ne.0.and.jer.ne.1) then
        call prinf('jer = *',jer,1)
        ier=10
        return
        end if
c
        tsout(1)=0
c
        ij0=1
c
        do 1200 i=2,n0
c
        t1=tn(i-1)*curvelen
        t2=tn(i)*curvelen
c
        call findint(wright,nints,t2,ij0,ijk)
c
        if (ij0.eq.ijk) then
        call anarsblg2(t1,t2,funcurv_fast,
     1       w,w3,xgs,whts,mm,tout)
        goto 1300
        end if
c
        call anarsblg2(t1,wright(ij0),funcurv_fast,
     1       w,w3,xgs,whts,mm,tout1)
c
        do 1500 kk=ij0+1,ijk-1
        tout1=tout1+wvals(kk)
 1500   continue
c
        call anarsblg2(wright(ijk-1),t2,funcurv_fast,
     1       w,w3,xgs,whts,mm,tout3)
c
        tout=tout1+tout3
c
 1300   continue
c
        tsout(i)=tsout(i-1)+tout
        ij0=ijk
c
 1200   continue
c
        do 1100 i=1,nlarge+1
        call funcurv_fast(ts(i),w,w3,zs(1,i),zs(2,i),
     1       der1(1,i),der1(2,i),curv)
 1100   continue
c
c       to account for change of variable to arc-length,
c       normalize the tangents
c
        do 1400 i=1,nlarge+1
        dtn=der1(1,i)**2+der1(2,i)**2
        dtn=sqrt(dtn)
        der1(1,i)=der1(1,i)/dtn
        der1(2,i)=der1(2,i)/dtn
 1400   continue
c
c       compute tangent angles and use that to compute bandlimit
c
        call rsfis(nlarge,der1,fis)
        call blapp(fis,nlarge,ww,wsave,nbl)
c
        curvelen=rltot
c
        return
c
        end
c
c
c
c
c
        subroutine getnodes(ier,xs,ys,n0,nmin,nlarge,
     1       curvelen,derr,w,xn,yn,tn,zs2,ts2)
c
        implicit real *8 (a-h,o-z)
        real *8 xs(n0),ys(n0),w(1),xn(n0),yn(n0),tn(n0),
     1       zs2(2,nlarge+1),ts2(nlarge+1),tang(10),curv(5),z3(10)
        real *8 der22(2)
c
c       subroutine that finds the points on the filtered curve
c       closest to the given data (xs,ys)
c
        h2=curvelen/nlarge
c
        do 100 i=1,nlarge+1
        ts2(i)=(i-1)*h2
 100    continue
c
        do 200 i=1,nlarge+1
        call rsrespnt(ts2(i),zs2(1,i),nlarge,tang,der22,w)
 200    continue
c
        derr=0
c
        del=1.0d-5*curvelen
cccc        del=1.0d-5
c
        do 300 i=1,n0
c
c       find the point on the curve
c       which is closest to xs,ys
c
c       to initialize, first find the closest
c       point among zs2; for the first point,
c       initialize at 0 to avoid situation where
c       it is initialized at approximately curvelen
c
        if (i.eq.1) then
        tt=0
        jk=1
        goto 900
        end if

        ddr=1.0d20
        do 400 j=1,nlarge
        dx=xs(i)-zs2(1,j)
        dy=ys(i)-zs2(2,j)
        dd0=dx*dx+dy*dy
        if (dd0.lt.ddr) then
        ddr=dd0
        tt=(j-1)*h2
        jk=j
        end if
 400    continue
c
 900    continue
c
c       find closest point via newton
c
        call getdist(jer,xs(i),ys(i),nlarge,w,tt,del,
     1       xn(i),yn(i),dn,tn(i))
c
        if (jer.eq.0) goto 1000
c
c       if newton does not converge to minimum
c       upsample region around minimum
c       and try again
c
        t1=(jk-2)*h2
        t2=jk*h2
c
        do 500 ii=1,3
c
        niter=10000
        do 600 j=1,niter+1
        t3=t1+(t2-t1)*(j-1)/niter
c
        call rsrespnt(t3,z3,nlarge,tang,der22,w)
c
        dx=xs(i)-z3(1)
        dy=ys(i)-z3(2)
        dd0=dx*dx+dy*dy
c
        if (dd0.lt.ddr) then
        ddr=dd0
        tt=t3
        end if
 600    continue
c
        call getdist(jer,xs(i),ys(i),nlarge,w,tt,del,
     1       xn(i),yn(i),dn,tn(i))
c
        if (jer.eq.0) goto 1000
c
        tt1=t1
        tt2=t2
        t1=tt-(tt2-tt1)/niter
        t2=tt+(tt2-tt1)/niter
c
 500    continue
c
c       newton failed to converge
c
        ier=4
        return
c
 1000   continue
c
        derr=derr+dn
c
 300    continue
c
        derr=derr/n0
cccc        call prin2_long('derr = *',derr,1)
c
c       normalize tn and verify it is in proper order
c
        tn(1)=tn(1)/curvelen
        do 800 i=2,n0
        tn(i)=tn(i)/curvelen
        if (tn(i).lt.tn(i-1)) ier=2000
 800    continue
c
        return
c
        end
c
c
c
c
c
        subroutine getdist(ier,xs,ys,nlarge,w,t2,del,xs2,ys2,dist,tout)
c
        implicit real *8 (a-h,o-z)
        real *8 w(1),tang(10),curv(5),zt(10)
        real *8 der22(2)

        real *8 witer(200 000),wdist(200 000)


c
c       subroutine to find the shortest distance between
c       xs,ys and the curve contained in w using newton
c
        ier=0
        ti=t2
c
c       use finite-difference newton
c       to find the shortest distance
c
cccc        eps=1.0d-10
        eps=1.0d-8
c
        do 100 i=1,50
c
        tt=ti
        call rsrespnt(tt,zt,nlarge,tang,der22,w)
        dx=zt(1)-xs
        dy=zt(2)-ys
        ff0=dx*dx+dy*dy
c
        tt=ti+del
        call rsrespnt(tt,zt,nlarge,tang,der22,w)
        dx=zt(1)-xs
        dy=zt(2)-ys
        ffu=dx*dx+dy*dy
c
        tt=ti-del
        call rsrespnt(tt,zt,nlarge,tang,der22,w)
        dx=zt(1)-xs
        dy=zt(2)-ys
        ffd=dx*dx+dy*dy
c
        ff1=(ffu-ffd)/(2*del)
cccc        call prin2('ff1 = *',ff1,1)
c
        if (abs(ff1).le.eps) goto 200
c
        ff2=(ffu-2*ff0+ffd)/(del*del)
cccc        call prin2('ff2 = *',ff2,1)
c
        ti=ti-ff1/ff2
cccc        call prin2('ti = *',ti,1)
c
 100    continue
c
c       newton failed to converge
c
        ier=1
        return
c
 200    continue
c
c       do three more iterations of newton
c       for good measure
c
        do 300 i=1,3
c
        tt=ti
        call rsrespnt(tt,zt,nlarge,tang,der22,w)
        dx=zt(1)-xs
        dy=zt(2)-ys
        ff0=dx*dx+dy*dy
c
        tt=ti+del
        call rsrespnt(tt,zt,nlarge,tang,der22,w)
        dx=zt(1)-xs
        dy=zt(2)-ys
        ffu=dx*dx+dy*dy
c
        tt=ti-del
        call rsrespnt(tt,zt,nlarge,tang,der22,w)
        dx=zt(1)-xs
        dy=zt(2)-ys
        ffd=dx*dx+dy*dy
c
        ff1=(ffu-ffd)/(2*del)
cccc        call prin2('ff1 = *',ff1,1)
c
        ff2=(ffu-2*ff0+ffd)/(del*del)
cccc        call prin2('ff2 = *',ff2,1)
c
        ti=ti-ff1/ff2
cccc        call prin2('ti = *',ti,1)
c
 300    continue
c
        tt=ti
        call rsrespnt(tt,zt,nlarge,tang,der22,w)
        dx=zt(1)-xs
        dy=zt(2)-ys
        ff0=dx*dx+dy*dy
c
c       outputs...
c
        xs2=zt(1)
        ys2=zt(2)
        dist=sqrt(ff0)
        tout=ti
c
c       verify minimum
c
        tt=ti+del
        call rsrespnt(tt,zt,nlarge,tang,der22,w)
        dx=zt(1)-xs
        dy=zt(2)-ys
        ffu=dx*dx+dy*dy
c
        tt=ti-del
        call rsrespnt(tt,zt,nlarge,tang,der22,w)
        dx=zt(1)-xs
        dy=zt(2)-ys
        ffd=dx*dx+dy*dy
c
        ff1=(ffu-ffd)/(2*del)
cccc        call prin2('final ff1 = *',ff1,1)
c
        ff2=(ffu-2*ff0+ffd)/(del*del)
cccc        call prin2('final ff2 = *',ff2,1)
c
c       ensure newton converged to a local maximum
c
        if (ff2.lt.0) ier=2
cccc        call prin2('ff2 = *',ff2,1)
c
        return
c
        end
c
c
c
c
c
        subroutine pertgauss(ier,xs,ys,n0,nmin,nlarge,
     1       curvelen,nbl,derr,w,w3,w2,xn,yn,tn,
     2       a,wdata,rhx,rhy,numit,work)
c
        implicit real *8 (a-h,o-z)
        real *8 xs(n0),ys(n0),w(*),w3(*),w2(*),
     1       xn(n0),yn(n0),tn(n0),
     2       a(n0*n0),wdata(2*n0),rhx(n0),rhy(n0),
     3       work(*),tang(10),curv(5),errs(1000)
        external funcurv,funcurv_fast,matvec
c
c       subroutine to add the appropriate gaussians to the bandlimited
c       curve to force it to pass through the user input data
c

c
c       memory management for w2 array
c
        w2(1)=n0
        w2(2)=nlarge
        w2(3)=curvelen
c
        it1=3
        it2=it1+n0
        it3=it2+n0
        it4=it3+n0
c
c       if magnitude of gaussian is 1 at the point
c       being shifted, the sharpness of the gaussian
c       is determined by requiring it be less/equal
c       than eps0 at both neighboring points
c
c       eps0 is chosen so that the resulting linear
c       system is diagonally dominant (by a factor rr)
c
        rr=2.0d0
        rf=1.0d0/rr
        eps1=1d-15
c
        do 300 j=1,n0
c
        if (j.eq.1) then
        t1=tn(n0)-1.0d0
        t2=tn(2)
        end if
c
        if (j.eq.n0) then
        t1=tn(n0-1)
        t2=tn(1)+1.0d0
        end if
c
        if (j.gt.1.and.j.lt.n0) then
        t1=tn(j-1)
        t2=tn(j+1)
        end if
c
        eps0=1.0d0
        do 400 k=1,20
c
        eps0=eps0/2.0d0
c
        bb=-log(eps0)
c
        b1=bb/((t1-tn(j))**2)
        b2=bb/((t2-tn(j))**2)
        b=b1
        if (b2.gt.b) b=b2
c
        sd=0
        do 500 i=1,n0
c
        dt=abs(tn(i)-tn(j))
        if (dt.gt.1.0d0) dt=dt-1.0d0
        if (dt.gt.0.5d0) dt=1.0d0-dt
c
        sd=sd+exp(-b*dt**2)
c
 500    continue
c
        sd=sd-1.0d0
c
        if (sd.lt.rf) goto 600
c
 400    continue
c
 600    continue
c
cccc        call prin2('eps0 = *',eps0,1)
c
c       ensure b is sufficiently large so that
c       exp(-b*0.5^2)<eps1**2; we do this because
c       diagonal dominance does not necessarily imply
c       that gaussians are sufficiently sharp, i.e., 
c       exp(-b*0.5)<eps1 -> stronger condition can speed
c       up computation but might negatively impact final result;
c       might want to look into this more...
c
c       note: in some cases, particularly when input data is not
c       roughly uniformly sampled, diagonal dominance is a bad
c       criteria - a more uniform choice of sigmas might be better;
c       need to investigate this further
c
        bb2=2*-log(eps1)/0.25d0
cccc        bb2=-log(eps1)/0.25d0
        if (bb2.gt.b) b=bb2
c
        w2(it1+j)=b
c
 300    continue
c
cccc        call prin2('sds = *',w2(it1+1),n0)
c
c       construct matrix of linear system
c       use sparse matrix data structure
c
        eps0=1.0d-15
        iind=1
c
c       initialize sparse data structure
c
        do 1500 i=1,2*n0
        wdata(i)=0
 1500   continue
c
        do 700 j=1,n0
c
        nsize=0
        iflag=0
c
        do 800 i=1,n0
c
        dt=abs(tn(i)-tn(j))
        if (dt.gt.1.0d0) dt=dt-1.0d0
        if (dt.gt.0.5d0) dt=1.0d0-dt
        dd=exp(-w2(it1+j)*dt**2)
c
c       determine the range of nonzero elements in
c       a particular row. if the range is [i1,i2], the
c       structure computes i2 and the number of elements
c       in the range (permitting wrapping as the matrix
c       is banded with nonzero corners)
c
        if (dd.ge.eps0) then
        a(iind)=dd
        iind=iind+1
        nsize=nsize+1
        iflag=1
        end if
c
        if (dd.lt.eps0.and.iflag.eq.1) then
        wdata(2*j-1)=i-1
        iflag=2
        end if
c
 800    continue
c
c       if the entire column is > eps0
c
        if (wdata(2*j-1).eq.0) wdata(2*j-1)=n0
c
c       size of column band
c
        wdata(2*j)=nsize
c
 700    continue
cccc        call prin2('wdata = *',wdata,2*n0)
c
c       construct right hand side
c
        do 1000 i=1,n0
        rhx(i)=xs(i)-xn(i)
        rhy(i)=ys(i)-yn(i)
 1000   continue
c
c       solve the linear system
c
        eps1=1.0d-15
        call dumb_conres(matvec,a,wdata,rhx,n0,eps1,numit,
     1       w2(it2+1),errs,niter,work)
cccc        call prin2('errsx = *',errs,niter)
        call dumb_conres(matvec,a,wdata,rhy,n0,eps1,numit,
     1       w2(it3+1),errs,niter,work)
cccc        call prin2('errsy = *',errs,niter)
c
cccc        call prin2('x soln = *',w2(it2+1),n0)
cccc        call prin2('y soln = *',w2(it3+1),n0)
c
c       last bit of info needed for funcurv
c
        do 1300 i=1,n0
        w2(it4+i)=tn(i)
 1300   continue
c
c       verify that the curve goes through the
c       original data
c
        dd=0
c
        do 1200 i=1,n0
c
        tt=tn(i)*curvelen
        call funcurv(tt,w,w2,xt,yt,dxt,dyt)
c
        dx=xs(i)-xt
        dy=ys(i)-yt
c
        dd=dd+sqrt(dx*dx+dy*dy)
c
 1200   continue
c
        dd=dd/n0
c
cccc        call prin2('after gaussians, l1 diff = *',dd,1)
c
        eps2=eps1*100
        if (dd.gt.eps2) then
        ier=2200
        return
        end if
c
c       resample the new curve at equally spaced
c       points using anarsbl
c
        w3(1)=n0
        w3(2)=nlarge
        w3(3)=curvelen
        w3(4)=tn(1)
c
        its=5
        lts=nlarge+4
c 
        izs=its+lts
        lzs=2*nlarge+10
c 
        ider1=izs+lzs
        lder1=2*nlarge+10
c 
        ider2=ider1+lder1
        lder2=2*nlarge+10
c
        ncoefs=20
        w3(ider2+lder2)=ncoefs+0.1
        icoefs=ider2+lder2+1
        lcoefs=ncoefs
c
        h=curvelen/nlarge
        do 2400 i=1,nlarge+1
        w3(its+i-1)=(i-1)*h
 2400   continue
c
        call rslagrin(w3(its),ncoefs,w3(icoefs))
c
        do 2500 i=1,nlarge+1
        call fungauss(w3(its+i-1),w2,w3(izs+2*i-2),w3(izs+2*i-1),
     1       w3(ider1+2*i-2),w3(ider1+2*i-1),w3(ider2+2*i-2),
     2       w3(ider2+2*i-1))
 2500   continue
c
c       shift the parametrization, t, so that the curve
c       at t=0 gives xs(1),ys(1); note that this shift
c       is accounted for in funcurv_fast, but not in
c       funcurv
c
        tn1=tn(1)
        do 1600 i=1,n0
        tn(i)=tn(i)-tn1
 1600   continue
c
        dd=0
        do 1400 i=1,n0
        tt=tn(i)*curvelen
        call funcurv_fast(tt,w,w3,xt,yt,dxt,dyt,curv)
        dx=xs(i)-xt
        dy=ys(i)-yt
        dd=dd+sqrt(dx*dx+dy*dy)
 1400   continue
        derr=dd/n0
cccc        call prin2('with interpolation, l1 diff = *',dd,1)
c
        return
c
        end
c
c
c
c
c
        subroutine eval_curve(ier,t,w,x,y,dxdt,dydt,curv)
c
        implicit real *8 (a-h,o-z)
        real *8 w(1),w2(1)
        external funcurv_fast
c
        nlarge=w(1)
        iw=w(2)
        iw2=w(3)
        its=w(4)
        iwr=w(5)
        iwv=w(6)
        nints=w(7)
        h=w(8)
        eps=w(9)
        m=w(10)
        ixgs=w(11)
        iwhts=w(12)
        curvelen = w(iw2+3-1)
c
        call anaptbl(ier,t,nlarge,h,w(its),funcurv_fast,w(iw),w(iw2),
     1       eps,m,w(ixgs),w(iwhts),w(iwr),w(iwv),nints,curvelen,tout,
     2       x,y,dxdt,dydt,curv)
c
        return
c
        end
c
c
c
c
c
        subroutine funcurv_fast(t,w,w2,x,y,dxdt,dydt,curv)
c
        implicit real *8 (a-h,o-z)
        real *8 w(1),w2(1),tang(10),curv,z(10),curv0,der22(2)
        real *8 d2xdt2,d2ydt2
c
        n0=w2(1)
        nlarge=w2(2)
        curvelen=w2(3)
        tn1=w2(4)
        iw=5
c
c       undo the shift introduced in pertgauss that
c       ensures that t=0 corresponds to the x0(1),y0(1)
c       where x0,y0 is the user input data to rsblcurve
c
        t2=t+tn1*curvelen
c
c       evaluate the curve
c
        curv = 0
        call rsrespnt(t2,z,nlarge,tang,der22,w)
c
        x=z(1)
        y=z(2)
        dxdt=tang(1)
        dydt=tang(2)

        dtn = dxdt**2 + dydt**2
        dtn = sqrt(dtn)
        dxdt = dxdt/dtn
        dydt = dydt/dtn
        d2xdt2 = der22(1)
        d2ydt2 = der22(2)

c
c       evaluate the guassian perturbations
c
        call rsrespnt(t2,z,nlarge,tang,der22,w2(iw))


        x=x+z(1)
        y=y+z(2)
        dxdt=dxdt+tang(1)
        dydt=dydt+tang(2)

        d2xdt2 = d2xdt2 + der22(1)
        d2ydt2 = d2ydt2 + der22(2)

        rr = sqrt(dxdt**2 + dydt**2)
        curv = 0 
c
        return
c
        end
c
c
c
c
c
        subroutine funcurv(t,w,w2,x,y,dxdt,dydt)
c
        implicit real*8 (a-h,o-z)
        real *8 w(1), w2(1),
     1       tang(10), curv(5), z(10),der22(2)
c
c       subroutine to evaluate the curve which is then adjusted
c       by gaussians
c
        n0=w2(1)
        nlarge=w2(2)
        curvelen=w2(3)
c
c
c       evaluate the x,y point along the curve
c       as well as the derivatives dx/dt and dy/dt
c
        call rsrespnt(t,z,nlarge,tang,der22,w)
c
        x=z(1)
        y=z(2)
        dxdt=tang(1)
        dydt=tang(2)
c
c       now add the gaussians
c
        call fungauss(t,w2,x2,y2,dxdt2,dydt2,d2xdt2,d2ydt2)
c
        x=x+x2
        y=y+y2
        dxdt=dxdt+dxdt2
        dydt=dydt+dydt2
c
        return
c
        end
c
c
c
c
c
        subroutine fungauss(t,w,x,y,dxdt,dydt,d2xdt2,d2ydt2)
c
        implicit real *8 (a-h,o-z)
        real *8 w(1)
c
        n0=w(1)
        nlarge=w(2)
        curvelen=w(3)
c
        it1=4
        lt1=n0
c
        it2=it1+lt1
        lt2=n0
c
        it3=it2+lt2
        lt3=n0
c
        it4=it3+lt3
        lt4=n0
c
        t2=t/curvelen
        tf=w(it4)
        tl=1.0d0+tf
c
        x=0
        y=0
        dxdt=0
        dydt=0
        d2xdt2 = 0
        d2ydt2 = 0
c
        do 100 j=1,n0
c
        t3=t2
c
        tn=w(it4+j-1)
        dt2=t3-tn
        dt=abs(dt2)
c
        if (dt.gt.1.0d0) dt=dt-1.0d0
        if (dt.gt.0.5d0) dt=1.0d0-dt
c
        a1=w(it2+j-1)
        a2=w(it3+j-1)
        sig=w(it1+j-1)
        gau=exp(-sig*dt**2)
c
        dx=a1*gau
        dy=a2*gau
c
c       sign of dt matters...
c
        if (dt2.gt.0.5d0) sgn=-1
        if (dt2.le.0.5d0.and.dt2.ge.0) sgn=1
        if (dt2.lt.0.and.dt2.gt.-0.5d0) sgn=-1
        if (dt2.le.-0.5d0) sgn=1
c
        ddx=-2*a1*sig*gau*dt/curvelen*sgn
        ddy=-2*a2*sig*gau*dt/curvelen*sgn

        ddx2 = -2*a1*sig*gau/curvelen**2 + 
     1     4*a1*gau*(sig*dt/curvelen)**2
        ddy2 = -2*a2*sig*gau/curvelen**2 + 
     1     4*a2*gau*(sig*dt/curvelen)**2
        

c
        x=x+dx
        y=y+dy
c
        dxdt=dxdt+ddx
        dydt=dydt+ddy

        d2xdt2 = d2xdt2 + ddx2
        d2ydt2 = d2ydt2 + ddy2
c
 100    continue
c
        return
c
        end
c
c
c
c
c
        subroutine matvec(a,n,x,y,ww)
c
        implicit real *8 (a-h,o-z)
        real *8 a(1), x(n), y(n), ww(1)
c
        iind=1
c
        do 500 j=1,n
        y(j)=0
 500    continue
c
        do 100 j=1,n
c
        itot=ww(2*j)
        iend=ww(2*j-1)
        ist=iend-itot+1
c
c       if ist>0 then non-zero index range for a column is consecutive,
c       i.e. [ist,iend]; otherwise range wraps around the end of the
c       column (e.g. if column has 12 elements with bandwidth 5, then
c       [ist,iend]=[-2,2] which is [1,2] and [10,12] with standard
c       indexing)
c
        if (ist.gt.0) then
        is0=ist
c
c       skip 300 loop if range is consecutive
c
        ist=n+1
        goto 400
        end if

        is0=1
        ist=ist+n

 400    continue
c
        do 200 i=is0,iend
c
        y(i)=y(i)+a(iind)*x(j)
        iind=iind+1
c
 200    continue
c
        do 300 i=ist,n
c
        y(i)=y(i)+a(iind)*x(j)
        iind=iind+1
c
 300    continue
c
 100    continue
c
        return
c
        end
c
c
c
c
c
        subroutine matvec22(a,n,x,y,w)
c
        implicit real *8 (a-h,o-z)
        real *8 a(n,n),x(n),y(n),w(1)
c
        do 100 i=1,n
        dd=0
        do 200 j=1,n
        dd=dd+a(i,j)*x(j)
 200    continue
        y(i)=dd
 100    continue
c
        return
c
        end
c
c
c
c
c 
        subroutine rsblcurve0(ier,x0,y0,n0,nmin,nlarge,
     1       curvelen,err,w,lenw,lsave,lused,wsave)
        implicit real *8 (a-h,o-z)
        real *8 x0(n0),y0(n0),w(*),wsave(*),
     1       z(2),tang(2),z2(2),tang2(2),der22(2)
c 
c       this entry constructs an equispaced discretization of the
c       user-supplied curve, to be used to obtain arbitrary points
c       on that curve, specified by their arc-length distance from
c       the origin. the user provides the curve in the form of the
c       nodes (x0,y0) on the curve; the subroutine smoothes and
c       resamples the curve in an equispaced manner; in reality, it
c       feeds the subroutine rsresa (see elsewhere) that
c       resamples the curve. please note that the actual evaluations
c       are performed by the entry rsrespnt of this subroutine (see
c       below). this entry has no use as a stand-alone device.
c 
c                   input parameters:
c 
c  x0,y0 - user-supplied curve to be resampled
c  n0 - the number of nodes in (x0,y0). note that
c       the input curve will be treated as closed, i.e.
c       x0(n0+1)=x0(1)
c       y0(n0+1)=y0(1)
c  nmin - the highest order of a fourier node of the
c       of the user's curve with penalty coefficient one in the
c       supergain
c  nlarge - the number of nodes into which the curve is to be resampled;
c       must be large
c  lenw - the length (in real *8 words) of the array w provided by
c       the user
c 
c                   output parameters:
c 
c  ier - error return code;
c       ier=0  means normal conclusion
c       ier=8  means that the newton process closing
c              the curve failed to converge to the prescribed
c              precision. more specifically, the discrepancy
c              failed to go down after an iteration. this is
c              not a fatal error, but a careful examination
c              of the vector errs (see below) is recommended
c       ier=12 means that the newton process closing
c              the curve failed to converge to the prescribed
c              precision. more specifically, the discrepancy
c              failed to decay below eps after 30 iterations,
c              while it kept going down on every step. this is
c              not a fatal error, but a careful examination
c              of the vector errs (see below) is recommended.
c              furhthermore, this is an extremely suspicious
c              situation that has been never encountered in
c              our experiments
c       ier=32 means that n is less than or equal to nmin*4.
c              this is a fatal error
c       ier=64 means that n is less than or equal to n0
c              this is a fatal error
c       ier=96 means that both of the preceeding conditions
c              have occured; this is a fatal error
c       ier=16000 means that the length of the user-provided array
c              is insufficient. this is a fatal error
c  curvelen - the length of the filtered and resampled curve
c  err - the accuracy to which (in the opinion of the subroutine) the
c        curve could be resampled by nlarge7 nodes with the parameter
c        nmin as supplied by the user.
c  w - the array to be used by the entry rsrespnt (see below);
c         the first lsave elements of w should not be altered between
c         the call to this entry and the subsequent calls to rsrespnt
c  lsave - the number of elements of w that should not be changed
c         between the call to this entry and the subsequent calls
c         to rsrespnt
c 
c       . . . allocate memory for the resampling of the user-supplied
c             curve at an ungodly number of points
c
        its=1
        lts=nlarge+4
c 
        izs=its+lts
        lzs=2*nlarge+10
c 
        ider1=izs+lzs
        lder1=2*nlarge+10
c 
        ider2=ider1+lder1
        lder2=2*nlarge+10
c
        ncoefs=20
        w(ider2+lder2)=ncoefs+0.1
        icoefs=ider2+lder2+1
        lcoefs=ncoefs
c 
        ifis=icoefs+lcoefs
        lfis=nlarge+4
c
        lsave=ifis-1
c
        lused=ifis+lfis+5 
c
        if(lused.le.lenw) goto 1100
        ier=16000
        return
 1100   continue
c
        iw=lused+1
        lenw2=lenw-lused
c 
c       resample the user-specified curve at an ungodly
c       number of points
c 
        eps=1.0d-14
        nders=2

        h = 0
        acc = 0
        err = 0




c
        call rsresa(ier,x0,y0,n0,nmin,nders,
     1     w(izs),w(ider1),w(ider2),nlarge,w(ifis),h,
     2     acc,err,w(iw),lenw2,lused2,wsave)
        if (ier.ne.0) return
cccc        if (ier.ne.0.and.ier.ne.8.and.ier.ne.12) return
c
        lused=lused+lused2
c 
cccc        call prin2('after rsresa, acc=*',acc,1)
cccc        call prin2('after rsresa, err=*',err,1)
c
        err=err+acc
c 
        do 1200 i=1,nlarge+1
        w(its+i-1)=(i-1)*h 
 1200   continue
c
        call rslagrin(w(its),ncoefs,w(icoefs))
        curvelen=nlarge*h
c
        return
c
        end
c
c 
c 
c 
c 
        subroutine rsrespnt(t,z,nlarge,tang,der22,w)
c
        implicit real *8 (a-h,o-z)
        real *8 z(2),tang(2),w(1),der22(2)
c 
c       this entry finds the location of a point on the curve,
c       given the user-supplied distance of this point from the
c       beginning of the curve, along the arc-length. in order
c       for this entry to work, the user must have called the
c       rsblcurve (see above).
c 
c                    input parameters:
c 
c  t - the distance of the point from the beginning of the curve,
c       along the arc-length
c 
c  w - array produced by a prior call to the entry rsblcurve0 (see above)
c 
c                     output parameters:
c 
c  z - the location of the point in r^2
c  tang - the tangent vector at the point z (not normalized)
c  curv - the curvature of the curve at the point t. please note that
c          curv can be either positive or negative; it is positive when
c          the curve is convex

c 
c        use interpolation to find the value of the original parameter
c        corresponding to the arclength value t
c
        its=1
        lts=nlarge+4
c 
        izs=its+lts
        lzs=2*nlarge+10
c 
        ider1=izs+lzs
        lder1=2*nlarge+10
c 
        ider2=ider1+lder1
        lder2=2*nlarge+10
c
        ncoefs=w(ider2+lder2)
        icoefs=ider2+lder2+1
        lcoefs=ncoefs
c
        call rsresper(w(its),w(izs),w(ider1),w(ider2),
     1      nlarge,t,w(icoefs),ncoefs,z,tang,der22)
c 
  
c 
        return
        end
  
c 
c 
c 
c 
c 
        subroutine rsresper(ts,zs,tangs,der2,m,x,coefs,n,zout,tangout,
     1      der2out)
        implicit real *8 (a-h,o-z)
        save
        real *8 ts(*),coefs(*),coefsx(100),zs(2,*),tangs(2,*),
     1      tsnew(100),tangsnew(2,100),zout(2),tangout(2),
     2      zsnew(2,100),der2(2,*),der2out(2),der2new(2,100)
c 
c       this subroutine interpolates (via standard lagrange
c       interpolation) a user-specified periodic function
c       from the user-provided equispaced grid to the user-supplied
c       point x on the line. the point is permitted to be outside the
c       interval where the discretization is located, but by no more
c       than the length of the interval on which the function is
c       originally defined (i.e. by no more than ts(m)-ts(1)). please
c       note that the subpoutine assumes that ts(m) .neq. ts(1),
c       just like standard fft routines do; conceptually,
c       ts(m+1)=ts(1).
c 
c 
c                    input parameters:
c 
c  ts - the equispaced nodes at which the user-supplied function is
c       tabulated
c  zs - the values of the interpolant at the nodes ts
c  m - the number of nodes in the array ts
c  x - the point where the function is to be interpolated
c  ifinit - an integer parameter telling the subroutine whether the
c       lagrange interpolation routine has to be initialized.
c          ifinit=1 will cause the initialization to be performed
c          ifinit=0 will cause the initialization to be skipped.
c       please note that ifinit has to be set to one each time the
c       sampling distance in the array ts changes
c 
c                     output parameters:
c 
c  zout - the interpolated value of the function at the point x
c 

c 
c       use interpolation to find the value of the function
c       for the user-specified x
c 
c       . . . find on which subinterval x lives
c 
        ier=0
c
        h=ts(2)-ts(1)
        iint=(x-ts(1))/h
        istart=iint-n/2+2
        iend=istart+n-1
c 
       do 2300 i=1,n
c 
       tsnew(i)=ts(1)+(istart+i-2)*h
cccc       fsnew(i)=fs(irsrwrap(m,istart+i-1))
       zsnew(1,i)=zs(1,irsrwrap(m,istart+i-1))
       zsnew(2,i)=zs(2,irsrwrap(m,istart+i-1))
c 
       tangsnew(1,i)=tangs(1,irsrwrap(m,istart+i-1))
       tangsnew(2,i)=tangs(2,irsrwrap(m,istart+i-1))
c 
       der2new(1,i)=der2(1,irsrwrap(m,istart+i-1))
       der2new(2,i)=der2(2,irsrwrap(m,istart+i-1))
c 
 2300 continue
c 
cccc        call prin2('in rsresper, x=*',x,1)
c
       call rslagrco(tsnew,coefs,n,x,coefsx)
c 
cccc        call prin2('in rsresper, coefsx=*',coefsx,n)
cccc        call prin2('in rsresper, coefs=*',coefs,n)
c 
        zout(1)=0
        zout(2)=0
        tangout(1)=0
        tangout(2)=0
  
        der2out(1)=0
        der2out(2)=0
c 
        do 2400 i=1,n
        zout(1)=zout(1)+coefsx(i)*zsnew(1,i)
        zout(2)=zout(2)+coefsx(i)*zsnew(2,i)
c 
        tangout(1)=tangout(1)+coefsx(i)*tangsnew(1,i)
        tangout(2)=tangout(2)+coefsx(i)*tangsnew(2,i)
c 
        der2out(1)=der2out(1)+coefsx(i)*der2new(1,i)
        der2out(2)=der2out(2)+coefsx(i)*der2new(2,i)
c 
 2400 continue
c
        return
        end
c 
c 
c 
c 
c 
        function irsrwrap(m,i)
c 
        if( (i .le. m) .and. (i .ge. 1)) then
           irsrwrap=i
           return
        endif
c 
        if(i .lt. 1) irsrwrap=m+i
        if(i .gt. m) irsrwrap=i-m
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine rslagrin(t,n,coefs)
c
        implicit real *8 (a-h,o-z)
        real *8 t(1),coefs(1)
c 
c       this subroutine performs initialization for the entry
c       rslagrco (see below). its output is the array coefs
c       to be used by rslagrco to construct interpolation
c       coefficients for the interpolation on the real line.
c 
c                         input parameters:
c 
c  t - the nodes on which the interpolation is based.
c  n - the number of elements in array t
c 
c                         output parameters:
c 
c  coefs -  an array of length n to be used by the entry lagreva
c 
c        . . . initialize the interpolation routine based on
c              user-specified nodes t(i)
c 
        done=1
        do 1400 i=1,n
        cd=1
        do 1200 j=1,n
        if(j.eq.i) goto 1200
        cd=cd*(t(i)-t(j))
 1200   continue
        coefs(i)=done/cd
 1400   continue
c
        return
c
        end
c 
c 
c 
c 
        subroutine rslagrco(t,coefs,n,x,coefsx)
c
        implicit real *8 (a-h,o-z)
        real *8 t(1),coefs(1),coefsx(1)
c 
c       this entry constructs interpolation coiefficients
c       connecting the values of a function f tabulated at the
c       nodes t(i) on the real line  with the value of f at
c       the real point x
c 
c                       input parameters:
c 
c  t - the nodes on which the interpolation is based.
c  coefs -  an array of length n produced the entry rslagrin (see)
c  n - the number of elements in arrays t, coefs
c  x - the point in the complex plane where f is to be evaluated
c 
c                         output parameters:
c 
c  coefsx - the array of interpolation coefficients
c 
c       . . . find the node nearest to this point
c 
        knear=1
        cd=x-t(1)
        dd=cd**2
        dist=dd
        do 2400 i=2,n
        cd=x-t(i)
        dd=cd**2
        if(dd .gt. dist) goto 2400
        dist=dd
        knear=i
 2400 continue
c 
c       calculate the coefficient in front of the sum
c 
        coesum=1
        do 4200 i=1,n
        coesum=coesum*(x-t(i))
 4200 continue
c 
c       calculate the sum excluding the term corresponding
c       to knear
c 
        do 4400 i=1,n
        if(i .eq. knear) goto 4400
c 
        coefsx(i)=coefs(i)*coesum/(x-t(i))
 4400 continue
c 
c       account for the term knear
c 
        snear=1
        do 4600 i=1,n
        if(i .eq. knear) goto 4600
        snear=snear*(x-t(i))
 4600 continue
        coefsx(knear)=snear*coefs(knear)
c 
        return
        end
c 
c 
c 
c 
c 
  
        subroutine rsresa(ier,x0,y0,n0,nmin,nders,
     1     z,der1,der2,n,fis,h,acc,err,w,lenw,lused,wsave)
        implicit real *8 (a-h,o-z)
        save
        dimension x0(*),y0(*),z(2,*),der1(2,*),der2(2,*),
     1       fis(*),w(*),wsave(*)
c 
c        this subroutine resamples a user-specified curve in an
c        equispaced manner with respect to the wavelength.
c        in addition, it produces several derivatives (1, 2, or 3)
c        of the curve with respect to the arc-length, tabulated
c        at the nodes of the discretization. as a matter of fact,
c        this is a memory management routine. all actual
c        processing is performed by the routine rsres2 (see).
c 
c          important note:
c 
c      the curve to be resampled must be closed.
c 
c                      input parameters:
c 
c  x0,y0 - user-supplied curve to be resampled
c  n0 - the number of nodes in (x0,y0). note that
c       the input curve will be treated as closed, i.e.
c       x0(n0+1)=x0(1)
c       y0(n0+1)=y0(1)
c   nmin - the highest order of a fourier node of the
c        curvature of the user's curve left unchanged by the
c        filter.
c   nders - the number of derivatives (with respect to
c        arc length) of the resampled curve to be produced
c        by the subroutine. permitted values: 1, 2, 3.
c        the values of the parameters that were not requested
c        are not returned, and these parameters can be
c        replaced by dummies.
c  n - the number of nodes into which the curve is to be resampled
c 
c                          output parameters:
c 
c   ier -error return code.
c        ier=0 means succcessful conclusion.
c        ier=8 means that the newton process closing
c              the curve failed to converge to the prescribed
c              precision. more specifically, the discrepancy
c              failed to go down after in iteration. this is
c              not a fatal error, but a careful examination
c              of the vector errs (see below) is recommended.
c        ier=12 means that the newton process closing
c              the curve failed to converge to the prescribed
c              precision. more specifically, the discrepancy
c              failed to decay below eps after 30 iterations,
c              while it kept going down on every step. this is
c              not a fatal error, but a careful examination
c              of the vector errs (see below) is recommended.
c              furhthermore, this is an extremely suspicious
c              situation that has been never encountered in
c              our experiments.
c        ier=32 means that n is less than or equal to nmin*4.
c              this is a fatal error
c        ier=64 means that n is less than or equal to n0.
c              this is a fatal error
c        ier=96 means that both of the preceeding conditions
c               have occured.
c              this is a fatal error
c 
c   z - resampled curve
c   der1 - the first derivative of the resampled curve
c   der2 - the second derivative of the resampled curve
c   fis -tangent angles for the smoothed curve at the
c           resampled nodes
c   h - the sampling interval on the resampled curve
c   acc -  the precision with which the curve has been resampled,
c           in the sense that acc is the estimated value of the
c           fourier series of z (as a function of arclength
c           of the curve) near n/2
c 
c   err - the accuracy with which the newton process
c           managed to close the curve
c 
c                       work arrays:
c 
c   w - must be at least 8*n+n0+100 real *8 elements long
c 
c       . . . verify that the user-supplied integer parameters
c             are compatible with each other
c
        ier=0
        if(n.le.n0) ier=64
        if(n.le.nmin*2) ier=ier+32
        if(ier.ne.0) return
c 
c        allocate memory for the resampling routine
c 
        is=1
        ls=n0+4
c 
        iw=is+ls
        lw=2*n+10
c 
        ierrs=iw+lw
        lerrs=26
c 
        iz2=ierrs+lerrs
        lz2=2*n+10
c 
        lused=iz2+lz2+5
c
        if(lused.le.lenw) goto 1100
        ier=16000
        return
 1100   continue
c 
c       set the appropriate parameters
c 
        delta=1.0d-5
        eps=1.0d-9

c 
c       . . . resample
c
        call rsres2(ier,x0,y0,n0,w(is),fis,
     1    n,w(iw),wsave,nmin,z,w(iz2),eps,niter,
     2    w(ierrs),delta,acc,der1,der2,nders,h)
         err=w(ierrs+niter-1)
cccc         call prin2('after rsres2, acc=*',acc,1)
cccc         call prin2('after rsres2, h=*',h,1)
cccc         call prin2('after rsres2, errs=*',w(ierrs),niter)
cccc         call prin2('after rsres2, err=*',err,1)
c
         return
         end
c 
c 
c 
c 
c 
        subroutine rsres2(ier,x0,y0,n0,s,fis,
     1    n,w,wsave,nmin,z,z2,eps,niter,errs,delta,acc,
     2    der1,der2,nders,h)
        implicit real *8 (a-h,o-z)
        real *8 x0(n0),y0(n0),s(n0+3),fis(*),w(*),wsave(*),
     1     z(2,*),z2(2,*),errs(*),der1(2,*),der2(2,*)
c 
c        this subroutine resamples a user-specified curve in an
c        equispaced manner with respect to the wavelength.
c        in addition, it produces several derivatives (1, 2, or 3)
c        of the curve with respect to the arc-length, tabulated
c        at the nodes of the discretization.
c 
c                      input parameters:
c 
c  x0,y0 - user-supplied curve to be resampled
c  n0 - the number of nodes in (x0,y0). note that
c       the input curve will be treated as closed, i.e.
c       x0(n0+1)=x0(1)
c       y0(n0+1)=y0(1)
c  n - the number of nodes into which the curve is to be resampled
c   nmin - the highest order of a fourier node of the
c        curvature of the user's curve left unchanged by the
c        filter.
c   nders - the number of derivatives (with respect to
c        arc length) of the resampled curve to be produced
c        by the subroutine. permitted values: 1, 2, 3.
c        the values of the parameters that were not requested
c        are not returned, and these parameters can be
c        replaced by dummies.
c  eps - precision to which the newton iterations will be
c        conducted in order to close the curve
c  delta - the step used in the evaluation of derivatives
c        via the finite differences during the newton
c 
c                          output parameters:
c 
c   ier -error return code.
c        ier=0 means succcessful conclusion.
c        ier=8 means that the newton process closing
c              the curve failed to converge to the prescribed
c              precision. more specifically, the discrepancy
c              failed to go down after an iteration. this is
c              not a fatal error, but a careful examination
c              of the vector errs (see below) is recommended.
c        ier=12 means that the newton process closing
c              the curve failed to converge to the prescribed
c              precision. more specifically, the discrepancy
c              failed to decay below eps after 25 iterations,
c              while it kept going down on every step. this is
c              not a fatal error, but a careful examination
c              of the vector errs (see below) is recommended.
c              furhthermore, this is an extremely suspicious
c              situation that has been never encountered in
c              our experiments.
c   fis -tangent angles for the smoothed curve at the
c           resampled nodes
c   z - resampled curve
c   niter - the number of iterations taken by the newton
c        process to close the curve
c   errs - the array of errors produced by newton on its
c        niter iterations
c   der1 - the first derivative of the resampled curve
c   der2 - the second derivative of the resampled curve
c   acc -  the precision with which the curve has been resampled
c   h - the sampling interval on the resampled curve
c 
c                       work arrays:
c 
c   s - must be at least n0+3 real *8 elements long
c   w - must be at least 2*n+2 real *8 elements long
c   wsave - must be at least 4*n+20 real *8 elements long
c   z2 - must be at least 4*n+20 real *8 elements long
c 
        ier=0
        numit=25
c 
c        construct the smoothed curve (not properly closed)
c
        call rsres1(ier,x0,y0,n0,s,fis,n,
     1       w,wsave,nmin,z)
        if(ier.ne.0) return

c
c       conduct newton iterations to close the curve
c

        call evalerr(fis,n,dx,dy)
        err=sqrt(dx*dx+dy*dy)
cccc        call prin2('before newt, err = *',err,1)
c
        if(err.lt.eps) then
cccc        if(err.lt.eps.or.ier.eq.4) then
        call prinf('skipping newton*',1,0)
        goto 2200
        end if
c
        ifjump=0
        err2=1.0d20
        do 2000 j=1,numit

c       
        call newt(fis,n,err,nmin,z,wsave)

        niter=j
        errs(j)=err
cccc        call prin2('err = *',err,1)
c
        if(ifjump.eq.1) goto 2200
        if(err .lt. err2) goto 1900
        ier=8
        goto 2200
 1900   continue
        err2=err
        if(err2 .lt. eps) ifjump=1
 2000   continue
        ier=12
 2200   continue

c        
cccc        if (ier.ne.0) call prinf('after newt, ier=*',ier,1)
cccc        call prin2('after newt, err = *',err,1)
c 
c       from the obtained tangent angles, reconstruct the
c       curve
c 
        ifdtr=0
        if(nders .ge. 2) ifdtr=1
        call rsrec0(fis,z,n,wsave,z2,acc,ifdtr,w)
cccc        call prin2('in rsres2 after rsrec0, w=*',
cccc     1  w,n*2)
c

cccc        call prin2('after rsrec0, acc = *',acc,1)
 
        do 2400 i=1,n
        z(1,i)=z2(1,i)
        z(2,i)=z2(2,i)
 2400 continue
c 
c       construct the tangent vectors to the curve
c 
        do 2600 i=1,n
        der1(1,i)=dcos(fis(i))
        der1(2,i)=dsin(fis(i))
 2600 continue
c 
c       scale and shift the resampled curve
c
        call rsscale(z,n,x0,y0,n0,scale)
c 
c       calculate the sampling interval on the resampled
c       curve
c 
        done=1
        pi=datan(done)*4
        h=2*pi*scale/n
c 
c        if the user so requested - generate higher order
c        derivatives
c 
        if(nders. le. 1) return
        rlen=n*h
c 
ccc        call rsdiff(z,z2,n2,wsave,der,m,rlen,ider)
        ider=1
        call rsdiff(w,z2,n,wsave,der2,n,rlen,ider)
c
        return
        end
c 
c 
c 
c 
c 
        subroutine rsscale(z,n,x0,y0,n0,scale)
        implicit real *8 (a-h,o-z)
        save
        dimension z(2,1),x0(1),y0(1)
c 
c        this subroutine scales and shifts the resampled
c        curve so that its length is the same as that for
c        the user-specified polyginal one, and that their
c        center of mass coincide.
c 
c                     input parameters
c 
c  z - the curve to be rescaled
c  n - the number of elements in z
c  x0,y0 - user-supplied curve
c  n0 - the number of elements in each of arrays x0, y0
c 
c                     output parameters:
c 
c  z - the rescaled and shifted curve
c  scale - the coefficient by which z has been scaled
c
        done=1
        pi=4*atan(done)
c
c       . . . determine the length of the user-supplied polygon
c             and its center of mass
c 
        uslen=0
        cxu=0
        cyu=0
        do 1200 i=1,n0-1
        d=(x0(i+1)-x0(i))**2+(y0(i+1)-y0(i))**2
        d=sqrt(d)
        uslen=uslen+d
c 
        dx=(x0(i+1)+x0(i))/2
        dy=(y0(i+1)+y0(i))/2
c 
        cxu=cxu+dx*d
        cyu=cyu+dy*d
 1200 continue
c
        d=(x0(1)-x0(n0))**2+(y0(1)-y0(n0))**2
        d=sqrt(d)
        uslen=uslen+d
c
        dx=(x0(1)+x0(n0))/2
        dy=(y0(1)+y0(n0))/2
c
        cxu=cxu+dx*d
        cyu=cyu+dy*d
c
        cxu=cxu/uslen
        cyu=cyu/uslen
cccc        call prin2('uslen=*',uslen,1)
cccc        call prin2('cxu=*',cxu,1)
cccc        call prin2('cyu=*',cyu,1)
c 
c       polygonal length of the resampled curve is 2*pi
c
        reslen=2*pi
c 
c       rescale the resampled curve
c
        d=uslen/reslen
        scale=d
        do 1600 i=1,n+1
        z(1,i)=z(1,i)*d
        z(2,i)=z(2,i)*d
 1600 continue
c 
c       find the center of mass of the scaled resampled curve
c 
        cxr=0
        cyr=0
        do 1800 i=1,n
        cxr=cxr+z(1,i)
        cyr=cyr+z(2,i)
 1800 continue
        cxr=cxr/n
        cyr=cyr/n
c 
c       shift the scaled resampled curve
c 
        xshift=cxu-cxr
        yshift=cyu-cyr
        do 2000 i=1,n
        z(1,i)=z(1,i)+xshift
        z(2,i)=z(2,i)+yshift
 2000   continue
c
        return
        end
c 
c 
c 
c 
c 
        subroutine rsdiff(z,z2,n2,wsave,der,m,rlen,ider)
        implicit real *8 (a-h,o-z)
        save
        complex *16 z(1),z2(1),ima,der(1),wsave(1)
        data ima/(0.0d0,1.0d0)/
c 
c       this subroutine differentiates (in the frequency domain)
c       a user-specified curve and converts it into space domain.
c       it assumes that the curve is given to it in the
c       frequency domain. the number m of points at which the
c       derivative is returned does not have to be equal
c       to the number n of points where it has been evaluated,
c       but m must divide n.
c 
c                        input parameters:
c 
c   z - the curve to be differentiated (in the frequency domain)
c       attention!!! this parameter is changed by the subroutine!
c       see it in the list of output parameters.
c   n - the number of elements in array z
c   wsave - the array to be used by the fft routine dcfftb.
c       the user must have initialized it by calling the
c       routine dcffti (see).
c   m - the number of elements in the output array der.
c   rlen - the length of the curve being resampled.
c   ider - the order of derivaive being computed
c 
c                         output parameters:
c 
c   z - the curve differentiated in the frequency domain,
c       i.e. multiplied by i*k
c   der - the derivative of the input curve in the space domain.
c
        done=1
        pi=datan(done)*4
        nmult=n2/m
        do 1400 i=1,n2/2
        z2(i)=z(i)
        z2(n2-i+1)=z(n2-i+1)
        z2(i)=z(i)*ima*(i-1)
        z2(n2-i+1)=-z(n2-i+1)*ima*i
        z(i)=z2(i)
        z(n2-i+1)=z2(n2-i+1)
 1400 continue
        call dcfftb(n2,z2,wsave)
        d=done/n2/rlen*pi*2
        d=done/n2 *(pi*2/rlen)**ider
        do 1600 i=1,m
        j=(i-1)*nmult+1
        der(i)=z2(j)*d
 1600 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine rsres1(ier,x0,y0,n0,s,fis,n,
     1       w,wsave,nmin,z)
        implicit real *8 (a-h,o-z)
        real *8 x0(n0),y0(n0),s(n0+3),fis(1),wsave(1),z(1)
        complex *16 w(1),ima
c
        data ima/(0.0d0,1.0d0)/
c 
c       this subroutine performs a preliminary resampling,
c       taking a user-specified curve (given by a collection
c       of points in the plane) and producing a smooth (but
c       not closed) curve resampled at equispaced nodes.
c       the curve will have to be closed by subsequent
c       processing.
c 
c                  input parameters:
c 
c   x0,y0 - coordinates of the user-specified polygon
c   n0 - the number of elements in arrays x0,y0
c   n - the number of nodes in the equispaced resampling
c       of the curve
c   nmin - the parameter of the filtering subroutine rsfilt (see)
c 
c                  output parameters:
c 
c   ier - error return code.
c      ier=0 means successfull execution
c   s - the array of cumulative polygonal arc-lengths on the
c       user-supplied polygon (note that s(1)=0).
c   h - the sampling interval on the equispaced-resampled
c       curve
c   fis - the tangent angles of the resampled curve
c   wsave - the array initialized y the subroutine dcffti,
c       that can be fed into dcfftf, dcfftb
c 
c                  work arrays:
c 
c   w - must be at least n*2+4 real *8 locations long
c   z - must be at least n*2+4 real *8 locations long. this array
c       is only used if ifz2=1 has been specified (see above)
c       otherwise, it is not used.
c
        done=1
        pi=atan(done)*4
c 
c       construct the polygonal resampling of the curve
c
        call repoly(x0,y0,n0,s,n,h,fis,w)
cccc        call prin2('in rsres1 after repoly, s = *',s,n0+1)
cccc        call prin2('in rsres1 after repoly, fis = *',fis,n+1)
c
c       note that while the angles are periodic, this assumes
c       that fis=0 is the same as fis=2*pi;
c       fft is not aware of this and hence, we need to
c       "periodize" the angles
c
c       determine whether input data is given
c       clockwise or counterclockwise
c
        thd=2*pi
        d1=fis(n)-(fis(1)+2*pi)
        d2=fis(n)-(fis(1)-2*pi)
        if (abs(d2).lt.abs(d1)) thd=-2*pi
        htet=thd/n
c
        do 1200 i=1,n
        w(i)=fis(i)-(i-1)*htet
 1200   continue
cccc        call prin2('before dcfftf, w=*',w,n*2)
c 
        call dcfftf(n,w,wsave)
cccc        call prin2('before filtering, w=*',w,n)
c
cccc        call rsfilt(w,n,nmin)
        call rsfilg(w,n,nmin)
cccc        call rsfilk(w,n,nmin)
cccc        call prin2('after filtering, w = *',w,n)
c
        d=1
        d=d/n
        call dcfftb(n,w,wsave)
        do 1300 i=1,n
        w(i)=w(i)*d
 1300   continue
cccc         call prin2('filtered fis are*',w,n+1)
c 
c        construct the tangent vectors for the filtered
c        curve and reconstruct the curve from them
c 
        do 1400 i=1,n
        fis(i)=w(i)+(i-1)*htet
 1400   continue
        fis(n+1)=fis(1)+thd
c
        return
c
        end
c 
c 
c 
c 
c 
        subroutine rsrec0(thetas,z,n,wsave,z2,acc,ifdtr,
     1    dertrans)
        implicit real *8 (a-h,o-z)
        save
        dimension wsave(*),thetas(*),z(2,*),z2(2,*),
     1    dertrans(*)
c 
c       this subroutine reconstructs a closed curve in r^2
c       from its tangent angles
c 
c                     input parameters:
c 
c  thetas - the tangent angles
c  n - the number of nodes in the discretization of the boundary
c       (the number of elements in the array thetas)
c  wsave - the information array to be used y the subroutines
c       dcfftf, dcfftb. this array must have been precomputed
c       by a prior call to the suroutine dcffti
c 
c                     output parameters:
c 
c  z2 - the reconstructed curve
c  acc - the accuracy of the reconstruction, in the sense
c       that acc is the estimated value of the fourier series
c       of z (as a function of arclength of the curve) near n/2
c 
c                     work array:
c 
c  z - must be 2*n+2 real *8 locations long
c 
c 
c       . . . construct the tangents to the curve from the
c             tangent angles
c
        do 1400 i=1,n
        z(1,i)=dcos(thetas(i))
        z(2,i)=dsin(thetas(i))
 1400 continue
c 
c        integrate the tangents
c
        call rsintc(n,z2,z,wsave,acc,ifdtr,dertrans)
ccc         call prin2('after rsintc, acc=*',acc,1)
        return
        end
c 
c 
c 
c 
c 
       subroutine rsintc(n,xout,w,wsave,acc,ifdtr,dertrans)
       implicit real *8 (a-h,o-z)
        save
       complex *16 w(1),wsave(1),c,ima,xout(1),d,h,dertrans(1)
        data ima/(0.0d0,1.0d0)/
c 
c        this subroutine calculates the values of an indefinite
c        integral of a complex periodic function, given by its
c        fourier transform, at a bunch of equispaced nodes.
c 
c 
c                      input parameters:
c 
c  n - the number of elements in arrays w, xout (see below).
c        it is important that n be a power of two, or at least
c        be a product of small primes. otherwize, the subroutine
c        can be     v e r y    slow.
c  w - array containing the (complex)  function
c        whose integral is being evaluated. the
c        contents of w are destroyed by this subroutine.
c  wsave - the array containing parameters created by the
c        entry dcffti and used by the entries dcfftf, dcfftb.
c 
c                       output parameters:
c 
c  xout - the table of values of the indefinite integral
c  acc - the accuracy of the integration, in the sense
c       that acc is the estimated value of the fourier series
c       of z (as a function of arclength of the curve) near n/2
c 
c 
c       . . .  fourier-transform the input function
c 
        ifin=0
        iftran=1
        if(iftran .eq. 0) goto 1300
        if(ifin .ne.0) call dcffti(n,wsave)
        call dcfftf(n,w,wsave)
c 
c        if the user so requested - return to the user
c        the fourier transform of the input vector
c 
        if(ifdtr .eq. 0) goto 1300
        do 1200 i=1,n
        dertrans(i)=w(i)
 1200 continue
c 
 1300 continue
cccc         call prin2('in rsintc, ffted w is*',w,n*2)
c 
c       integrate the function in the dual domain
c 
        c=w(1)
ccc          call prin2('in rsintc, c=*',c,2)
        w(1)=0
        w(n)=-w(n)/ima
        do 1400 i=2,n/2
        w(i)=w(i)/(ima*(i-1))
        w(n-i+1)=-w(n-i+1)/(ima*i)
 1400 continue
c 
c        determine the accuracy of the discretization
c 
        acc=0
        do 1500 i=-3,3
        dd=cdabs(w(n/2+i))
        if(dd .gt. acc) acc=dd
 1500 continue
cccc        call prin2('in rsintc, acc=*',acc,1)
c 
c       fourier-transform the thing back
c 
        call dcfftb(n,w,wsave)
c 
c       scale and add in the linear term
c
        done=1
        pi=datan(done)*4
        coef=1
        coef=coef/ n
        rn=n
        h=2*pi*c/rn**2
c        h=2*pi*c/n**2
        do 1600 i=1,n
        xout(i)=w(i)*coef+(i-1)*h
 1600 continue
        xout(n+1)=2*pi*c/n + xout(1)
c 
c        now, the particular representation of the indefinite
c        integral is non-zero at zero. shift it so it will
c        be that way
c 
        d=xout(1)
        do 1800 i=1,n+1
        xout(i)=xout(i)-d
 1800 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine rsfilt(w,n,nmin)
        implicit real *8 (a-h,o-z)
        save
        complex *16 zero,w(1)
c 
c       this subroutine filters a function in the frequency
c       domain by applying the standard cosine filter. the filer
c       starts at nmin-th frequency and is nmin frequencies wide.
c 
c                        input parameters:
c 
c  w - the fourier transform of the function being filtered
c  n - the number of elements in w
c  nmin - the lowest frequency at which the filter is .neq. 1
c 
c                        output parameters:
c 
c  w - the filtered fourier transform of the functoin
c 
c 
c       . . . set up the parameters of the bell
c 
        done=1
        half=done/2
        pi=datan(done)*4
        h=pi/(nmin-1)
c 
c       impose the bell
c 
        do 1200 i=1,nmin
        d=i*h
        dd=(dcos(d)+done)*half
        w(nmin+i)=w(nmin+i)*dd
        w(n-nmin-i)=w(n-nmin-i)*dd
 1200 continue
c 
c        set to zero everything outside the bell
c 
        zero=0
        do 1400 i=2*nmin, n-2*nmin
        w(i)=zero
 1400 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine rsfilg(w,n,nmin)
        implicit real *8 (a-h,o-z)
cccc        complex *16 zero,w(1)
        complex *16 w(1)
  
c       this subroutine filters a function in the frequency
c       domain by applying the standard gaussian filter. the filter
c       starts at nmin-th frequency and is nmin frequencies wide.
c 
c                        input parameters:
c 
c  w - the fourier transform of the function being filtered
c  n - the number of elements in w
c  nmin - the lowest frequency at which the filter is .neq. 1
c 
c                        output parameters:
c 
c  w - the filtered fourier transform of the function
c 
c 
c       . . . set up the parameters of the bell
c 
        done=1
        d=10
        rnmin=nmin
        a=dlog(d)/rnmin**2
c 
c       impose the bell
c 
        dd0=1
        do 1200 i=1,n/2
        d=i
        dd=exp(-a*d**2)
        w(i+1)=w(i+1)*dd
        w(n-i+1)=w(n-i+1)*dd
        dd0=dd0+dd
 1200 continue
        if(2.ne.3) return
c 
        dd0=1/dd0
        do 1400 i=1,n
        w(i)=w(i)*dd0
 1400 continue
        return
        end
c
c
c
c
c
        subroutine rsfilk(w,n,nmin)
        implicit real *8 (a-h,o-z)
        complex *16 w(1)
c
c       this subroutine filters a function in the frequency
c       domain by applying the kaiser window. the filter is
c       defined as
c
c              I_{0} (\beta*\sqrt{1-(\frac{x}{m/2})^2}) \over    (1)
c                    I_{0} (\beta)
c
c       where I_{0} is the zeroth order modified Bessel function.
c
        done=1
        dc=0.1d0
c
c       determine beta; specifically, find beta so that (1)
c       evaluated at x=nmin equals to dc
c
cccc        beta=2.0d4
        bini=done
        rm=n/2
        rnm=nmin
        d0=rnm/rm
        d0=sqrt(1-d0**2)
c
        do 400 i=1,20
c
        call besseli(bini*d0,d1,d1p)
        call besseli(bini,d2,d2p)
        ds=exp(bini*(d0-1))
c
        ff=d1/d2*ds-dc
        if (abs(ff).le.1d-10) goto 500
c
        fp=(d0*d2*d1p-d1*d2p)/d2**2
        fp=fp*ds
c
        bini=bini-ff/fp
c
 400    continue
c
        call prinf('newton could not find beta*',1,0)
        stop
c
 500    continue
c
c       do three more iterations for good measure
c
        do 600 i=1,3
c
        call besseli(bini*d0,d1,d1p)
        call besseli(bini,d2,d2p)
        ds=exp(bini*(d0-1))
c
        ff=d1/d2*ds-dc
c
        fp=(d0*d2*d1p-d1*d2p)/d2**2
        fp=fp*ds
c
        bini=bini-ff/fp
c
 600    continue
c
        beta=bini
c
        do 1200 i=1,n/2
c
        d=i
        d0=d/rm
        d0=sqrt(1-d0**2)
c
c       compute numerator and denomenator
c
        if (i.eq.n/2) goto 1300

        call besseli(beta*d0,d1,d1p)
        call besseli(beta,d2,d2p)
        goto 1500
c
 1300   continue
        d1=done
c
 1500   continue
c
c       rescale e^-x term
c
        ds=exp(beta*(d0-1))
        dd=d1/d2*ds
c
        w(i+1)=w(i+1)*dd
        w(n-i+1)=w(n-i+1)*dd
c
 1200   continue
c
        return
c
        end
c
c
c
c
c
        subroutine besseli(x,y,yp)
        implicit real *8 (a-h,o-z)
        real *8 funs(10000)
c 
c       this subroutine evaluates the scaled zeroth order I-type
c       (modified) Bessel functions of the argument x. More
c       specifically, on exit, y and yp are
c 
c              e^{-x} \cdot I_{j} (x)                             (1)
c
c       for j=0,1 respectively.
c 
c                       input parameters:
c 
c  x - the argument of which the zeroth order modified Bessel functions
c       is to be evaluated.
c 
c                        output parameters:
c 
c  y - zeroth order modified Bessel function evaluated at x
c  yp - first order modified Bessel function evaluated at x


c
c       this code is not currently robust
c


c 
c       . . . recurse up till the thing becomes fairly large
c
        done=1
c 
        funs(1)=0
        funs(2)=1
        rlarge=1d40
c 
        do 1200 i=2,10000
        ni=i
        funs(i+1)=-2*(i-1)/x*funs(i)+funs(i-1)
        if(funs(i+1) .gt. rlarge) goto 1400
 1200   continue
        call prinf('in besseli0, max iter reached*',1,0)
        stop
 1400   continue
cccc        call prin2('funs after firt recursion*',funs,ni)
c 
c       now, starting with the position ni, recurse down
c 
        funs(ni+1)=0
        funs(ni)=1
c 
        do 1600 i=1,ni-1
        funs(ni-i)=2*(ni-i)/x*funs(ni-i+1)+funs(ni-i+2)
 1600   continue
c 
c       sum up the the recursed values, to get
c       the normalizing coefficient
c 
        sum=funs(1)
        do 1800 i=2,ni
        sum=sum+2*funs(i)
 1800   continue
        sum=done/sum
c 
c       . . . scale
c 
        do 2000 i=1,ni
        funs(i)=funs(i)*sum
 2000   continue
c
        y=funs(1)
        yp=funs(2)
c
        return
c
        end










c 
c 
c 
c 
c 
        subroutine repoly(x,y,n,s,m,h,fis,thetas)
        implicit real *8 (a-h,o-z)
        dimension x(1),y(1),s(1),fis(1),thetas(1)
c 
c       this subroutine for a user-specified polygonal curve
c       (x,y) constructs the array s of its cumulative arc-lengths,
c       and the array fis of its tangent angles;
c       the user-provided polygon will be treated as closed, i.e.
c       x(n+1)=x(1),
c       y(n+1)=y(1).
c       note also that m must be greater than n
c 
c                  input parameters:
c 
c   x,y - coordinats of the user-specified polygon.
c   n - the number of elements in arrays x,y
c   m - the number of nodes in the equispaced resampling
c       of the curve
c 
c                  output parameters:
c 
c   s - the array of cumulative polygonal arc-lengths on the
c       user-supplied polygon (note that s(1)=0).
c   h - the sampling interval on the equispaced-resampled
c       curve
c   fis - the tangent angles of the resampled curve
c 
c                  work arrays:
c 
c   thetas - must be at least m+2 real *8 locations long
c 
cccc        call prin2('in repoly, x=*',x,n)
cccc        call prin2('in repoly, y=*',y,n)
cccc        call prinf('in repoly, n=*',n,1)
cccc        call prinf('in repoly, m=*',m,1)
c 
c       construct the array of user-provided lengths for all
c       nodes of the user-provided polygon
c 
        s(1)=0
        do 1200 i=2,n
        d=(x(i)-x(i-1))**2+(y(i)-y(i-1))**2
        s(i)=s(i-1)+sqrt(d)
 1200   continue
c
        d=(x(1)-x(n))**2+(y(1)-y(n))**2
        s(n+1)=s(n)+sqrt(d)
cccc         call prin2('inside repoly, s=*',s,n+1)
c 
c       determine the sampling interval on the resampled curve
c 
        h=s(n+1)/m
c 
c        construct the angles between the x axis and the
c        various sides of the polygon
c 
c       . . . construct the angles between the consecutive
c             sides of the polygon
c 
        call reancr(x,y,n,fis)
ccc         call prin2('angles between sides of polygon are*',fis,n+1)
c 
c       find the angle between the first side of the polygon and
c       the x axis
c 
        wx0=x(1)-1
        wy0=y(1)
c 
        call reangl(wx0,wy0,x(1),y(1),x(2),y(2),thetas(1))
c 
cccc        subroutine reangl(x1,y1,x2,y2,x3,y3,fi)
c 
c        construct the angles
c 
         do 2200 i=2,n+1
         thetas(i)=thetas(i-1)+fis(i)
 2200    continue
cccc         call prin2('tangent angles of original polygon are*',
cccc     1    thetas,n+1)
c 
c       construct the angles between the x axis and the various
c       sides of the polygon, but tabulated at the nodes of
c       the finer discretization
c 
        j=1
ccc         call prin2('in repoly, h=*',h,1)
        do 2400 i=1,m
        d=(i-1)*h
        if(d.gt.s(j+1)) j=j+1
        fis(i)=thetas(j)
 2400   continue
        fis(m+1)=thetas(n+1)
c 
        return
c
        end
c 
c 
c 
c 
c 
        subroutine reancr(x,y,n,fis)
        implicit real *8 (a-h,o-z)
        save
        dimension x(1),y(1),fis(1)
c 
c       for the user-specified polygon, this subroutine
c       constructs the angles made by the consecutive sides
c       of the polygon. the polygon is assumed to be closed.
c 
c                  input parameters:
c 
c   x,y - the coordinates of the vertices of the user-supplied
c         polygon
c   n - the number of elements in arrays x,y
c 
c                  output parameters:
c 
c   fis - the angles made by consecutive sides of the polygon
c 
c      remark:
c       the polygon is treated as closed, i.e.
c       x(n+1)=x(1),
c       y(n+1)=y(1).
c 
c       . . . initialize the angle-computing routine and
c             calculate the angle between the last and the first
c             sides. this angle is declared to be at the
c             first vertex.
c 
        call reangi(ier)
        call reangl(x(n),y(n),x(1),y(1),x(2),y(2),fis(1))
c 
c       calculate the rest of the angles
c 
        do 1200 i=2,n-1
        call reangl(x(i-1),y(i-1),x(i),y(i),x(i+1),y(i+1),fis(i))
cccc        call prin2('after reangl, fis(i)  is*',fis(i),1)
 1200   continue
c
        call reangl(x(n-1),y(n-1),x(n),y(n),x(1),y(1),fis(n))
        fis(n+1)=fis(1)
c
        return
c
        end
c 
c 
c 
c 
c 
        subroutine reangl(x1,y1,x2,y2,x3,y3,fi)
        implicit real *8 (a-h,o-z)
        save
        data eps/1.0d-30/
c 
c 
c       this subroutine finds the angle fi between the
c       line segments  [(x1,y1),(x2,y2)] and  [(x2,y2),(x3,y3)]
c       initialization entry point reangi (see below) has to be
c       called prior to the invocation of this main entry
c 
c        . . . find the (unscaled) sin and cos of the angle
c
        ux=x2-x1
        uy=y2-y1
c 
        vx=x3-x2
        vy=y3-y2
c 
        sinfi=ux*vy-vx*uy
        cosfi=ux*vx+uy*vy
c 
c       if the angle is not right - find it and exit
c 
          if(dabs(cosfi) .lt. eps) goto 2000
        fi=datan(sinfi/cosfi)
cccc    if(cosfi .le. 0) fi=fi+pi
        if( (cosfi .le. 0) .and. (sinfi.ge.0) ) fi=fi+pi
        if( (cosfi .le. 0) .and. (sinfi.le.0) ) fi=fi-pi
        return
 2000 continue
c 
c       the angle is right. find it
c 
cccc    call prin2('in reangl, the angle is right. cosfi is*',cosfi,1)
        fi=pi2
        if(sinfi .lt. 0.0d0) fi=-fi
ccc     if(sinfi .gt. 0.0d0) fi=-fi
        return
c 
c 
c 
c 
        entry reangi(ier)
c 
c       this is initialization entry point. it has to be called
c       prior to any calls to the main entry reangl (above).
c 
        done=1
        pi=datan(done)*4
        pi2=pi/2
        return
        end
c
c
c
c
c
        subroutine blapp(fis,n,w,wsave,nbl)
c
        implicit real *8 (a-h,o-z)
        real *8 fis(n),wsave(1)
        complex *16 w(n)
c
c       compute supergain given angles
c
        done=1
        pi=atan(done)*4
c
c       determine whether input data is given
c       clockwise or counterclockwise
c
        thd=2*pi
        d1=fis(n)-(fis(1)+2*pi)
        d2=fis(n)-(fis(1)-2*pi)
        if (abs(d2).lt.abs(d1)) thd=-2*pi
        htet=thd/n
c
        do 100 i=1,n
        w(i)=fis(i)-(i-1)*htet
 100    continue
c
        call dcfftf(n,w,wsave)
c
c       need to adjust for fft not being unitary
c
        rn=n
        rn=sqrt(rn)
        do 200 i=1,n
        w(i)=w(i)/rn
 200    continue
c
c       compute approximate band-limit
c       skipping first element seems to output
c       bandlimit approximately equal to nmin
c
        s2=0
        do 300 i=2,n
c        wr=w(i)
c        wi=-ima*w(i)
c        s2=s2+(wr*wr+wi*wi)
        s2=s2+abs(w(i))
 300   continue
c
cccc        call prin2_long('s2 = *',s2,1)
c
        rperc=.99
        d=rperc*s2
c
cccc        call prin2_long('desired s1 = *',d,1)
c
c       use bisection to find bandlimit such that
c       the energy is rperc of the total energy
c
cccc        m=72
        nu=n/2+1
        nl=1
c
        m=nu/2
c
        rnu=nu
        riter=log(rnu)/log(2.0d0)
        niter=riter+2
c
        do 600 j=1,niter
c
        s1=0
        do 400 i=2,m
        s1=s1+abs(w(i))
 400    continue
c
        do 500 i=2,m
        s1=s1+abs(w(i))
 500    continue
c
cccc        call prin2_long('s1-d = *',s1-d,1)
c
        if (s1-d.gt.0) then
        nu=m
        m=(nl+m)/2
        end if
c
        if (s1-d.le.0) then
        nl=m
        m=(m+nu)/2
        end if
c
 600    continue
c
        nbl=m
c
        return
c
        end
c
c
c
c
c
        subroutine newt(thetas,n,err,nmin,z,wsave)
c
        implicit real *8 (a-h,o-z)
        real *8 thetas(*),z(2,*),wsave(*)
c
c       compute the perturbing tangent angles
c

        do 200 i=1,n
          z(1,i)=cos(thetas(i))
          z(2,i)=sin(thetas(i))
 200    continue

 
c 
c       filter the perturbing tangent angles
c 
        call dcfftf(n,z,wsave)
c
cccc        call rsfilt(z,n,nmin)
        call rsfilg(z,n,nmin)

c
        d=1
        d=d/n
        call dcfftb(n,z,wsave)
c
        do 400 i=1,n
        z(1,i)=z(1,i)*d
        z(2,i)=z(2,i)*d
 400    continue
c
c       perform newton iteration
c
        call evalerr(thetas,n,dx,dy)
c
c       construct 2x2 matrix
c
        a11=0
        a12=0
        a21=0
        a22=0
        do 1000 i=1,n
        a11=a11+sin(thetas(i))*z(1,i)
        a12=a12+sin(thetas(i))*z(2,i)
        a21=a21-cos(thetas(i))*z(1,i)
        a22=a22-cos(thetas(i))*z(2,i)
 1000   continue

c
cccc        call prin2('a11 = *',a11,1)
cccc        call prin2('a12 = *',a12,1)
cccc        call prin2('a21 = *',a21,1)
cccc        call prin2('a22 = *',a22,1)
c
c       find inverse of 2x2 matrix
c
        deta=a11*a22-a12*a21
        b11=a22/deta
        b12=-a12/deta
        b21=-a21/deta
        b22=a11/deta
c
        delx=b11*dx+b12*dy
        dely=b21*dx+b22*dy
c
cccc        call prin2('delx = *',delx,1)
cccc        call prin2('dely = *',dely,1)
c
c       update thetas
c
        do 1400 i=1,n
        thetas(i)=thetas(i)+delx*z(1,i)+dely*z(2,i)
 1400   continue
c

        call evalerr(thetas,n,dx,dy)
        err=sqrt(dx*dx+dy*dy)
c
        return
c
        end




        subroutine evalerr(thetas,n,dx,dy)
c
        implicit real *8 (a-h,o-z)
        real *8 thetas(*)
c
        dx=0
        dy=0
        do 1000 i=1,n
          dx=dx+cos(thetas(i))
          dy=dy+sin(thetas(i))
 1000 continue
c
        return
c
        end
c
c
c
c
c
        subroutine rsfis(n,dxy,fis)
c
        implicit real *8 (a-h,o-z)
        real *8 dxy(2,n), fis(n)
c
        done=1
        pi=atan(done)*4
        eps=1.0d-30
c
c       compute the first angle
c
c
c       if the angle is not right - find it
c
        i=1
        if(abs(dxy(1,i)).lt.eps) goto 550
        th=dxy(2,i)/dxy(1,i)
        fis(i)=datan(th)
        if( (dxy(1,i).le.0) .and. (dxy(2,i).ge.0) ) fis(i)=fis(i)+pi
        if( (dxy(1,i).le.0) .and. (dxy(2,i).le.0) ) fis(i)=fis(i)-pi
        goto 350
 550    continue
c
c       the angle is right - find it
c
        fis(i)=pi/2
        if(dxy(2,i).lt.0.0d0) fis(i)=-fis(i)
 350    continue
c
c       compute the rest of the angles
c
        do 200 i=2,n
c
c       if the angle is not right - find it
c
        if(abs(dxy(1,i)).lt.eps) goto 500
        th=dxy(2,i)/dxy(1,i)
        fis(i)=datan(th)
        if( (dxy(1,i).le.0) .and. (dxy(2,i).ge.0) ) fis(i)=fis(i)+pi
        if( (dxy(1,i).le.0) .and. (dxy(2,i).le.0) ) fis(i)=fis(i)-pi
        goto 300
 500    continue
c
c       the angle is right - find it
c
        fis(i)=pi/2
        if(dxy(2,i).lt.0.0d0) fis(i)=-fis(i)

 300    continue
c
c       since pi and -pi are the same (but are not the same
c       to fft), we enforce continuity by shifting angle by
c       +-2*pi and choosing the one with the smallest
c       difference to previous point
c
        dd=abs(fis(i)-fis(i-1))
        dd0=abs(fis(i)+2*pi-fis(i-1))
        if(dd0.lt.dd) then
        fis(i)=fis(i)+2*pi
        dd=dd0
        end if

        dd0=abs(fis(i)-2*pi-fis(i-1))
        if(dd0.lt.dd) then
        fis(i)=fis(i)-2*pi
        end if
c
 200    continue
c
        return
c
        end
c
c
c
c
c
        subroutine anarsbl(ier,funcurve,par1,par2,rl,n,eps,t,h,rltot,
     1       m,xgs,whts,wright,wvals,lenwr,lenwv,nints)
        implicit real *8 (a-h,o-z)
        save
        dimension t(1),par1(1),par2(1),xgs(1),whts(1),
     1       wright(1),wvals(1)
        external funcurve
c 
c       This subroutine produces an equispaced (in terms
c       of the arc length) resampling of a user-specified
c       curve. The curve is specified in the form of the
c       subroutine funcurve providing some (generally,
c       not equispaced) parametrization of the curve, and
c       the length rl of the interval of parametrization.
c       Explanation: this subroutine assumes that the curve
c       to be resampled is defined by the mapping
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
c       ier=16000 means that the length of the user-provided wright
c               or wvals array is insufficient. this is a fatal error
c  t - the equispaced discretization of the curve. note that the
c       subroutine will return n+1 values in array t. that is, if the
c       curve is a closed one, then t(n+1) will be equal to t(1).
c  h - the sampling interval along the arc of the resampling returned
c       in array t
c  m - number of gaussian nodes used for integration
c  xgs - gaussian nodes
c  whts - gaussian weights
c 
c       . . . determine the total length of the curve
c 
        ier=0
        jer=0
        ker=0
        a=0
        b=rl
cccc        m=12
        m=6
cccc        epsgauss=eps*1000
        epsgauss=eps
        epsnewt=dsqrt(eps)/1000
c
c       compute gaussian nodes and weights
c
        itype=1
        call legeexps(itype,m,xgs,unil,vnil,whts)
c
        call anarsblga(ier,a,b,funcurve,par1,par2,m,xgs,whts,epsgauss,
     1      rltot,wright,wvals,lenwr,lenwv,maxrec,numint,nints)
        if (ier.ne.0) return
cccc        call prin2('after first anarsblga, rltot=*',rltot,1)
cccc        call prinf('after first anarsblga, maxrec=*',maxrec,1)
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
c
        do 3000 k=1,n
cccc         call prinf('in anarsbl, k=*',k,1)
c 
c       determine the dl/dt at the point t(k)
c 
        call funcurve(t(k),par1,par2,x,y,dxdt,dydt,curv)
c 
c       get the initial approximation for the next point
c 
        dldt=dsqrt(dxdt**2+dydt**2)
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
        call findint(wright,nints,t(k)+dt,ij0,ijk)
c
        if (ij0.eq.ijk) then
        call anarsblg2(t(k),t(k)+dt,funcurve,
     1       par1,par2,xgs,whts,m,htrue2)
        goto 2300
        end if
c
        call anarsblg2(t(k),wright(ij0),funcurve,
     1       par1,par2,xgs,whts,m,htrue2)
c
        do 2200 kk=ij0+1,ijk-1
        htrue2=htrue2+wvals(kk)
 2200   continue
c
        call anarsblg2(wright(ijk-1),t(k)+dt,funcurve,
     1       par1,par2,xgs,whts,m,htrue3)
        htrue2=htrue2+htrue3
c
 2300   continue
        htrue=htrue2
c 
c       determine the dl/dt at that point and make a
c       newton step
c 
        call funcurve(t(k)+dt,par1,par2,x,y,dxdt,dydt,curv)
        dldt=dsqrt(dxdt**2+dydt**2)
c 
        derr=htrue-h
cccc        call prin2('htrue-h=*',derr,1)
c 
        if(dabs(derr) .lt. dabs(derrold) ) goto 1300
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
        if(dabs(derr) .lt. epsnewt) goto 2000
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
cccc        call prinf('exiting anarsbl, ier=*',ier,1)
c
        return
c
        end
c
c
c
c
c
        subroutine findint(w,n,t,ii,ij)
c
        implicit real *8 (a-h,o-z)
        real *8 w(1)
c
        do 100 i=ii,n
        ij=i
        if (t .le. w(i)) return
 100    continue
c
        return
c
        end
c 
c 
c 
c 
c 
         subroutine anaptbl(ier,x,n,h,ts,funcurve,par1,par2,eps,
     1       mm,xgs,whts,wright,wvals,nints,rl,tout,
     2       xout,yout,dxdtout,dydtout,curvout)
         implicit real *8 (a-h,o-z)
         real *8 ts(1),xgs(1),whts(1),wright(1),wvals(1)
         real *8 tz(10),xz(10),cz(10)
         external funcurve
c 
c       This subroutine finds the location (both on the
c       parametrizing interval [0,rl] and in R^2) of the
c       point specified by the user on the curve. The
c       curve has to have been preprocessed by a preceding
c       call to the subroutine anarsbl (see), and the point
c       is specified by its distance (in terms of arc-length)
c       from the beginning of the curve. PLEASE NOTE THAT THIS
C       SUBROUTINE DEPENDS ON THE SUBROUTINE ANARSBL FOR THE
C       PARAMETERS H,TS; THIS SUBROUTINE HAS NO FUNCTION AS
C       A STAND-ALONE DEVICE.
c 
c                        Input parameters:
c
c  x - the location (in terms of arc-length) on the curve of the
c       point to be found
c  n - the number of nodes in the equispaced discretization of the
c       curve produced by a preceding call to anarsbl
c  h - the sampling interval along the arc of the resampling returned
c       in array ts by the subroutine anarsbl
c  ts - the equispaced discretization of the curve produced by the
c       subroutine anarsbl
c  funcurve - the subroutine providing a parametrization
c        of the curve. the calling sequence of funcurve
c        must be
c 
c                  funcurve(t,par1,par2,x,y,dxdt,dydt,curv),
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
c       curv - the curvature
c  eps - the precision to which the point is to be found. Generally, it
c       is a good idea to make this parameter the same as in the
c       preceding call to anarsbl
c 
c                         Output parameters:
c 
c  tout - the value of the curve parameter corresponding to the location
c       x along the arc-length
c  xout, yout - the coordinates in R^2 of the point x on the curve
c  dxdt, dydt - the coordinates of the tangent (of length 1) to the
c       curve at the point (xout,yout)
c

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
        if(m .gt. n) ier=512
        if(m .gt. n) return

        if(m.eq.n) tk = rl
        if(m.eq.n) goto 3000

        t1 = ts(m+1)
        x1 = m*h
c
c       use cubic interpolation to find the approximate location
c       of the curve parameter (as opposed to the arc-length)
c       corresponding to the user-supplied arc-length
c
        if (m.ge.1) tz(1)=ts(m)
        if (m.eq.0) tz(1)=ts(n)-ts(n+1)
        tz(2)=ts(m+1)
        tz(3)=ts(m+2)
        if(m.le.n-2) tz(4) = ts(m+3)
        if(m.eq.n-1) tz(4) = ts(n+1)+ts(2)
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
cccc        epsgauss=eps*1000
        epsgauss=eps
        epsnewt=dsqrt(eps)/1000
        maxit=20
        fk=h
c
        ij0=1
        call findint(wright,nints,t1,ij0,ija)
c 
        do 1400 i=1,maxit
c 
c       evaluate the function whose root we are seeking
c 
        step=1
        do 1200 j=1,20
c
        call findint(wright,nints,tk,ija,ijb)
c
c
        if (ija.eq.ijb) then
        call anarsblg2(t1,tk,funcurve,
     1       par1,par2,xgs,whts,mm,rltot)
        goto 2300
        end if
c
        call anarsblg2(t1,wright(ija),funcurve,
     1       par1,par2,xgs,whts,mm,rltot)
c
        do 2200 kk=ija+1,ijb-1
        rltot=rltot+wvals(kk)
 2200   continue
c
        call anarsblg2(wright(ijb-1),tk,funcurve,
     1       par1,par2,xgs,whts,mm,rltot2)
        rltot=rltot+rltot2
c
 2300   continue
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
        call funcurve(tk,par1,par2,xk,yk,dxdt,dydt,curv)
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

 3000 continue
c
        tout=tk
        curvout = 0
        call funcurve(tout,par1,par2,xout,yout,dxdtout,dydtout,curvout)
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
        subroutine anarsblga(ier,a,b,fun,par1,par2,m,t,w,eps,
     1      rint,wright,wvals,lenwr,lenwv,maxrec,numint,nints)
c
        implicit real *8 (a-h,o-z)
        dimension t(100),w(100),stack(400),vals(200),idep(200),
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
c        fun(t,x,y,dxdt,dydt,curv)
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
        call anarsblg1(ier,stack,a,b,fun,par1,par2,t,w,m,
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
        subroutine anarsblg1(ier,stack,a,b,fun,par1,par2,t,w,m,
     1      vals,nnmax,eps,rint,wright,wvals,lenwr,lenwv,nints,idepth,
     2      maxdepth,maxrec,numint)
        implicit real *8 (a-h,o-z)
        save
        dimension stack(2,1),t(1),w(1),vals(1),par1(1),par2(1),
     1       wright(1),wvals(1),idepth(1)
c 
c       start the recursion
c 
        stack(1,1)=a
        stack(2,1)=b
        idepth(1)=0
        call anarsblg2(a,b,fun,par1,par2,t,w,m,vals(1))
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
        if(idepth(j) .gt. maxrec) maxrec=idepth(j)
cccc        call prinf('j=*',j,1)
c 
c       subdivide the current subinterval
c 
        c=(stack(1,j)+stack(2,j))/2
        call anarsblg2(stack(1,j),c,fun,
     1       par1,par2,t,w,m,value2)
c 
        call anarsblg2(c,stack(2,j),fun,
     1       par1,par2,t,w,m,value3)
c 
        dd=dabs(value2+value3-vals(j))
cccc         call prin2('in anarsblg1, dd=*',dd,1)
c
        ifdone=0
        if(dd .le. eps) ifdone=1
c 
c       if the function on this subinterval has been
c       integrated with sufficient accuracy - add the
c       value to that of the global integral and move up
c       in the stack
c 
        if(ifdone .eq. 0) goto 2000
c
        nints=nints+1
        if(nints.le.lenwr.and.nints.le.lenwv) goto 1200
        ier=16000
        return
 1200   continue
        wright(nints)=stack(2,j)
        wvals(nints)=value2+value3
c 
        rint=rint+value2+value3
        j=j-1
c 
c       if the whole thing has been integrated - return
c 
        if(j .eq. 0) return
        goto 3000
c
 2000 continue
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
        if(j .le. maxdepth) goto 3000
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
        subroutine anarsblg2(a,b,fun,par1,par2,t,w,m,rint)
        implicit real *8 (a-h,o-z)
        save
        dimension t(1),w(1),par1(1),par2(1)
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
        call fun(tt,par1,par2,x,y,dxdt,dydt,curv)
c
        rint=rint+w(i)*dsqrt(dxdt**2+dydt**2)
c
 1200   continue
c
        rint=rint*u
c
        return
c
        end
