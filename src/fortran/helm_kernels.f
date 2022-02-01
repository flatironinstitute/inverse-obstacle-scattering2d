      subroutine slp(srcinfo,targinfo,dpars,zk,ipars,u)
c
c       single layer interaction kernel
c
c       input:
c         srcinfo(2) - double
c           x,y location of source
c         targinfo(2) - double
c           x,y location of target
c         dpars - double
c           dummy parameter
c         zk - complex
c           Helmholtz parameter
c         ipars - integer
c           dummy parameter
c
c       output:
c         u = i/4 H_{0}(zk*|r|) \, ,
c         where r is the distance between source and target
c          
c
      implicit real *8 (a-h,o-z)
      integer ipars
      real *8 dpars,srcinfo(2),targinfo(2)
      complex *16 u,h0,ima,zs,z,zk,h1
      data ima/(0.0d0,1.0d0)/
      data zs/(0.0d0,0.25d0)/

      rr2 = (srcinfo(1)-targinfo(1))**2 + (srcinfo(2)-targinfo(2))**2
      rr = dsqrt(rr2)
      z = zk*rr
      ifexpon = 1

      call hank103(z,h0,h1,ifexpon)
      u = zs*h0
      

      return
      end


c
c
c
c
c
c
      subroutine sprime(srcinfo,targinfo,dpars,zk,ipars,u)
      implicit real *8 (a-h,o-z)
      integer ipars
      real *8 dpars,srcinfo(2),targinfo(4)
      complex *16 u,h0,ima,zs,z,zk,h1,gx,gy,ztmp
      data ima/(0.0d0,1.0d0)/
      data zs/(0.0d0,0.25d0)/

      rr2 = (srcinfo(1)-targinfo(1))**2 + (srcinfo(2)-targinfo(2))**2
      rr = dsqrt(rr2)
      rinv = 1/rr
      z = zk*rr
      ifexpon = 1
      call hank103(z,h0,h1,ifexpon)

      ztmp = -zs*zk*h1*rinv

      gx = ztmp*(targinfo(1)-srcinfo(1))
      gy = ztmp*(targinfo(2)-srcinfo(2))
      
      u = gx*targinfo(3) + gy*targinfo(4)

      return
      end
c
c
c
c
c
      subroutine dlp(srcinfo,targinfo,dpars,zk,ipars,u)
      implicit real *8 (a-h,o-z)
      integer ipars
      real *8 dpars,srcinfo(4),targinfo(2)
      complex *16 u,h0,ima,zs,z,zk,h1,gx,gy,ztmp
      data ima/(0.0d0,1.0d0)/
      data zs/(0.0d0,0.25d0)/

      rr2 = (srcinfo(1)-targinfo(1))**2 + (srcinfo(2)-targinfo(2))**2
      rr = dsqrt(rr2)
      rinv = 1/rr
      z = zk*rr
      ifexpon = 1
      call hank103(z,h0,h1,ifexpon)

      ztmp = -zs*zk*h1*rinv

      gx = ztmp*(srcinfo(1)-targinfo(1))
      gy = ztmp*(srcinfo(2)-targinfo(2))
      
      u = gx*srcinfo(3) + gy*srcinfo(4)

      return
      end
c
c
c
c
c
      subroutine comb(srcinfo,targinfo,dpars,zpars,ipars,u)
      implicit real *8 (a-h,o-z)
      integer ipars
      real *8 dpars,srcinfo(4),targinfo(2)
      complex *16 zpars(3),u,h0,ima,zs,z,zk,h1,gx,gy,ztmp
      data ima/(0.0d0,1.0d0)/
      data zs/(0.0d0,0.25d0)/

      zk = zpars(1)
      ima = dcmplx(0.0d0,1.0d0)
      zs = ima/4.0d0
      rr2 = (srcinfo(1)-targinfo(1))**2 + (srcinfo(2)-targinfo(2))**2
      rr = dsqrt(rr2)
      rinv = 1/rr
      z = zk*rr
      ifexpon = 1
      call hank103(z,h0,h1,ifexpon)

      ztmp = -zs*zk*h1*rinv

      gx = ztmp*(srcinfo(1)-targinfo(1))
      gy = ztmp*(srcinfo(2)-targinfo(2))
      
      u = zpars(3)*(gx*srcinfo(3) + gy*srcinfo(4)) + zpars(2)*zs*h0


      return
      end
c
c
c        transmission kernels
c
c
      subroutine transmission_dir(srcinfo,targinfo,dpars,zpars,ipars,u)
c
c
c         The kernel of interaction is given by
c           alpha S_{k1} + beta S_{k2} + gamma D_{k1} + delta D_{k2}
c         
c          zpars(1) = k1
c          zpars(2) = k2
c          zpars(3:6) = alpha,beta,gamma,delta
c          
c

      implicit real *8 (a-h,o-z)
      integer ipars
      real *8 dpars,srcinfo(4),targinfo(2)
      complex *16 zpars(6),u,h0,ima,zs,z,zk,h1,gx,gy,ztmp
      data ima/(0.0d0,1.0d0)/
      data zs/(0.0d0,0.25d0)/

      zk = zpars(1)
      zk2 = zpars(2)
      
      rr2 = (srcinfo(1)-targinfo(1))**2 + (srcinfo(2)-targinfo(2))**2
      rr = dsqrt(rr2)

      rinv = 1/rr


      z = zk*rr
      ifexpon = 1
      call hank103(z,h0,h1,ifexpon)

      ztmp = -zs*zk*h1*rinv

      gx = ztmp*(srcinfo(1)-targinfo(1))
      gy = ztmp*(srcinfo(2)-targinfo(2))
      
      u = zpars(5)*(gx*srcinfo(3) + gy*srcinfo(4)) + zpars(3)*zs*h0

      z = zk2*rr
      ifexpon = 1
      call hank103(z,h0,h1,ifexpon)

      ztmp = -zs*zk2*h1*rinv

      gx = ztmp*(srcinfo(1)-targinfo(1))
      gy = ztmp*(srcinfo(2)-targinfo(2))
      
      u = u+zpars(6)*(gx*srcinfo(3) + gy*srcinfo(4)) + zpars(4)*zs*h0


      return
      end
c
c
c
c
c
      subroutine transmission_neu(srcinfo,targinfo,dpars,zpars,ipars,u)
c
c
c         The kernel of interaction is given by
c           alpha S_{k1}' + beta S_{k2}' + gamma D_{k1}' + delta D_{k2}'
c         
c          zpars(1) = k1
c          zpars(2) = k2
c          zpars(3:6) = alpha,beta,gamma,delta
c          
c

      implicit real *8 (a-h,o-z)
      integer ipars
      real *8 dpars,srcinfo(4),targinfo(4)
      complex *16 zpars(6),u,h0,ima,zs,z,zk,h1,gx,gy,h2,zk2
      complex *16 d2gdx2,d2gdy2,d2gdxdy,ztmp
      complex *16 gd0,gs0,gd1,gs1
      data ima/(0.0d0,1.0d0)/
      data zs/(0.0d0,0.25d0)/

      zk = zpars(1)
      zk2 = zpars(2)

      xd = targinfo(1) - srcinfo(1)
      yd = targinfo(2) - srcinfo(2)
      
      rr2 = (srcinfo(1)-targinfo(1))**2 + (srcinfo(2)-targinfo(2))**2
      rr = dsqrt(rr2)
      rinv = 1/rr

      xd = xd*rinv
      yd = yd*rinv


      z = zk*rr
      ifexpon = 1
      call hank103(z,h0,h1,ifexpon)
      h2 = (2/z*h1 - h0)*zk

      d2gdx2 = (-h1*rinv + h2*xd*xd)*zk*zs
      d2gdxdy = h2*xd*yd*zk*zs
      d2gdy2 = (-h1*rinv+h2*yd*yd)*zk*zs

      gd0 = -(d2gdx2*srcinfo(3)*targinfo(3) +
     1    d2gdxdy*(srcinfo(3)*targinfo(4) + srcinfo(4)*targinfo(3)) + 
     2    d2gdy2*srcinfo(4)*targinfo(4))

      gx = -zs*zk*h1*(targinfo(1)-srcinfo(1))*rinv
      gy = -zs*zk*h1*(targinfo(2)-srcinfo(2))*rinv

      gs0 = gx*targinfo(3) + gy*targinfo(4)


      u = zpars(3)*gs0 + zpars(5)*gd0

      z = zk2*rr
      ifexpon = 1
      call hank103(z,h0,h1,ifexpon)
      h2 = (2/z*h1 - h0)*zk2

      d2gdx2 = (-h1*rinv + h2*xd*xd)*zk2
      d2gdxdy = h2*xd*yd*zk2
      d2gdy2 = (-h1*rinv+h2*yd*yd)*zk2

      gd1 = -zs*(d2gdx2*srcinfo(3)*targinfo(3) +
     1    d2gdxdy*(srcinfo(3)*targinfo(4) + srcinfo(4)*targinfo(3)) + 
     2    d2gdy2*srcinfo(4)*targinfo(4))

      gx = -zs*zk2*h1*(targinfo(1)-srcinfo(1))*rinv
      gy = -zs*zk2*h1*(targinfo(2)-srcinfo(2))*rinv
      
      gs1 = gx*targinfo(3) + gy*targinfo(4)


      u = u+zpars(4)*gs1 + zpars(6)*gd1


      return
      end
c
c
c
c
c
c
      subroutine dprime0_diff(srcinfo,targinfo,dpars,zpars,ipars,u)
c
c
c         The kernel of interaction is given by
c           D_{k1}' - D_{0}'
c         
c          zpars(1) = k1
c          
c

      implicit real *8 (a-h,o-z)
      integer ipars
      real *8 dpars,srcinfo(4),targinfo(4),over4pi
      complex *16 zpars(1),u,h0,ima,zs,z,zk,h1,gx,gy,h2,zk2
      complex *16 d2gdx2,d2gdy2,d2gdxdy,ztmp
      complex *16 gd0,gs0,gd1,gs1
      data ima/(0.0d0,1.0d0)/
      data zs/(0.0d0,0.25d0)/
      data over2pi/0.15915494309189535d0/

      zk = zpars(1)

      xd = targinfo(1) - srcinfo(1)
      yd = targinfo(2) - srcinfo(2)
      
      rr2 = (srcinfo(1)-targinfo(1))**2 + (srcinfo(2)-targinfo(2))**2
      rr = dsqrt(rr2)
      rinv = 1/rr

      xd = xd*rinv
      yd = yd*rinv


      z = zk*rr
      ifexpon = 1
      call hank103(z,h0,h1,ifexpon)
      h2 = (2/z*h1 - h0)*zk

      d2gdx2 = (-h1*rinv + h2*xd*xd)*zk*zs
      d2gdxdy = h2*xd*yd*zk*zs
      d2gdy2 = (-h1*rinv+h2*yd*yd)*zk*zs

      gd0 = -(d2gdx2*srcinfo(3)*targinfo(3) +
     1    d2gdxdy*(srcinfo(3)*targinfo(4) + srcinfo(4)*targinfo(3)) + 
     2    d2gdy2*srcinfo(4)*targinfo(4))

      u = gd0

      d2gdx2 = rinv**2*(1.0d0 - 2.0d0*xd*xd)
      d2gdy2 = rinv**2*(1.0d0 - 2.0d0*yd*yd)
      d2gdxdy = -rinv**2*xd*yd*2.0d0

      gd1 = over2pi*(d2gdx2*srcinfo(3)*targinfo(3) +
     1    d2gdxdy*(srcinfo(3)*targinfo(4) + srcinfo(4)*targinfo(3)) + 
     2    d2gdy2*srcinfo(4)*targinfo(4))

      u = u-gd1



      return
      end
c
c
c
c
c
c
      subroutine helm_c_p(n,srcinfo,m,targinfo,zk,u)
c
c       single layer interaction kernel
c
c       input:
c         srcinfo(2,n) - double
c           x,y location of source
c         targinfo(2,m) - double
c           x,y location of target
c         zk - complex
c           Helmholtz parameter
c
c       output:
c         u = i/4 H_{0}(zk*|r|) \, ,
c         where r is the distance between source and target
c          
c
      implicit real *8 (a-h,o-z)
      integer ipars
      real *8 dpars,srcinfo(2,n),targinfo(2,m)
      complex *16 u(m,n),h0,ima,zs,z,zk,h1
      data ima/(0.0d0,1.0d0)/
      data zs/(0.0d0,0.25d0)/

      ifexpon = 1
      do i=1,n
        do j=1,m
          rr2 = (srcinfo(1,i)-targinfo(1,j))**2 + 
     1        (srcinfo(2,i)-targinfo(2,j))**2
          rr = dsqrt(rr2)
          z = zk*rr

          call hank103(z,h0,h1,ifexpon)
          u(j,i) = zs*h0
        enddo
      enddo
      

      return
      end


c
c
c
c
c
c
      subroutine helm_c_gn(n,srcinfo,m,targinfo,zk,u)
      implicit real *8 (a-h,o-z)
      integer ipars
      real *8 dpars,srcinfo(2,n),targinfo(4,m)
      complex *16 u(m,n),h0,ima,zs,z,zk,h1,gx,gy,ztmp
      data ima/(0.0d0,1.0d0)/
      data zs/(0.0d0,0.25d0)/


      ifexpon = 1
      do i=1,n
        do j=1,m
           rr2 = (srcinfo(1,i)-targinfo(1,j))**2 + 
     1         (srcinfo(2,i)-targinfo(2,j))**2
           rr = dsqrt(rr2)
           rinv = 1/rr
           z = zk*rr
           call hank103(z,h0,h1,ifexpon)

           ztmp = -zs*zk*h1*rinv

           gx = ztmp*(targinfo(1,j)-srcinfo(1,i))
           gy = ztmp*(targinfo(2,j)-srcinfo(2,i))
      
           u(j,i) = gx*targinfo(3,j) + gy*targinfo(4,j)
         enddo
       enddo

      return
      end
c
c
c
c
c
c
c
      subroutine helm_c_g(n,srcinfo,m,targinfo,zk,ux,uy)
      implicit real *8 (a-h,o-z)
      integer ipars
      real *8 dpars,srcinfo(2,n),targinfo(2,m)
      complex *16 ux(m,n),uy(m,n),h0,ima,zs,z,zk,h1,gx,gy,ztmp
      data ima/(0.0d0,1.0d0)/
      data zs/(0.0d0,0.25d0)/


      ifexpon = 1
      do i=1,n
        do j=1,m
           rr2 = (srcinfo(1,i)-targinfo(1,j))**2 + 
     1         (srcinfo(2,i)-targinfo(2,j))**2
           rr = dsqrt(rr2)
           rinv = 1/rr
           z = zk*rr
           call hank103(z,h0,h1,ifexpon)

           ztmp = -zs*zk*h1*rinv

           ux(j,i) = ztmp*(targinfo(1,j)-srcinfo(1,i))
           uy(j,i) = ztmp*(targinfo(2,j)-srcinfo(2,i))
      
         enddo
       enddo

      return
      end
c
c
c
c
c
      subroutine helm_d_p(n,srcinfo,m,targinfo,zk,u)
      implicit real *8 (a-h,o-z)
      integer ipars
      real *8 dpars,srcinfo(4,n),targinfo(2,m)
      complex *16 u(m,n),h0,ima,zs,z,zk,h1,gx,gy,ztmp
      data ima/(0.0d0,1.0d0)/
      data zs/(0.0d0,0.25d0)/

      ifexpon = 1
      do i=1,n
        do j=1,m
          rr2 = (srcinfo(1,i)-targinfo(1,j))**2 + 
     1        (srcinfo(2,i)-targinfo(2,j))**2
          rr = dsqrt(rr2)
          rinv = 1/rr
          z = zk*rr
          call hank103(z,h0,h1,ifexpon)

          ztmp = -zs*zk*h1*rinv

          gx = ztmp*(srcinfo(1,i)-targinfo(1,j))
          gy = ztmp*(srcinfo(2,i)-targinfo(2,j))
      
          u(j,i) = gx*srcinfo(3,i) + gy*srcinfo(4,i)
        enddo
      enddo

      return
      end
c
c
c
c
c
      subroutine helm_d_gn(n,srcinfo,m,targinfo,zk,u)
      implicit real *8 (a-h,o-z)
      integer ipars
      real *8 dpars,srcinfo(4,n),targinfo(4,m),xd,yd
      complex *16 u(m,n),h0,ima,zs,z,zk,h1,gx,gy,ztmp
      complex *16 h2,d2gdx2,d2gdxdy,d2gdy2,dgn
      data ima/(0.0d0,1.0d0)/
      data zs/(0.0d0,0.25d0)/


      ifexpon = 1
      do i=1,n
        do j=1,m
           xd = targinfo(1,j) - srcinfo(1,i)
           yd = targinfo(2,j) - srcinfo(2,i)
           rr2 = (srcinfo(1,i)-targinfo(1,j))**2 + 
     1         (srcinfo(2,i)-targinfo(2,j))**2
           rr = dsqrt(rr2)
           rinv = 1/rr
           z = zk*rr
           call hank103(z,h0,h1,ifexpon)
           h2 = (2/z*h1 - h0)*zk

           d2gdx2 = (-h1*rinv + h2*xd*xd)*zk*zs
           d2gdxdy = h2*xd*yd*zk*zs
           d2gdy2 = (-h1*rinv + h2*yd*yd)*zk*zs

           dgn = -(d2gdx2*srcinfo(3,i)*targinfo(3,j) + 
     1       d2gdxdy*(srcinfo(3,i)*targinfo(4,j) +
     2       srcinfo(4,i)*targinfo(3,j)) + 
     3       d2gdy2*srcinfo(4,i)*targinfo(4,j))

      
           u(j,i) = dgn 
         enddo
       enddo

      return
      end
c
c
c
c
c
