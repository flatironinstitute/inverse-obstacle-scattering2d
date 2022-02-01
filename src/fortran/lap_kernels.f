      subroutine lap_dlp(srcinfo,targinfo,dpars,zpars,ipars,u)
      implicit real *8 (a-h,o-z)
      integer ipars
      real *8 dpars
      real *8 srcinfo(4),targinfo(2)
      complex *16 u,zpars
      data over2pi/0.15915494309189535d0/

      
      rr2 = (srcinfo(1)-targinfo(1))**2 + (srcinfo(2)-targinfo(2))**2

      xd = targinfo(1) - srcinfo(1)
      yd = targinfo(2) - srcinfo(2)

      gx = xd/rr2 
      gy = yd/rr2

      u = over2pi*(gx*srcinfo(3) + gy*srcinfo(4))


      return
      end
c
c
c
c
c
c
      subroutine lap_c_p(n,srcinfo,m,targinfo,u)
c
c       single layer interaction kernel
c
c       input:
c         srcinfo(2,n) - double
c           x,y location of source
c         targinfo(2,m) - double
c           x,y location of target
c
c       output:
c         u = -1/2pi log(|r|) \, ,
c         where r is the distance between source and target
c          
c
      implicit real *8 (a-h,o-z)
      integer ipars
      real *8 dpars,srcinfo(2,n),targinfo(2,m),over2pi
      complex *16 u(m,n)
      data over2pi/0.15915494309189535d0/

      do i=1,n
        do j=1,m
          rr2 = (srcinfo(1,i)-targinfo(1,j))**2 + 
     1        (srcinfo(2,i)-targinfo(2,j))**2
          u(j,i) = -log(rr2)/2*over2pi 
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
      subroutine lap_c_gn(n,srcinfo,m,targinfo,u)
      implicit real *8 (a-h,o-z)
      integer ipars
      real *8 dpars,srcinfo(2,n),targinfo(4,m)
      complex *16 u(m,n)
      data over2pi/0.15915494309189535d0/


      do i=1,n
        do j=1,m
           rr2 = (srcinfo(1,i)-targinfo(1,j))**2 + 
     1         (srcinfo(2,i)-targinfo(2,j))**2
           dx = targinfo(1,j)-srcinfo(1,i)
           dy = targinfo(2,j)-srcinfo(2,i)
      
           u(j,i) = -over2pi*(dx/rr2*targinfo(3,j) + 
     1         dy/rr2*targinfo(4,j))
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
      subroutine lap_d_p(n,srcinfo,m,targinfo,u)
      implicit real *8 (a-h,o-z)
      integer ipars
      real *8 dpars,srcinfo(4,n),targinfo(2,m)
      complex *16 u(m,n)
      data over2pi/0.15915494309189535d0/


      do i=1,n
        do j=1,m
           rr2 = (srcinfo(1,i)-targinfo(1,j))**2 + 
     1         (srcinfo(2,i)-targinfo(2,j))**2
           dx = targinfo(1,j)-srcinfo(1,i)
           dy = targinfo(2,j)-srcinfo(2,i)
      
           u(j,i) = over2pi*(dx/rr2*srcinfo(3,i) + 
     1         dy/rr2*srcinfo(4,i))
         enddo
       enddo

      return
      end
c
c
c
c
c
