c
c
c     this file contains various routines for forming 
c      different matrices that arise in solving 
c      helmholtz dirichlet/Neumann/transmission problems
c
c      slp_mat - single layer matrix
c      dlp_ext_mat - double layer exterior limit matrix
c      comb_ext_mat - combined field exterior limit matrix
c      trans_mat - matrix corresponding to the transmission 
c                  problem
c      
c
      subroutine slp_mat(n,norder,h,srcinfo,zk,xmat)
c
c        this subroutine constructs the discretization
c        matrix corresponding to the single layer
c
c        input:
c         n - number of discretization points
c         norder - order of alpert quadrature
c         h - 2*pi/n
c           spacing in parameter space
c         srcinfo - double (5,n)
c           boundary discretization info
c           srcinfo(1:2,:) - xy locations
c           srcinfo(3:4,:) - normal info
c           srcinfo(5,:) - dsdt
c           srcinfo(6,:) - curvature
c         zk - complex *16
c           Helmholtz parameter
c 
c        output:
c          xmat(n,n) - complex *16
c            helmholtz single layer matrix

      implicit none
      integer n,norder
      real *8 srcinfo(6,n),dpars,h
      integer ipars
      complex *16 zk,xmat(n,n)
      external slp

      call formmatbac(xmat,norder,n,srcinfo,h,slp,dpars,zk,ipars)


      return
      end
c
c
c
c
      subroutine dlp_ext_mat(n,norder,h,srcinfo,zk,xmat)
c
c        this subroutine constructs the discretization
c        matrix corresponding to the exterior limit
c        of the double layer
c
c        input:
c         n - number of discretization points
c         norder - order of alpert quadrature
c         h - 2*pi/n
c           spacing in parameter space
c         srcinfo - double (5,n)
c           boundary discretization info
c           srcinfo(1:2,:) - xy locations
c           srcinfo(3:4,:) - normal info
c           srcinfo(5,:) - dsdt
c           srcinfo(6,:) - curvature
c         zk - complex *16
c           Helmholtz parameter
c 
c        output:
c          xmat(n,n) - complex *16
c            exterior limit of Helmholtz double layer matrix

      implicit none
      integer n,i,norder
      real *8 srcinfo(6,n),dpars,h
      integer ipars
      complex *16 zk,xmat(n,n)
      external dlp

      call formmatbac(xmat,norder,n,srcinfo,h,dlp,dpars,zk,ipars)

      do i=1,n
        xmat(i,i) = xmat(i,i) + 0.5d0
      enddo


      return
      end
c
c
c
c

      subroutine comb_ext_mat(n,norder,h,srcinfo,zpars,xmat)
c
c        this subroutine constructs the discretization
c        matrix corresponding to the exterior limit
c        of the combined field representation
c
c        input:
c         n - number of discretization points
c         norder - order of alpert quadrature
c         h - 2*pi/n
c           spacing in parameter space
c         srcinfo - double (5,n)
c           boundary discretization info
c           srcinfo(1:2,:) - xy locations
c           srcinfo(3:4,:) - normal info
c           srcinfo(5,:) - dsdt
c           srcinfo(6,:) - curvature
c         zpars(3) - complex *16
c           zpars(1) = zk = Helmholtz parameter
c           zpars(2) - single layer strength
c           zpars(3) - double layer strength
c 
c        output:
c          xmat(n,n) - complex *16
c            helmholtz combined field exterior limit 

      implicit none
      integer n,norder,i
      real *8 srcinfo(6,n),dpars,h
      integer ipars
      complex *16 zpars(3),xmat(n,n)
      external comb

      call formmatbac(xmat,norder,n,srcinfo,h,comb,dpars,zpars,ipars)

      do i=1,n
        xmat(i,i) = xmat(i,i) + 0.5d0*zpars(3)
      enddo
      

      return
      end
c
c
c
c
c
c

      subroutine sprime_ext_mat(n,norder,h,srcinfo,zpars,xmat)
c
c        this subroutine constructs the discretization
c        matrix corresponding to the exterior limit
c        of the neumann data for the single layer potential 
c
c        input:
c         n - number of discretization points
c         norder - order of alpert quadrature
c         h - 2*pi/n
c           spacing in parameter space
c         srcinfo - double (5,n)
c           boundary discretization info
c           srcinfo(1:2,:) - xy locations
c           srcinfo(3:4,:) - normal info
c           srcinfo(5,:) - dsdt
c           srcinfo(6,:) - curvature
c         zpars(1) - complex *16
c           zpars(1) = zk = Helmholtz parameter
c 
c        output:
c          xmat(n,n) - complex *16
c            helmholtz combined field exterior limit 

      implicit none
      integer n,norder,i
      real *8 srcinfo(6,n),dpars,h
      integer ipars
      complex *16 zpars,xmat(n,n)
      external sprime

      call formmatbac(xmat,norder,n,srcinfo,h,sprime,dpars,zpars,ipars)

      do i=1,n
        xmat(i,i) = xmat(i,i) - 0.5d0
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
c
c
c

      subroutine ddiff_neu_mat(n,norder,h,srcinfo,zpars,xmat)
c
c        this subroutine constructs the discretization
c        matrix corresponding to the exterior limit
c        of the neumann data for the difference of two double
c        layer potentials
c
c        D_{k1}-D_{k2}
c
c        input:
c         n - number of discretization points
c         norder - order of alpert quadrature
c         h - 2*pi/n
c           spacing in parameter space
c         srcinfo - double (5,n)
c           boundary discretization info
c           srcinfo(1:2,:) - xy locations
c           srcinfo(3:4,:) - normal info
c           srcinfo(5,:) - dsdt
c           srcinfo(6,:) - curvature
c         zpars(2) - complex *16
c           zpars(1) = k1 = Helmholtz parameter
c           zpars(2) = k2 = Helmholtz parameter 2 
c 
c        output:
c          xmat(n,n) - complex *16
c            helmholtz combined field exterior limit 

      implicit none
      integer n,norder,i
      real *8 srcinfo(6,n),dpars,h
      integer ipars
      complex *16 zpars(2),xmat(n,n)
      complex *16 zpars_tmp(6)
      external transmission_neu 


      zpars_tmp(1) = zpars(1)
      zpars_tmp(2) = zpars(2)
      zpars_tmp(3) = 0
      zpars_tmp(4) = 0
      zpars_tmp(5) = 1
      zpars_tmp(6) = -1
      call formmatbac(xmat,norder,n,srcinfo,h,transmission_neu,
     1   dpars,zpars_tmp,ipars)


      return
      end
c
c
c
c
c

      subroutine ddiff0_neu_mat(n,norder,h,srcinfo,zk,xmat)
c
c        this subroutine constructs the discretization
c        matrix corresponding to the exterior limit
c        of the neumann data for the difference 
c
c        D_{k}'-D_{0}'
c
c        input:
c         n - number of discretization points
c         norder - order of alpert quadrature
c         h - 2*pi/n
c           spacing in parameter space
c         srcinfo - double (5,n)
c           boundary discretization info
c           srcinfo(1:2,:) - xy locations
c           srcinfo(3:4,:) - normal info
c           srcinfo(5,:) - dsdt
c           srcinfo(6,:) - curvature
c         zk - complex *16
c           k = Helmholtz parameter
c 
c        output:
c          xmat(n,n) - complex *16
c            helmholtz combined field exterior limit 

      implicit none
      integer n,norder,i
      real *8 srcinfo(6,n),dpars,h
      integer ipars
      complex *16 zk,xmat(n,n)
      external dprime0_diff 

      call formmatbac(xmat,norder,n,srcinfo,h,dprime0_diff,
     1   dpars,zk,ipars)

      return
      end
c
c
c
c

      subroutine lap_dlp_mat(n,norder,h,srcinfo,xmat)
c
c        this subroutine constructs the discretization
c        matrix corresponding to the exterior limit
c        of the dirichlet data for 
c       
c        D_{0}
c
c        input:
c         n - number of discretization points
c         norder - order of alpert quadrature
c         h - 2*pi/n
c           spacing in parameter space
c         srcinfo - double (5,n)
c           boundary discretization info
c           srcinfo(1:2,:) - xy locations
c           srcinfo(3:4,:) - normal info
c           srcinfo(5,:) - dsdt
c           srcinfo(6,:) - curvature
c 
c        output:
c          xmat(n,n) - complex *16
c            Laplace double layer discretization 

      implicit none
      integer n,norder,i
      real *8 srcinfo(6,n),dpars,h
      integer ipars
      complex *16 xmat(n,n),zpars
      external lap_dlp 

      call formmatbac(xmat,norder,n,srcinfo,h,lap_dlp,
     1   dpars,zpars,ipars)

      do i=1,n
        xmat(i,i) = xmat(i,i) + 0.5d0
      enddo


      return
      end
c
c
c
c
c
c
      subroutine trans_mat(n,norder,h,srcinfo,zks,a,b,xmat)
c
c        this subroutine constructs the discretization
c        matrix corresponding to the following tranmission
c        problem
c
c        subscript 1 - interior
c        subscript 2 - exterior
c
c        [au]/q = f/q, [b du/dn] = g \, ,
c
c        The representation chosen for u is
c        u_{1} = -(1/b_{1}) S_{k_{1}} \sigma + (1/b_{1}) D_{k_{1}} \mu 
c
c        u_{2} = -(1/b_{2}) S_{k_{2}} \sigma + (1/b_{2}) D_{k_{2}} \mu
c 
c        and the constant q is given by
c        q = 0.5d0*(a_{1}/b_{1} + a_{2}/b_{2})
c
c        inputs are ordered as, sigma_{1},\mu_{1}, \sigma_{2},\mu_{2}...
c        outputs are ordered as [au]/q_{1}, [b du/dn]_{1}, [au]/q_{2},
c            [b du/dn]_{2}...
c
c
c        input:
c         n - number of discretization points
c         norder - order of alpert quadrature
c         h - 2*pi/n
c           spacing in parameter space
c         srcinfo - double (5,n)
c           boundary discretization info
c           srcinfo(1:2,:) - xy locations
c           srcinfo(3:4,:) - normal info
c           srcinfo(5,:) - dsdt
c           srcinfo(6,:) - curvature
c         zks(2) - complex *16
c           zks(1) = interior Helmholtz parameter
c           zks(2) - exterior Helmholtz parameter
c         a(2) - complex *16
c           scaling for jump in u
c         b(2) - complex *16
c           scaling for jump in du/dn
c 
c        output:
c          xmat(2*n,2*n) - complex *16
c            Helmholtz transmission matrix
c
c
      implicit none
      integer n,i,j,nsys,norder
      real *8 srcinfo(6,n),dpars,h
      integer ipars
      complex *16 zpars(6),zks(2),a(2),b(2),xmat(2*n,2*n),q
      complex *16, allocatable :: xmattmp(:,:)

      external transmission_dir,transmission_neu

      allocate(xmattmp(n,n))

      nsys = 2*n
      do i=1,nsys
        do j=1,nsys
          xmat(i,j) = 0
        enddo
      enddo

c
c
c        the unknowns are organized as
c    \sigma_{1},\mu_{1},\sigma_{2},\mu_{2}...
c

c
c     form (a/b sk - a0/b0 sk0)/q system matrix (k is interior k0 is 
c      exterior
c
      q = 0.5d0*(a(1)/b(1) + a(2)/b(2))
      zpars(1) = zks(1)
      zpars(2) = zks(2)
      zpars(3) = a(1)/b(1)/q
      zpars(4) = -a(2)/b(2)/q
      zpars(5) = 0
      zpars(6) = 0

      call formmatbac(xmattmp,norder,n,srcinfo,h,transmission_dir,
     1    dpars,zpars,ipars)

  
      do i=1,n
        do j=1,n
          xmat(2*i-1,2*j-1) = xmattmp(i,j) 
        enddo
      enddo

c
c     form (-a/b dk + a0/b0 dk0)/q system matrix (k is interior k0 is 
c      exterior
c
      zpars(1) = zks(1)
      zpars(2) = zks(2)
      zpars(3) = 0
      zpars(4) = 0
      zpars(5) = -a(1)/b(1)/q
      zpars(6) = a(2)/b(2)/q

      call formmatbac(xmattmp,norder,n,srcinfo,h,transmission_dir,
     1    dpars,zpars,ipars)
 
      do i=1,n
        do j=1,n
          xmat(2*i-1,2*j) = xmattmp(i,j)
        enddo
        xmat(2*i-1,2*i) = xmat(2*i-1,2*i)+1 
      enddo

c
c     form sk' - sk0' system matrix (k is interior k0 is 
c      exterior
c
      zpars(1) = zks(1)
      zpars(2) = zks(2)
      zpars(3) = 1.0d0
      zpars(4) = -1.0d0
      zpars(5) = 0
      zpars(6) = 0


      call formmatbac(xmattmp,norder,n,srcinfo,h,transmission_neu,
     1    dpars,zpars,ipars)
  
      do i=1,n
        do j=1,n
          xmat(2*i,2*j-1) = xmattmp(i,j) 
        enddo
        xmat(2*i,2*i-1) = xmat(2*i,2*i-1)+1 
      enddo


c
c     form -dk + dk0 system matrix (k is interior k0 is 
c      exterior
c
      zpars(1) = zks(1)
      zpars(2) = zks(2)
      zpars(3) = 0
      zpars(4) = 0
      zpars(5) = -1
      zpars(6) = 1
      call formmatbac(xmattmp,norder,n,srcinfo,h,transmission_neu,
     1    dpars,zpars,ipars)
 
      do i=1,n
        do j=1,n
          xmat(2*i,2*j) = xmattmp(i,j) 
        enddo
      enddo


c
c       end of generating matrix
c


      return
      end
c
c
c
c
