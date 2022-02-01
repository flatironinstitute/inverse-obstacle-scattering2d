cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
C        This file contains four user-callable subroutines: corrand4,
c        corrand_norm_4, corrand_integer_knuth_4, corrand_one_integer_4.
c        Following is a brief description of the said four subroutines.
c
c   corrand4 - returns to the user a pseudo-random vector distributed 
c        uniformly on the interval [0,1]. The algorithm used by this
c        subroutine is a pompous one: constructs nine pseudo-random 
c        sequences using nine different congruential generators, 
c        averages them, and uses the obvious transformation to bring it 
c        back to the uniform distribution. 
c
c   corrand_norm_4 - returns to the user two random vectors y1, y2, 
c        each of which is distributed normally with the 
c        distribution density
c
c        1/sqrt(2 * \pi) * e^(y^2/2)                         (1)
c   
c   corrand_integer_knuth_4 - for a user-specified n, returns to 
c        the user a random permutation of n integers 1,2,...,n
c
c   corrand_one_integer_4 - for a user-specified n, returns to 
c        the user a random integer on the interval [1,n]
c
C 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c 
c 
        subroutine corrand_one_integer_4(n,i)
        implicit real *8 (a-h,o-z)
        dimension tt(10)
        save
c       
c        This subroutine returns to the user a random integer
c        number i on the interval [0,n]
c
        call corrand4(1,tt(6))
c
        i=tt(6)*n
        i=i+1
c        
        return
        end
c
c
c
c 
c 
        subroutine corrand_integer_knuth_4(n,ixs)
        implicit real *8 (a-h,o-z)
        save
        dimension ixs(n),tt(12)
c
c        This subroutine returns to the user a random permutation 
c        of n integer numbers 1,2,...,n.
c
c              Input parameters:
c
c  n - random numbers in the array ixs will be the distributed
c        uniformly on the interval [1,n]
c
c              Output parameters:
c
c  ixs - a pseudo-random permutation of length n
c
        call corrand4(11,tt)
        call corrand4(11,tt)
        call corrand4(11,tt)
        do 1200 i=1,n
c
        ixs(i)=i
 1200 continue
c
        done=1
        do 1400 i=1,n-1
c
        call corrand4(1,tt)
c
        k=n-i+1

cccc        call prinf('k=*',k,1)

        h=done/k
        j=tt(1)/h+1

cccc        call prinf('and j=*',j,1)
c
        jj=ixs(k)
        ixs(k)=ixs(j)
        ixs(j)=jj
 1400 continue
c
        return
        end
c
c
c
c
c
        SUBROUTINE corrand_norm_4(N,Y1,y2)
        IMPLICIT REAL *8 (A-H,O-Z)
        save
        DIMENSION Y1(1),y2(1)
c
c        This subroutine returns to the user two random
c        vectors y1, y2, each of which is distiibuted
c        normally with the distribution density
c
c        1/sqrt(2 * \pi) * e^(y^2/2)                         (1)
c   
c
c              Input parameters:
c
c  n - the number of elements to be returned in each of the 
c        arrays y1, y2
c
c              Output parameters:
c  y1, y2 - two pseudo-random arrays distributed normally 
c        with the distribution density (1)
c
c
c        . . . construct vectors of variables distributed
c              uniformly on the interval [0,1]
c
        call corrand4(N,Y1)
        call corrand4(N,Y2)
c
c       combine the variables y1, y2 converting them
c       into variables distributed normally (Box-Muller 
c       algorithm)
c
        done=1
        pi=atan(done)*4
        do 1400 i=1,n
c
        z1=sqrt(-2*log(y1(i)))*cos(2*pi*y2(i))
        z2=sqrt(-2*log(y1(i)))*sin(2*pi*y2(i))
c
        y1(i)=z1
        y2(i)=z2
 1400 continue
c
        return
        end

        SUBROUTINE corrand4(n,y)
        IMPLICIT REAL *8 (A-H,O-Z)
        save
        dimension js21(10),js22(10),js23(10),y(1),
     1      ias(10),ixks(10),ds(10)
ccc        data ias/7,17,31,37,43,53,61,67,79,91/
        data ias/97,17,31,37,43,53,61,67,79,91/
        data ifcalled/0/

        done=1
c
c       retrieve the prime numbers and conduct preliminary
c       randomization
c
        if(ifcalled .ne. 0) goto 1350
c
        call corrand_primes_4(js21,js22,js23)
c
        do 1200 i=1,10
c
        ixks(i)=js21(i)
 1200 continue
        ifcalled=1
c
        do 1300 i=1,100
        do 1250 j=1,9    
c
        call corrand_onestep_4(js23(j),js22(j),ias(j),ixks(j))

 1250 continue
 1300 continue
c
 1350 continue
c
        do 2000 i=1,n
c
        do 1400 j=1,3
        call corrand_onestep_4(js23(j),js22(j),ias(j),ixks(j))
c
        ds(j)=ixks(j)*done/(js23(j)-1)
 1400 continue
c
        d=ds(1)+ds(2)+ds(3)
        call corrand_comp_4(d,dd1)
c
        do 1600 j=4,6
        call corrand_onestep_4(js23(j),js22(j),ias(j),ixks(j))
c
        ds(j)=ixks(j)*done/(js23(j)-1)
 1600 continue
c
        d=ds(4)+ds(5)+ds(6)
        call corrand_comp_4(d,dd2)

c
        do 1800 j=7,9
        call corrand_onestep_4(js23(j),js22(j),ias(j),ixks(j))
c
        ds(j)=ixks(j)*done/(js23(j)-1)
 1800 continue
c
        d=ds(7)+ds(8)+ds(9)
        call corrand_comp_4(d,dd3)
c
        d=dd1+dd2+dd3
        call corrand_comp_4(d,y(i))
c
 2000 continue
c
        return
        end
c
c
c
c
c
        SUBROUTINE corrand_comp_4(x,rint)
        IMPLICIT REAL *8 (A-H,O-Z)
        save
c
        if (x .lt. 1) then
            rint=x**3/6
            return
        endif
c      
        if ( (x .ge. 1) .and. (x .le. 2) ) then
            rint= 0.75d0*x-(x-1.5d0)**3/3 - 0.625d0
            return
        endif
c
        rint= (x-3)**3/6 +1
c
        return
        end
c
c
c
c
c
        SUBROUTINE corrand_onestep_4(m,ic,ia,ixk)
        IMPLICIT REAL *8 (A-H,O-Z)
        save
c
        jj=ia*ixk+ic  
        j=jj/m
        ixk=jj-j*m
c
        return
        END
c
c
c
c
c
        SUBROUTINE corrand_primes_4(js21,js22,js23)
        IMPLICIT REAL *8 (A-H,O-Z)
        save
        DIMENSION is21(10),is22(10),is23(10),
     1      js21(1),js22(1),js23(1)
c
        data is21/9,19,21,55,61,69,105,111,121,129/
        data is22/3,17,27,33,57,87,105,113,117,123/
        data is23/15,21,27,37,61,69,135,147,157,159/
c
        i21=2**21
        i22=i21*2
        i23=i22*2
c
        do 1200 i=1,10
c
        js21(i)=i21-is21(i)
        js22(i)=i22-is22(i)
        js23(i)=i23-is23(i)
 1200 continue
c
        return
        end
