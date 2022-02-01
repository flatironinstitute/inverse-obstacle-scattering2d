        subroutine qua_labels_conn(iw)
        implicit real *8 (a-h,o-z)
        save
        character *1 lb(2),file1(20),anum1(8),quo,gnlb(4),
     1      gnl(3)
        character *8 dummy,anum8,file8
c 
        equivalence (file1,file8),(anum1,anum8)
        data lb/'l','b'/,quo/'"'/,gnl/'g','n','l'/
c 
c        convert the user-specified Fortran unit number to
c        character format
c 
        if( (iw .ge. 0) .and. (iw .le. 9) )  write(dummy,1100) iw
 1100 format(i1)
c 
        if( (iw .ge. 10) .and. (iw .le. 99) ) write(dummy,1200) iw
 1200 format(i2)
c 
        if( (iw .ge. 100) .and. (iw .le. 999) ) write(dummy,1300) iw
 1300 format(i3)
c 
        if( (iw .ge. 1000) .and. (iw .le. 9999) ) write(dummy,1400) iw
 1400 format(i4)
c 
 2000 format(1a8)
        read(dummy,2000) anum8
c 
c        construct the file name on which the Gnuplot instructions
c        are to be written
c 
        file1(1)=gnl(1)
        file1(2)=gnl(2)
        file1(3)=gnl(3)
  
        do 2200 i=1,6
        file1(i+3)=anum1(i)
 2200 continue
c 
c        open the fortran file with the unit 87 and name file8
c 
        iun=87
        open(unit=iun,file=file8)
c 
c       write the "load" instructuins to the file we just opened
c 
        if( (iw .ge. 0) .and. (iw .le. 9) )  write(iun,2410)
     1     quo,anum1(1),quo
 2410 format(4x,'  load ',a1,'gnlb',2a1)
c 
        if( (iw .ge. 0) .and. (iw .le. 9) )  write(iun,2510)
     1     quo,anum1(1),quo
 2510 format(4x,'  load ',a1,'gn',2a1)
c 
        if( (iw .ge. 10) .and. (iw .le. 99) )  write(iun,2420)
     1     quo,anum1(1),anum1(2),quo
 2420 format(4x,'  load ',a1,'gnlb',3a1)
  
  
        if( (iw .ge. 10) .and. (iw .le. 99) )  write(iun,2520)
     1     quo,anum1(1),anum1(2),quo
 2520 format(4x,'  load ',a1,'gn',3a1)
c 
        if( (iw .ge. 100) .and. (iw .le. 999) )  then
            write(iun,2430) quo,anum1(1),anum1(2),anum1(3),
     1          quo
c 
 2430 format(4x,'  load ',a1,'gnlb',4a1)
c 
            write(iun,2530) quo,anum1(1),anum1(2),anum1(3),
     1          quo
c 
 2530 format(4x,'  load ',a1,'gn',4a1)
c 
        endif
c 
c 
c 
        if( (iw .ge. 1000) .and. (iw .le. 9999) )  then
            write(iun,2440) quo,anum1(1),anum1(2),anum1(3),
     1          anum1(4),quo
c 
 2440 format(4x,'  load ',a1,'gnlb',5a1)
c 
            write(iun,2540) quo,anum1(1),anum1(2),anum1(3),
     1          anum1(4),quo
c 
 2540 format(4x,'  load ',a1,'gn',5a1)
c 
        endif
c 
        close(iun)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine zqua_labels_int(iw,zs,labs,n)
        implicit real *8 (a-h,o-z)
        save
        integer *4 labs(1)
        dimension zs(2,1)
        character *1 lb(2),file1(20),anum1(20),quo,gnlb(4)
        character *8 dummy,anum8,file8
c 
        equivalence (file1,file8),(anum1,anum8)
        data lb/'l','b'/,quo/'"'/,gnlb/'g','n','l','b'/
c 
c        convert the user-specified Fortran unit number to
c        character format
c 
        if( (iw .ge. 0) .and. (iw .le. 9) )  write(dummy,1100) iw
 1100 format(i1)
c 
        if( (iw .ge. 10) .and. (iw .le. 99) ) write(dummy,1200) iw
 1200 format(i2)
c 
        if( (iw .ge. 100) .and. (iw .le. 999) ) write(dummy,1300) iw
 1300 format(i3)
c 
        if( (iw .ge. 1000) .and. (iw .le. 9999) ) write(dummy,1400) iw
 1400 format(i4)
c 
 2000 format(1a8)
        read(dummy,2000) anum8
c 
c        construct the file name on which the Gnuplot instructions
c        are to be written
c 
        file1(1)=gnlb(1)
        file1(2)=gnlb(2)
        file1(3)=gnlb(3)
        file1(4)=gnlb(4)
  
        do 2200 i=1,6
        file1(i+4)=anum1(i)
 2200 continue
c 
c        open the fortran file with the unit 87 and name file8
c 
        iun=87
        open(unit=iun,file=file8)
c 
c       write the labels to the file we just opened
c 
        do 2600 i=1,n
c 
 2410 format(4x,'set label ',i4,2x,1a1,i1,1a1,' at ',
     1    e11.5,', ',e11.5)
c 
 2420 format(4x,'set label ',i4,2x,1a1,i2,1a1,' at ',
     1    e11.5,', ',e11.5)
c 
 2430 format(4x,'set label ',i4,2x,1a1,i3,1a1,' at ',
     1    e11.5,', ',e11.5)
c 
 2440 format(4x,'set label ',i4,2x,1a1,i4,1a1,' at ',
     1    e11.5,', ',e11.5)
c 
 2450 format(4x,'set label ',i4,2x,1a1,i5,1a1,' at ',
     1    e11.5,', ',e11.5)
c 
        if( (labs(i) .ge. 0) .and. (labs(i) .le. 9) )
     1      write(iun,2410) i,quo,labs(i),quo,zs(1,i),zs(2,i)
c 
        if( (labs(i) .ge. 10) .and. (labs(i) .le. 99) )
     1      write(iun,2420) i,quo,labs(i),quo,zs(1,i),zs(2,i)
c 
        if( (labs(i) .ge. 100) .and. (labs(i) .le. 999) )
     1      write(iun,2430) i,quo,labs(i),quo,zs(1,i),zs(2,i)
c 
        if( (labs(i) .ge. 1000) .and. (labs(i) .le. 9999) )
     1      write(iun,2440) i,quo,labs(i),quo,zs(1,i),zs(2,i)
c 
        if( (labs(i) .ge. 10000) .and. (labs(i) .le. 99999) )
     1      write(iun,2450) i,quo,labs(i),quo,zs(1,i),zs(2,i)
c 
c 
 2600 continue
  
  
  
cccc    set label 1 "look" at 0.4, 0.4
  
cccc    set label 2 "fool" at 0.6, 0.4
  
        close(iun)
  
        call qua_labels_conn(iw)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine qua_labels_int(iw,x,y,labs,n)
        implicit real *8 (a-h,o-z)
        save
        integer *4 labs(1)
        dimension x(1),y(1)
        character *1 lb(2),file1(20),anum1(20),quo,gnlb(4)
        character *8 dummy,anum8,file8
c 
        equivalence (file1,file8),(anum1,anum8)
        data lb/'l','b'/,quo/'"'/,gnlb/'g','n','l','b'/
c 
c        convert the user-specified Fortran unit number to
c        character format
c 
        if( (iw .ge. 0) .and. (iw .le. 9) )  write(dummy,1100) iw
 1100 format(i1)
c 
        if( (iw .ge. 10) .and. (iw .le. 99) ) write(dummy,1200) iw
 1200 format(i2)
c 
        if( (iw .ge. 100) .and. (iw .le. 999) ) write(dummy,1300) iw
 1300 format(i3)
c 
        if( (iw .ge. 1000) .and. (iw .le. 9999) ) write(dummy,1400) iw
 1400 format(i4)
c 
 2000 format(1a8)
        read(dummy,2000) anum8
c 
c        construct the file name on which the Gnuplot instructions
c        are to be written
c 
        file1(1)=gnlb(1)
        file1(2)=gnlb(2)
        file1(3)=gnlb(3)
        file1(4)=gnlb(4)
  
        do 2200 i=1,6
        file1(i+4)=anum1(i)
 2200 continue
c 
c        open the fortran file with the unit 87 and name file8
c 
        iun=87
        open(unit=iun,file=file8)
c 
c       write the labels to the file we just opened
c 
        do 2600 i=1,n
c 
 2410 format(4x,'set label ',i4,2x,1a1,i1,1a1,' at ',
     1    e11.5,', ',e11.5)
c 
 2420 format(4x,'set label ',i4,2x,1a1,i2,1a1,' at ',
     1    e11.5,', ',e11.5)
c 
 2430 format(4x,'set label ',i4,2x,1a1,i3,1a1,' at ',
     1    e11.5,', ',e11.5)
c 
 2440 format(4x,'set label ',i4,2x,1a1,i4,1a1,' at ',
     1    e11.5,', ',e11.5)
c 
 2450 format(4x,'set label ',i4,2x,1a1,i5,1a1,' at ',
     1    e11.5,', ',e11.5)
c 
        if( (labs(i) .ge. 0) .and. (labs(i) .le. 9) )
     1      write(iun,2410) i,quo,labs(i),quo,x(i),y(i)
c 
        if( (labs(i) .ge. 10) .and. (labs(i) .le. 99) )
     1      write(iun,2420) i,quo,labs(i),quo,x(i),y(i)
c 
        if( (labs(i) .ge. 100) .and. (labs(i) .le. 999) )
     1      write(iun,2430) i,quo,labs(i),quo,x(i),y(i)
c 
        if( (labs(i) .ge. 1000) .and. (labs(i) .le. 9999) )
     1      write(iun,2440) i,quo,labs(i),quo,x(i),y(i)
c 
        if( (labs(i) .ge. 10000) .and. (labs(i) .le. 99999) )
     1      write(iun,2450) i,quo,labs(i),quo,x(i),y(i)
c 
c 
 2600 continue
  
  
  
  
        call qua_labels_conn(iw)
  
        close(iun)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine quaplot(iw,x,y,n,itype1,title)
        implicit real *8 (a-h,o-z)
        save
        dimension x(1),y(1),x2(1),y2(1),x3(1),y3(1)
        character *1 title(1)
c 
        inumgr=1
c 
        call quaplo0(iw,x,y,x2,y2,x3,y3,n,n2,n3,
     1      title,inumgr,itype1,itype2,itype3)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine quaplot2(iw,x,y,n,itype1,x2,y2,n2,itype2,title)
        implicit real *8 (a-h,o-z)
        save
        dimension x(1),y(1),x2(1),y2(1),x3(1),y3(1)
        character *1 title(1)
c 
        inumgr=2
c 
        call quaplo0(iw,x,y,x2,y2,x3,y3,n,n2,n3,
     1      title,inumgr,itype1,itype2,itype3)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine quaplot3(iw,x,y,n,itype1,x2,y2,n2,itype2,
     1      x3,y3,n3,itype3,title)
        implicit real *8 (a-h,o-z)
        save
        dimension x(1),y(1),x2(1),y2(1),x3(1),y3(1)
        character *1 title(1)
c 
        inumgr=3
c 
        call quaplo0(iw,x,y,x2,y2,x3,y3,n,n2,n3,
     1      title,inumgr,itype1,itype2,itype3)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine quagraph(iw,x,y,n,itype1,title)
        implicit real *8 (a-h,o-z)
        save
        dimension x(1),y(1),x2(1),y2(1),x3(1),y3(1)
        character *1 title(1)
c 
        inumgr=1
c 
        call quagrap0(iw,x,y,x2,y2,x3,y3,n,n2,n3,
     1      title,inumgr,itype1,itype2,itype3)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine quagraph2(iw,x,y,n,itype1,x2,y2,n2,itype2,title)
        implicit real *8 (a-h,o-z)
        save
        dimension x(1),y(1),x2(1),y2(1),x3(1),y3(1)
        character *1 title(1)
c 
        inumgr=2
c 
        call quagrap0(iw,x,y,x2,y2,x3,y3,n,n2,n3,
     1      title,inumgr,itype1,itype2,itype3)
c 
        return
        end
c 
c 
c 
c 
c 
c 
        subroutine quagraph3(iw,x,y,n,itype1,x2,y2,n2,itype2,
     1      x3,y3,n3,itype3,title)
        implicit real *8 (a-h,o-z)
        save
        dimension x(1),y(1),x2(1),y2(1),x3(1),y3(1)
        character *1 title(1)
c 
        inumgr=3
c 
        call quagrap0(iw,x,y,x2,y2,x3,y3,n,n2,n3,
     1      title,inumgr,itype1,itype2,itype3)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine quagrap0(iw,x,y,x2,y2,x3,y3,n,n2,n3,
     1      title,inumgr,itype1,itype2,itype3)
        implicit real *8 (a-h,o-z)
        save
        dimension x(1),y(1),x2(1),y2(1),x3(1),y3(1)
        character *1 title(1),gn(2),file1(8),anum1(8),line(82),
     1      blank,quo,backslash,anum1c(8),ccc,file1c(8)
        character *8 dummy,anum8,file8,anum8c,file8c
c 
        equivalence (file1,file8),(anum1,anum8),(anum1c,anum8c)
        equivalence (file1c,file8c)
        data gn/'g','n'/,blank/' '/,quo/'"'/,backslash/'\\'/,ccc/'c'/
c 
c        convert the user-specified Fortran unit number to
c        character format
c 
        if( (iw .ge. 0) .and. (iw .le. 9) )  write(dummy,1100) iw
 1100 format(i1)
c 
        if( (iw .ge. 10) .and. (iw .le. 99) ) write(dummy,1200) iw
 1200 format(i2)
c 
        if( (iw .ge. 100) .and. (iw .le. 999) ) write(dummy,1300) iw
 1300 format(i3)
c 
        if( (iw .ge. 1000) .and. (iw .le. 9999) ) write(dummy,1400) iw
 1400 format(i4)
c 
 2000 format(1a8)
        read(dummy,2000) anum8
        anum8c=anum8

        if( (iw .ge. 0) .and. (iw .le. 9) ) anum1c(2)=ccc
        if( (iw .ge. 10) .and. (iw .le. 99) ) anum1c(3)=ccc
        if( (iw .ge. 100) .and. (iw .le. 999) ) anum1c(4)=ccc
        if( (iw .ge. 1000) .and. (iw .le. 9999) ) anum1c(5)=ccc
c 
c        construct the file name on which the Gnuplot instructions
c        are to be written
c 
        file1(1)=gn(1)
        file1(2)=gn(2)
        file1c(1)=gn(1)
        file1c(2)=gn(2)
        do 2200 i=1,6
        file1(i+2)=anum1(i)
        file1c(i+2)=anum1c(i)
 2200 continue
c 
c        open the fortran file with the unit 87 and name file8
c 
        iun=87
        open(unit=iun,file=file8)
c 
c 
        iunc=97
        open(unit=iunc,file=file8c)
c 
 2250 format('   # set terminal postscript',/,
     1    '   # set output "plot.ps"', /,
     2    '   # set size 0.76, 1.0 ')
c 
 2255 format('    set terminal postscript',/,
     1    '    set output "plot.ps"')
c 
        write(iun,2250)
        write(iunc,2255)
c 
c        generate the title for the plot
c 
        line(1)=blank
        line(2)=blank
        line(3)=quo
c 
        call quamesslen(title,nchar,line(4))
c 
        line(nchar+4)=quo
 2300 format('  set title ',80a1)
        write(iun,2300) (line(i),i=1,nchar+4)
        write(iunc,2300) (line(i),i=1,nchar+4)
c 
 2350 format('   show title')
        write(iun,2350)
        write(iunc,2350)
c 
c        write the instructions
c 
 2400 format(i6)
        iun2=iw+100000
        write (dummy,2400) iun2
        read(dummy,2000) anum8
c 
        do 2600 i=1,6
        file1(i+2)=anum1(i)
 2600 continue
c 
 2800 format('plot ','"',1a8,'"     ','notitle  with dots')
 2803 format('plot ','"',1a8,'"     ','notitle  with points')
 2805 format('plot ','"',1a8,'"     ','notitle  with lines')
c 
        if( (inumgr .eq. 1) .and. (itype1 .eq. 1) )
     1      write(iun,2800) file8
        if( (inumgr .eq. 1) .and. (itype1 .eq. 1) )
     1      write(iunc,2800) file8
c 
        if( (inumgr .eq. 1) .and. (itype1 .eq. 2) )
     1      write(iun,2803) file8
        if( (inumgr .eq. 1) .and. (itype1 .eq. 2) )
     1      write(iunc,2803) file8
c 
        if( (inumgr .eq. 1) .and. (itype1 .eq. 3) )
     1      write(iun,2805) file8
        if( (inumgr .eq. 1) .and. (itype1 .eq. 3) )
     1      write(iunc,2805) file8
  
 2830 format('plot ','"',1a8,'"     ','notitle  with dots, ',1a1)
 2831 format('plot ','"',1a8,'"     ','notitle  with points, ',1a1)
 2832 format('plot ','"',1a8,'"     ','notitle  with lines, ',1a1)
c 
        if( (inumgr .ne. 1) .and. (itype1 .eq. 1) )
     1      write(iun,2830) file8,backslash
        if( (inumgr .ne. 1) .and. (itype1 .eq. 1) )
     1      write(iunc,2830) file8,backslash
c 
        if( (inumgr .ne. 1) .and. (itype1 .eq. 2) )
     1       write(iun,2831) file8,backslash
        if( (inumgr .ne. 1) .and. (itype1 .eq. 2) )
     1       write(iunc,2831) file8,backslash
c 
        if( (inumgr .ne. 1) .and. (itype1 .eq. 3) )
     1      write(iun,2832) file8,backslash
        if( (inumgr .ne. 1) .and. (itype1 .eq. 3) )
     1      write(iunc,2832) file8,backslash
c 
c        write the first data file to be plotted
c 
        iun22=88
        open(unit=iun22,file=file8)
c 
        write(iun22,3000) (x(i),y(i),i=1,n)
c 
        close(iun22)
c 
c       if the user so requested - write the instructions for the
c       plotting the second data file
c 
        if(inumgr .eq. 1) close(iun)
        if(inumgr .eq. 1) close(iunc)
        if(inumgr .eq. 1) return
  
        iun2=iw+200000
        write (dummy,2400) iun2
        read(dummy,2000) anum8
c 
        do 2840 i=1,6
        file1(i+2)=anum1(i)
 2840 continue
c 
 2850 format('     ','"',1a8,'"     ','notitle  with dots')
 2851 format('     ','"',1a8,'"     ','notitle  with points')
 2852 format('     ','"',1a8,'"     ','notitle  with lines')
c 
        if( (inumgr .eq. 2) .and. (itype2 .eq. 1) )
     1      write(iun,2850) file8
        if( (inumgr .eq. 2) .and. (itype2 .eq. 1) )
     1      write(iunc,2850) file8
c 
        if( (inumgr .eq. 2) .and. (itype2 .eq. 2) )
     1      write(iun,2851) file8
        if( (inumgr .eq. 2) .and. (itype2 .eq. 2) )
     1      write(iunc,2851) file8
c 
        if( (inumgr .eq. 2) .and. (itype2 .eq. 3) )
     1      write(iun,2852) file8
        if( (inumgr .eq. 2) .and. (itype2 .eq. 3) )
     1      write(iunc,2852) file8
c 
c 
        if( (inumgr .eq. 3) .and. (itype2 .eq. 1) )
     1      write(iun,2855) file8,backslash
        if( (inumgr .eq. 3) .and. (itype2 .eq. 1) )
     1      write(iunc,2855) file8,backslash
c 
        if( (inumgr .eq. 3) .and. (itype2 .eq. 2) )
     1      write(iun,2856) file8,backslash
        if( (inumgr .eq. 3) .and. (itype2 .eq. 2) )
     1      write(iunc,2856) file8,backslash
c 
        if( (inumgr .eq. 3) .and. (itype2 .eq. 3) )
     1      write(iun,2857) file8,backslash
        if( (inumgr .eq. 3) .and. (itype2 .eq. 3) )
     1      write(iunc,2857) file8,backslash
c 
 2855 format('     ','"',1a8,'"     ','notitle  with dots, ',1a1)
 2856 format('     ','"',1a8,'"     ','notitle  with points, ',1a1)
 2857 format('     ','"',1a8,'"     ','notitle  with lines, ',1a1)
c 
c        write the second data file to be plotted
c 
        iun22=88
        open(unit=iun22,file=file8)
c 
        write(iun22,3000) (x2(i),y2(i),i=1,n2)
c 
        close(iun22)
c 
c       if the user so requested - write the instructions for the
c       plotting the third data file
c 
        if(inumgr .eq. 2) close(iun)
        if(inumgr .eq. 2) return
  
        iun2=iw+300000
        write (dummy,2400) iun2
        read(dummy,2000) anum8
c 
        do 2860 i=1,6
        file1(i+2)=anum1(i)
 2860 continue
c 
        if(itype3 .eq. 1) write(iun,2870) file8
        if(itype3 .eq. 2) write(iun,2871) file8
        if(itype3 .eq. 3) write(iun,2872) file8
        if(itype3 .eq. 1) write(iunc,2870) file8
        if(itype3 .eq. 2) write(iunc,2871) file8
        if(itype3 .eq. 3) write(iunc,2872) file8
c 
 2870 format('     ','"',1a8,'"     ','notitle  with dots')
 2871 format('     ','"',1a8,'"     ','notitle  with points')
 2872 format('     ','"',1a8,'"     ','notitle  with lines')
c 
c        write the third data file to be plotted
c 
        iun22=88
        open(unit=iun22,file=file8)
c 
 3000 format(2x,e11.5,2x,e11.5)
c 
        write(iun22,3000) (x3(i),y3(i),i=1,n3)
c 
        close(iun22)
        close(iun)
        close(iunc)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine quaplo0(iw,x,y,x2,y2,x3,y3,n,n2,n3,
     1      title,inumgr,itype1,itype2,itype3)
        implicit real *8 (a-h,o-z)
        save
        dimension x(1),y(1),x2(1),y2(1),x3(1),y3(1)
        character *1 title(1),gn(2),file1(8),anum1(8),line(82),
     1      blank,quo,backslash,anum1c(8),ccc,file1c(8)
        character *8 dummy,anum8,file8,anum8c,file8c
c 
        equivalence (file1,file8),(anum1,anum8),(anum1c,anum8c)
        equivalence (file1c,file8c)
        data gn/'g','n'/,blank/' '/,quo/'"'/,backslash/'\\'/,ccc/'c'/
c 
c        convert the user-specified Fortran unit number to
c        character format
c 
        if( (iw .ge. 0) .and. (iw .le. 9) )  write(dummy,1100) iw
 1100 format(i1)
c 
        if( (iw .ge. 10) .and. (iw .le. 99) ) write(dummy,1200) iw
 1200 format(i2)
c 
        if( (iw .ge. 100) .and. (iw .le. 999) ) write(dummy,1300) iw
 1300 format(i3)
c 
        if( (iw .ge. 1000) .and. (iw .le. 9999) ) write(dummy,1400) iw
 1400 format(i4)
c 
 2000 format(1a8)
        read(dummy,2000) anum8
        anum8c=anum8
c
        if( (iw .ge. 0) .and. (iw .le. 9) ) anum1c(2)=ccc
        if( (iw .ge. 10) .and. (iw .le. 99) ) anum1c(3)=ccc
        if( (iw .ge. 100) .and. (iw .le. 999) ) anum1c(4)=ccc
        if( (iw .ge. 1000) .and. (iw .le. 9999) ) anum1c(5)=ccc
c 
c        construct the file name on which the Gnuplot instructions
c        are to be written
c 
        file1(1)=gn(1)
        file1(2)=gn(2)
        file1c(1)=gn(1)
        file1c(2)=gn(2)
        do 2200 i=1,6
        file1(i+2)=anum1(i)
        file1c(i+2)=anum1c(i)
 2200 continue
c 
c        open the fortran file with the unit 87 and name file8
c 
        iun=87
        open(unit=iun,file=file8)
c 
        iunc=97
        open(unit=iunc,file=file8c)
c 
cc 2250 format('   # set terminal postscript',/,
cc     1    '   # set output "plot.ps"')
cc
 2250 format('   # set terminal postscript',/,
     1    '   # set output "plot.ps"', /,
     2    '   # set size 0.76, 1.0 ')
c 

c 
        write(iun,2250)
c 
c        generate the title for the plot
c 
c 
 2255 format('    set terminal postscript',/,
     1    '    set output "plot.ps"')
c 
        write(iunc,2255)

        line(1)=blank
        line(2)=blank
        line(3)=quo
c 
        call quamesslen(title,nchar,line(4))
c 
        line(nchar+4)=quo
 2270 format('  set title ',80a1)
c 
        write(iun,2270) (line(i),i=1,nchar+4)
        write(iunc,2270) (line(i),i=1,nchar+4)
c 
 2280 format('   show title')
        write(iun,2280)
        write(iunc,2280)
c 
c        find the limits for both x and y
c 
        xmin=1.0d30
        ymin=1.0d30
        xmax=-1.0d20
        ymax=-1.0d20
c 
        do 2300 i=1,n
        if(x(i) .lt. xmin) xmin=x(i)
        if(y(i) .lt. ymin) ymin=y(i)
        if(x(i) .gt. xmax) xmax=x(i)
        if(y(i) .gt. ymax) ymax=y(i)
 2300 continue
c 
        if(inumgr .eq. 1) goto 2340
c 
        do 2310 i=1,n2
        if(x2(i) .lt. xmin) xmin=x2(i)
        if(y2(i) .lt. ymin) ymin=y2(i)
        if(x2(i) .gt. xmax) xmax=x2(i)
        if(y2(i) .gt. ymax) ymax=y2(i)
 2310 continue
        if(inumgr .eq. 2) goto 2340
c 
        do 2320 i=1,n3
        if(x3(i) .lt. xmin) xmin=x3(i)
        if(y3(i) .lt. ymin) ymin=y3(i)
        if(x3(i) .gt. xmax) xmax=x3(i)
        if(y3(i) .gt. ymax) ymax=y3(i)
 2320 continue
c 
 2340 continue
c 
        xcenter=(xmin+xmax)/2
        ycenter=(ymax+ymin)/2
c 
        xsize=(xmax-xmin)
        ysize=(ymax-ymin)
        size=xsize
        if(ysize .gt. size) size=ysize
        size=size*1.2
c 
        xmin=xcenter-size/2
        xmax=xcenter+size/2
        ymin=ycenter-size/2
        ymax=ycenter+size/2
c 
c        set the size of the stupid thing
c 
cccc 2350 format(2x,' set size 0.75,1.0')
 2350 format(2x,' set size 0.76, 1.0')
        write(iunc,2350)
c 
c        write the instructions
c 
 2400 format(i6)
        iun2=iw+100000
        write (dummy,2400) iun2
        read(dummy,2000) anum8
c 
        do 2600 i=1,6
        file1(i+2)=anum1(i)
 2600 continue
c 
 2800 format('plot ', '[',  e8.2, ':', e8.2,'] ',
     1    '[',  e8.2, ':', e8.2,'] ',
     2    '"',1a8,'" ','notitle with dots')
  
 2803 format('plot ', '[',  e8.2, ':', e8.2,'] ',
     1    '[',  e8.2, ':', e8.2,'] ',
     2    '"',1a8,'" ','notitle with points')
  
 2805 format('plot ', '[',  e8.2, ':', e8.2,'] ',
     1    '[',  e8.2, ':', e8.2,'] ',
     2    '"',1a8,'" ','notitle with lines')
c 
        if( (inumgr .eq. 1) .and. (itype1 .eq. 1) )
     1      write(iun,2800) xmin,xmax,ymin,ymax,file8
c 
        if( (inumgr .eq. 1) .and. (itype1 .eq. 2) )
     1      write(iun,2803) xmin,xmax,ymin,ymax,file8
c 
        if( (inumgr .eq. 1) .and. (itype1 .eq. 3) )
     1      write(iun,2805) xmin,xmax,ymin,ymax,file8
c 
c 
        if( (inumgr .eq. 1) .and. (itype1 .eq. 1) )
     1      write(iunc,2800) xmin,xmax,ymin,ymax,file8
c 
        if( (inumgr .eq. 1) .and. (itype1 .eq. 2) )
     1      write(iunc,2803) xmin,xmax,ymin,ymax,file8
c 
        if( (inumgr .eq. 1) .and. (itype1 .eq. 3) )
     1      write(iunc,2805) xmin,xmax,ymin,ymax,file8
c 

 2830 format('plot ', '[',  e8.2, ':', e8.2,'] ',
     1    '[',  e8.2, ':', e8.2,'] ',
     2    '"',1a8,'" ','notitle with dots, ',1a1)
c 
 2833 format('plot ', '[',  e8.2, ':', e8.2,'] ',
     1    '[',  e8.2, ':', e8.2,'] ',
     2    '"',1a8,'" ','notitle with points, ',1a1)
c 
 2835 format('plot ', '[',  e8.2, ':', e8.2,'] ',
     1    '[',  e8.2, ':', e8.2,'] ',
     2    '"',1a8,'" ','notitle with lines, ',1a1)
c 
        if( (inumgr .ne. 1) .and. (itype1 .eq. 1) )
     1      write(iun,2830) xmin,xmax,ymin,ymax,file8,backslash
c 
        if( (inumgr .ne. 1) .and. (itype1 .eq. 2) )
     1      write(iun,2833) xmin,xmax,ymin,ymax,file8,backslash
c 
        if( (inumgr .ne. 1) .and. (itype1 .eq. 3) )
     1      write(iun,2835) xmin,xmax,ymin,ymax,file8,backslash

        if( (inumgr .ne. 1) .and. (itype1 .eq. 1) )
     1      write(iunc,2830) xmin,xmax,ymin,ymax,file8,backslash
c 
        if( (inumgr .ne. 1) .and. (itype1 .eq. 2) )
     1      write(iunc,2833) xmin,xmax,ymin,ymax,file8,backslash
c 
        if( (inumgr .ne. 1) .and. (itype1 .eq. 3) )
     1      write(iunc,2835) xmin,xmax,ymin,ymax,file8,backslash
c 
c        write the first data file to be plotted
c 
        iun22=88
        open(unit=iun22,file=file8)
c 
        write(iun22,3000) (x(i),y(i),i=1,n)
c 
        close(iun22)
c 
c       if the user so requested - write the instructions for the
c       plotting the second data file
c 
        if(inumgr .eq. 1) close(iun)
        if(inumgr .eq. 1) close(iunc)
        if(inumgr .eq. 1) return
c 
        iun2=iw+200000
        write (dummy,2400) iun2
        read(dummy,2000) anum8
c 
        do 2840 i=1,6
        file1(i+2)=anum1(i)
 2840 continue
c 
 2850 format('     ', '"',1a8,'" ','notitle with dots')
 2851 format('     ', '"',1a8,'" ','notitle with points')
 2852 format('     ', '"',1a8,'" ','notitle with lines')
c 
         if( (inumgr .eq. 2) .and. (itype2 .eq. 1) )
     1       write(iun,2850) file8
c 
         if( (inumgr .eq. 2) .and. (itype2 .eq. 2) )
     1       write(iun,2851) file8
c 
         if( (inumgr .eq. 2) .and. (itype2 .eq. 3) )
     1       write(iun,2852) file8
c 
c 
         if( (inumgr .eq. 2) .and. (itype2 .eq. 1) )
     1       write(iunc,2850) file8
c 
         if( (inumgr .eq. 2) .and. (itype2 .eq. 2) )
     1       write(iunc,2851) file8
c 
         if( (inumgr .eq. 2) .and. (itype2 .eq. 3) )
     1       write(iunc,2852) file8
c 
 2855 format('     ', '"',1a8,'" ','notitle with dots, ',1a1)
 2856 format('     ', '"',1a8,'" ','notitle with points, ',1a1)
 2857 format('     ', '"',1a8,'" ','notitle with lines, ',1a1)
c 
         if( (inumgr .ne. 2) .and. (itype2 .eq. 1) )
     1       write(iun,2855) file8,backslash
c 
         if( (inumgr .ne. 2) .and. (itype2 .eq. 2) )
     1       write(iun,2856) file8,backslash
c 
         if( (inumgr .ne. 2) .and. (itype2 .eq. 3) )
     1       write(iun,2857) file8,backslash
c
c 
         if( (inumgr .ne. 2) .and. (itype2 .eq. 1) )
     1       write(iunc,2855) file8,backslash
c 
         if( (inumgr .ne. 2) .and. (itype2 .eq. 2) )
     1       write(iunc,2856) file8,backslash
c 
         if( (inumgr .ne. 2) .and. (itype2 .eq. 3) )
     1       write(iunc,2857) file8,backslash
c 
c        write the second data file to be plotted
c 
        iun22=88
        open(unit=iun22,file=file8)
c 
        write(iun22,3000) (x2(i),y2(i),i=1,n2)
c 
        close(iun22)
c 
c       if the user so requested - write the instructions for the
c       plotting the third data file
c 
        if(inumgr .eq. 2) close(iun)
        if(inumgr .eq. 2) close(iunc)
        if(inumgr .eq. 2) return
c  
        iun2=iw+300000
        write (dummy,2400) iun2
        read(dummy,2000) anum8
c 
        do 2860 i=1,6
        file1(i+2)=anum1(i)
 2860 continue
c 
         if(itype3 .eq. 1) write(iun,2870) file8
         if(itype3 .eq. 2) write(iun,2871) file8
         if(itype3 .eq. 3) write(iun,2872) file8
c 
         if(itype3 .eq. 1) write(iunc,2870) file8
         if(itype3 .eq. 2) write(iunc,2871) file8
         if(itype3 .eq. 3) write(iunc,2872) file8

 2870 format('     ','"',1a8,'"','notitle  with dots')
 2871 format('     ','"',1a8,'"','notitle  with points')
 2872 format('     ','"',1a8,'"','notitle  with lines')
c 
c        write the third data file to be plotted
c 
        iun22=88
        open(unit=iun22,file=file8)
c 
 3000 format(2x,e11.5,2x,e11.5)
c 
        write(iun22,3000) (x3(i),y3(i),i=1,n3)
c 
        close(iun22)
        close(iun)
        close(iunc)
c 
        return
        end
c 
c 
c 
c 
c 
        SUBROUTINE quamesslen(MES,nchar,line)
        save
        CHARACTER *1 MES(1),AST,line(1)
        DATA AST/'*'/
C 
C         DETERMINE THE LENGTH OF THE MESSAGE
C 
        I=0
        DO 1400 I=1,10000
        IF(MES(I).EQ.AST) GOTO 1600
        I1=I
 1400 CONTINUE
 1600 CONTINUE
c 
        nchar=i1
        do 1800 i=1,nchar
        line(i)=mes(i)
 1800 continue
         RETURN
         END
  
  
  
  
c 
c 
c 
c 
c 
        subroutine zquaplot(iw,z,n,itype1,title)
        implicit real *8 (a-h,o-z)
        save
        dimension z(2,1),z2(2,1),z3(2,1)
        character *1 title(1)
c 
        inumgr=1
c 
        call zquaplo0(iw,z,z2,z3,n,n2,n3,
     1      title,inumgr,itype1,itype2,itype3)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine zquaplot2(iw,z,n,itype1,z2,n2,itype2,title)
        implicit real *8 (a-h,o-z)
        save
        dimension z(2,1),z2(2,1),z3(2,1)
        character *1 title(1)
c 
        inumgr=2
c 
        call zquaplo0(iw,z,z2,z3,n,n2,n3,
     1      title,inumgr,itype1,itype2,itype3)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine zquaplot3(iw,z,n,itype1,z2,n2,itype2,
     1      z3,n3,itype3,title)
        implicit real *8 (a-h,o-z)
        save
        dimension z(2,1),z2(2,1),z3(2,1)
        character *1 title(1)
c 
        inumgr=3
c 
        call zquaplo0(iw,z,z2,z3,n,n2,n3,
     1      title,inumgr,itype1,itype2,itype3)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine zquagraph(iw,z,n,itype1,title)
        implicit real *8 (a-h,o-z)
        save
        dimension z(2,1),z2(2,1),z3(2,1)
        character *1 title(1)
c 
        inumgr=1
c 
cccc        call zquagrap0(iw,x,y,x2,y2,x3,y3,n,n2,n3,
        call zquagrap0(iw,z,z2,z3,n,n2,n3,
     1      title,inumgr,itype1,itype2,itype3)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine zquagraph2(iw,z,n,itype1,z2,n2,itype2,title)
        implicit real *8 (a-h,o-z)
        save
        dimension z(2,1),z2(2,1),z3(2,1)
        character *1 title(1)
c 
        inumgr=2
c 
        call zquagrap0(iw,z,z2,z3,n,n2,n3,
     1      title,inumgr,itype1,itype2,itype3)
c 
        return
        end
c 
c 
c 
c 
c 
c 
        subroutine zquagraph3(iw,z,n,itype1,z2,n2,itype2,
     1      z3,n3,itype3,title)
        implicit real *8 (a-h,o-z)
        save
        dimension z(2,1),z2(2,1),z3(2,1)
        character *1 title(1)
c 
        inumgr=3
c 
        call zquagrap0(iw,z,z2,z3,n,n2,n3,
     1      title,inumgr,itype1,itype2,itype3)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine zquagrap0(iw,z,z2,z3,n,n2,n3,
     1      title,inumgr,itype1,itype2,itype3)
        implicit real *8 (a-h,o-z)
        save
        dimension z(2,1),z2(2,1),z3(2,1)
        character *1 title(1),gn(2),file1(8),anum1(8),line(82),
     1      blank,quo,backslash
        character *8 dummy,anum8,file8
c 
        equivalence (file1,file8),(anum1,anum8)
        data gn/'g','n'/,blank/' '/,quo/'"'/,backslash/'\\'/
c 
c        convert the user-specified Fortran unit number to
c        character format
c 
        if( (iw .ge. 0) .and. (iw .le. 9) )  write(dummy,1100) iw
 1100 format(i1)
c 
        if( (iw .ge. 10) .and. (iw .le. 99) ) write(dummy,1200) iw
 1200 format(i2)
c 
        if( (iw .ge. 100) .and. (iw .le. 999) ) write(dummy,1300) iw
 1300 format(i3)
c 
        if( (iw .ge. 1000) .and. (iw .le. 9999) ) write(dummy,1400) iw
 1400 format(i4)
c 
 2000 format(1a8)
        read(dummy,2000) anum8
c 
c        construct the file name on which the Gnuplot instructions
c        are to be written
c 
        file1(1)=gn(1)
        file1(2)=gn(2)
        do 2200 i=1,6
        file1(i+2)=anum1(i)
 2200 continue
c 
c        open the fortran file with the unit 87 and name file8
c 
        iun=87
        open(unit=iun,file=file8)
c 
 2250 format('   # set terminal postscript',/,
     1    '   # set output "plot.ps"')
c 
        write(iun,2250)
c 
 2255 format('    set terminal postscript',/,
     1    '    set output "plot.ps"')
c 
c        generate the title for the plot
c 
        line(1)=blank
        line(2)=blank
        line(3)=quo
c 
        call quamesslen(title,nchar,line(4))
c 
        line(nchar+4)=quo
 2300 format('  set title ',80a1)
        write(iun,2300) (line(i),i=1,nchar+4)
c 
 2350 format('   show title')
        write(iun,2350)
c 
c        write the instructions
c 
 2400 format(i6)
        iun2=iw+100000
        write (dummy,2400) iun2
        read(dummy,2000) anum8
c 
        do 2600 i=1,6
        file1(i+2)=anum1(i)
 2600 continue
c 
 2800 format('plot ','"',1a8,'"     ','notitle  with dots')
 2803 format('plot ','"',1a8,'"     ','notitle  with points')
 2805 format('plot ','"',1a8,'"     ','notitle  with lines')
c 
        if( (inumgr .eq. 1) .and. (itype1 .eq. 1) )
     1      write(iun,2800) file8
c 
        if( (inumgr .eq. 1) .and. (itype1 .eq. 2) )
     1      write(iun,2803) file8
c 
        if( (inumgr .eq. 1) .and. (itype1 .eq. 3) )
     1      write(iun,2805) file8
  
 2830 format('plot ','"',1a8,'"     ','notitle  with dots, ',1a1)
 2831 format('plot ','"',1a8,'"     ','notitle  with points, ',1a1)
 2832 format('plot ','"',1a8,'"     ','notitle  with lines, ',1a1)
c 
        if( (inumgr .ne. 1) .and. (itype1 .eq. 1) )
     1      write(iun,2830) file8,backslash
c 
        if( (inumgr .ne. 1) .and. (itype1 .eq. 2) )
     1       write(iun,2831) file8,backslash
c 
        if( (inumgr .ne. 1) .and. (itype1 .eq. 3) )
     1      write(iun,2832) file8,backslash
c 
c        write the first data file to be plotted
c 
        iun22=88
        open(unit=iun22,file=file8)
c 
cccc        write(iun22,3000) (x(i),y(i),i=1,n)
        write(iun22,3000) (z(1,i),z(2,i),i=1,n)
c 
        close(iun22)
c 
c       if the user so requested - write the instructions for the
c       plotting the second data file
c 
        if(inumgr .eq. 1) close(iun)
        if(inumgr .eq. 1) close(iunc)
        if(inumgr .eq. 1) return
  
        iun2=iw+200000
        write (dummy,2400) iun2
        read(dummy,2000) anum8
c 
        do 2840 i=1,6
        file1(i+2)=anum1(i)
 2840 continue
c 
 2850 format('     ','"',1a8,'"     ','notitle  with dots')
 2851 format('     ','"',1a8,'"     ','notitle  with points')
 2852 format('     ','"',1a8,'"     ','notitle  with lines')
c 
        if( (inumgr .eq. 2) .and. (itype2 .eq. 1) )
     1      write(iun,2850) file8
c 
        if( (inumgr .eq. 2) .and. (itype2 .eq. 2) )
     1      write(iun,2851) file8
c 
        if( (inumgr .eq. 2) .and. (itype2 .eq. 3) )
     1      write(iun,2852) file8
c 
c 
        if( (inumgr .eq. 3) .and. (itype2 .eq. 1) )
     1      write(iun,2855) file8,backslash
c 
        if( (inumgr .eq. 3) .and. (itype2 .eq. 2) )
     1      write(iun,2856) file8,backslash
c 
        if( (inumgr .eq. 3) .and. (itype2 .eq. 3) )
     1      write(iun,2857) file8,backslash
c 
 2855 format('     ','"',1a8,'"     ','notitle  with dots, ',1a1)
 2856 format('     ','"',1a8,'"     ','notitle  with points, ',1a1)
 2857 format('     ','"',1a8,'"     ','notitle  with lines, ',1a1)
c 
c        write the second data file to be plotted
c 
        iun22=88
        open(unit=iun22,file=file8)
c 
cccc        write(iun22,3000) (x2(i),y2(i),i=1,n2)
        write(iun22,3000) (z2(1,i),z2(2,i),i=1,n2)
c 
        close(iun22)
c 
c       if the user so requested - write the instructions for the
c       plotting the third data file
c 
        if(inumgr .eq. 2) close(iun)
        if(inumgr .eq. 2) close(iunc)
        if(inumgr .eq. 2) return
  
        iun2=iw+300000
        write (dummy,2400) iun2
        read(dummy,2000) anum8
c 
        do 2860 i=1,6
        file1(i+2)=anum1(i)
 2860 continue
c 
        if(itype3 .eq. 1) write(iun,2870) file8
        if(itype3 .eq. 2) write(iun,2871) file8
        if(itype3 .eq. 3) write(iun,2872) file8
c 
 2870 format('     ','"',1a8,'"     ','notitle  with dots')
 2871 format('     ','"',1a8,'"     ','notitle  with points')
 2872 format('     ','"',1a8,'"     ','notitle  with lines')
c 
c        write the third data file to be plotted
c 
        iun22=88
        open(unit=iun22,file=file8)
c 
 3000 format(2x,e11.5,2x,e11.5)
c 
cccc        write(iun22,3000) (x3(i),y3(i),i=1,n3)
        write(iun22,3000) (z3(1,i),z3(2,i),i=1,n3)
c 
        close(iun22)
        close(iun)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine zquaplo0(iw,z,z2,z3,n,n2,n3,
     1      title,inumgr,itype1,itype2,itype3)
        implicit real *8 (a-h,o-z)
        save
        dimension z(2,1),z2(2,1),z3(2,1)
        character *1 title(1),gn(2),file1(8),anum1(8),line(82),
     1      blank,quo,backslash
        character *8 dummy,anum8,file8
c 
        equivalence (file1,file8),(anum1,anum8)
        data gn/'g','n'/,blank/' '/,quo/'"'/,backslash/'\\'/
c 
c        convert the user-specified Fortran unit number to
c        character format
c 
        if( (iw .ge. 0) .and. (iw .le. 9) )  write(dummy,1100) iw
 1100 format(i1)
c 
        if( (iw .ge. 10) .and. (iw .le. 99) ) write(dummy,1200) iw
 1200 format(i2)
c 
        if( (iw .ge. 100) .and. (iw .le. 999) ) write(dummy,1300) iw
 1300 format(i3)
c 
        if( (iw .ge. 1000) .and. (iw .le. 9999) ) write(dummy,1400) iw
 1400 format(i4)
c 
 2000 format(1a8)
        read(dummy,2000) anum8
c 
c        construct the file name on which the Gnuplot instructions
c        are to be written
c 
        file1(1)=gn(1)
        file1(2)=gn(2)
        do 2200 i=1,6
        file1(i+2)=anum1(i)
 2200 continue
c 
c        open the fortran file with the unit 87 and name file8
c 
        iun=87
        open(unit=iun,file=file8)
c 
 2250 format('   # set terminal postscript',/,
     1    '   # set output "plot.ps"')
c 
        write(iun,2250)
c 
c        generate the title for the plot
c 
        line(1)=blank
        line(2)=blank
        line(3)=quo
c 
        call quamesslen(title,nchar,line(4))
c 
        line(nchar+4)=quo
 2270 format('  set title ',80a1)
c 
        write(iun,2270) (line(i),i=1,nchar+4)
c 
 2280 format('   show title')
        write(iun,2280)
c 
c        find the limits for both x and y
c 
        xmin=1.0d30
        ymin=1.0d30
        xmax=-1.0d20
        ymax=-1.0d20
c 
        do 2300 i=1,n
cccc        if(x(i) .lt. xmin) xmin=x(i)
cccc        if(y(i) .lt. ymin) ymin=y(i)
cccc        if(x(i) .gt. xmax) xmax=x(i)
cccc        if(y(i) .gt. ymax) ymax=y(i)
        if(z(1,i) .lt. xmin) xmin=z(1,i)
        if(z(2,i) .lt. ymin) ymin=z(2,i)
        if(z(1,i) .gt. xmax) xmax=z(1,i)
        if(z(2,i) .gt. ymax) ymax=z(2,i)
 2300 continue
c 
        if(inumgr .eq. 1) goto 2340
c 
        do 2310 i=1,n2
cccc        if(x2(i) .lt. xmin) xmin=x2(i)
cccc        if(y2(i) .lt. ymin) ymin=y2(i)
cccc        if(x2(i) .gt. xmax) xmax=x2(i)
cccc        if(y2(i) .gt. ymax) ymax=y2(i)
  
        if(z2(1,i) .lt. xmin) xmin=z2(1,i)
        if(z2(2,i) .lt. ymin) ymin=z2(2,i)
        if(z2(1,i) .gt. xmax) xmax=z2(1,i)
        if(z2(2,i) .gt. ymax) ymax=z2(2,i)
  
 2310 continue
  
  
cccc        call prin2('z2=*',z2,n2*2)
  
  
  
        if(inumgr .eq. 2) goto 2340
c 
        do 2320 i=1,n3
cccc        if(x3(i) .lt. xmin) xmin=x3(i)
cccc        if(y3(i) .lt. ymin) ymin=y3(i)
cccc        if(x3(i) .gt. xmax) xmax=x3(i)
cccc        if(y3(i) .gt. ymax) ymax=y3(i)
  
  
        if(z3(1,i) .lt. xmin) xmin=z3(1,i)
        if(z3(2,i) .lt. ymin) ymin=z3(2,i)
        if(z3(1,i) .gt. xmax) xmax=z3(1,i)
        if(z3(2,i) .gt. ymax) ymax=z3(2,i)
  
  
 2320 continue
c 
 2340 continue
c 
        xcenter=(xmin+xmax)/2
        ycenter=(ymax+ymin)/2
c 
        xsize=(xmax-xmin)
        ysize=(ymax-ymin)
        size=xsize
        if(ysize .gt. size) size=ysize
        size=size*1.2
c 
        xmin=xcenter-size/2
        xmax=xcenter+size/2
        ymin=ycenter-size/2
        ymax=ycenter+size/2
c 
c        set the size of the stupid thing
c 
cccc 2350 format(2x,' set size 0.75,1.0')
 2350 format(2x,' set size 1.0, 1.0')
        write(iun,2350)
c 
c        write the instructions
c 
 2400 format(i6)
        iun2=iw+100000
        write (dummy,2400) iun2
        read(dummy,2000) anum8
c 
        do 2600 i=1,6
        file1(i+2)=anum1(i)
 2600 continue
c 
 2800 format('plot ', '[',  e8.2, ':', e8.2,'] ',
     1    '[',  e8.2, ':', e8.2,'] ',
     2    '"',1a8,'" ','notitle with dots')
  
 2803 format('plot ', '[',  e8.2, ':', e8.2,'] ',
     1    '[',  e8.2, ':', e8.2,'] ',
     2    '"',1a8,'" ','notitle with points')
  
 2805 format('plot ', '[',  e8.2, ':', e8.2,'] ',
     1    '[',  e8.2, ':', e8.2,'] ',
     2    '"',1a8,'" ','notitle with lines')
c 
        if( (inumgr .eq. 1) .and. (itype1 .eq. 1) )
     1      write(iun,2800) xmin,xmax,ymin,ymax,file8
c 
        if( (inumgr .eq. 1) .and. (itype1 .eq. 2) )
     1      write(iun,2803) xmin,xmax,ymin,ymax,file8
c 
        if( (inumgr .eq. 1) .and. (itype1 .eq. 3) )
     1      write(iun,2805) xmin,xmax,ymin,ymax,file8
c 
 2830 format('plot ', '[',  e8.2, ':', e8.2,'] ',
     1    '[',  e8.2, ':', e8.2,'] ',
     2    '"',1a8,'" ','notitle with dots, ',1a1)
c 
 2833 format('plot ', '[',  e8.2, ':', e8.2,'] ',
     1    '[',  e8.2, ':', e8.2,'] ',
     2    '"',1a8,'" ','notitle with points, ',1a1)
c 
 2835 format('plot ', '[',  e8.2, ':', e8.2,'] ',
     1    '[',  e8.2, ':', e8.2,'] ',
     2    '"',1a8,'" ','notitle with lines, ',1a1)
c 
        if( (inumgr .ne. 1) .and. (itype1 .eq. 1) )
     1      write(iun,2830) xmin,xmax,ymin,ymax,file8,backslash
c 
        if( (inumgr .ne. 1) .and. (itype1 .eq. 2) )
     1      write(iun,2833) xmin,xmax,ymin,ymax,file8,backslash
c 
        if( (inumgr .ne. 1) .and. (itype1 .eq. 3) )
     1      write(iun,2835) xmin,xmax,ymin,ymax,file8,backslash
c 
c        write the first data file to be plotted
c 
        iun22=88
        open(unit=iun22,file=file8)
c 
cccc        write(iun22,3000) (x(i),y(i),i=1,n)
        write(iun22,3000) (z(1,i),z(2,i),i=1,n)
c 
        close(iun22)
c 
c       if the user so requested - write the instructions for the
c       plotting the second data file
c 
        if(inumgr .eq. 1) close(iun)
        if(inumgr .eq. 1) return
c 
        iun2=iw+200000
        write (dummy,2400) iun2
        read(dummy,2000) anum8
c 
        do 2840 i=1,6
        file1(i+2)=anum1(i)
 2840 continue
c 
 2850 format('     ', '"',1a8,'" ','notitle with dots')
 2851 format('     ', '"',1a8,'" ','notitle with points')
 2852 format('     ', '"',1a8,'" ','notitle with lines')
c 
         if( (inumgr .eq. 2) .and. (itype2 .eq. 1) )
     1       write(iun,2850) file8
c 
         if( (inumgr .eq. 2) .and. (itype2 .eq. 2) )
     1       write(iun,2851) file8
c 
         if( (inumgr .eq. 2) .and. (itype2 .eq. 3) )
     1       write(iun,2852) file8
c 
 2855 format('     ', '"',1a8,'" ','notitle with dots, ',1a1)
 2856 format('     ', '"',1a8,'" ','notitle with points, ',1a1)
 2857 format('     ', '"',1a8,'" ','notitle with lines, ',1a1)
c 
         if( (inumgr .ne. 2) .and. (itype2 .eq. 1) )
     1       write(iun,2855) file8,backslash
c 
         if( (inumgr .ne. 2) .and. (itype2 .eq. 2) )
     1       write(iun,2856) file8,backslash
c 
         if( (inumgr .ne. 2) .and. (itype2 .eq. 3) )
     1       write(iun,2857) file8,backslash
c 
c        write the second data file to be plotted
c 
        iun22=88
        open(unit=iun22,file=file8)
c 
cccc        write(iun22,3000) (x2(i),y2(i),i=1,n2)
        write(iun22,3000) (z2(1,i),z2(2,i),i=1,n2)
c 
        close(iun22)
c 
c       if the user so requested - write the instructions for the
c       plotting the third data file
c 
        if(inumgr .eq. 2) close(iun)
        if(inumgr .eq. 2) return
  
        iun2=iw+300000
        write (dummy,2400) iun2
        read(dummy,2000) anum8
c 
        do 2860 i=1,6
        file1(i+2)=anum1(i)
 2860 continue
c 
         if(itype3 .eq. 1) write(iun,2870) file8
         if(itype3 .eq. 2) write(iun,2871) file8
         if(itype3 .eq. 3) write(iun,2872) file8
 2870 format('     ','"',1a8,'"','notitle  with dots')
 2871 format('     ','"',1a8,'"','notitle  with points')
 2872 format('     ','"',1a8,'"','notitle  with lines')
c 
c        write the third data file to be plotted
c 
        iun22=88
        open(unit=iun22,file=file8)
c 
 3000 format(2x,e11.5,2x,e11.5)
c 
cccc        write(iun22,3000) (x3(i),y3(i),i=1,n3)
        write(iun22,3000) (z3(1,i),z3(2,i),i=1,n3)
c 
        close(iun22)
        close(iun)
c 
        return
        end
  
c
c
c
c
c
        subroutine exograph2d(iw,ts,n,vals,isub,title)
        implicit real *8 (a-h,o-z)
        save
        dimension ts(n),vals(n,n)
        character *8 file8,anum8,dummy
        character *1 title(1),file1(8),anum1(8),line(82),
     1      blank,quo
        equivalence (file1,file8),(anum1,anum8)
        data blank/' '/,quo/'"'/
c
c
c        convert the user-specified Fortran unit number to 
c        character format
c
 2000 format(1a8)
c
        call qua_filename(iw,file8)
c
c        open the fortran file with the unit 87 and name file8
c 
        iun=87
        open(unit=iun,file=file8)
c
 2250 format('   # set terminal postscript',/,
     1    '   # set output "plot.ps"')
c
        write(iun,2250)
c
c        generate the title for the plot
c
        line(1)=blank
        line(2)=blank
        line(3)=quo
c
        call quamesslen(title,nchar,line(4))
c
        line(nchar+4)=quo
 2270 format('  set title ',80a1)
c
        write(iun,2270) (line(i),i=1,nchar+4)
c
 2280 format('   show title')
        write(iun,2280)
c
c
c        write the instructions 
c
 2400 format(i6)
        iun2=iw+100000
        write (dummy,2400) iun2
        read(dummy,2000) anum8
c
        do 2600 i=1,6
        file1(i+2)=anum1(i)
 2600 continue
c
 2800 format('splot ','"',1a8,'"     ','notitle  with lines')

        write(iun,2800) file8
c
c       create trhe data file
c
c
c        write the data file to be plotted
c
        iun22=88
        iss=isub
        if(iss .eq. 0) iss=1
c
        open(unit=iun22,file=file8)
c
        do 3400 i=1,n,iss
        do 3200 j=1,n,iss
c
        ii=ii+1
        write(iun22,3000) ts(i),ts(j),vals(i,j)
c
 3000 format(2x,e11.5,2x,e11.5,2x,e11.5)
 3100 format(20x)
c
 3200 continue
c
        write(iun22,3100)
 3400 continue
c
        close(iun22)
c

        return
        end
c
c
c
c
c
        subroutine qua_filename(iw,file88)
        implicit real *8 (a-h,o-z)
        character *1 gn(2),file1(8),anum1(8)
        character *8 dummy,anum8,file8,file88
c
        equivalence (file1,file8),(anum1,anum8)
        data gn/'g','n'/
c
c        convert the user-specified Fortran unit number to 
c        character format
c
        if( (iw .ge. 0) .and. (iw .le. 9) )  write(dummy,1100) iw
 1100 format(i1)
c
        if( (iw .ge. 10) .and. (iw .le. 99) ) write(dummy,1200) iw
 1200 format(i2)
c
        if( (iw .ge. 100) .and. (iw .le. 999) ) write(dummy,1300) iw
 1300 format(i3)
c
        if( (iw .ge. 1000) .and. (iw .le. 9999) ) write(dummy,1400) iw
 1400 format(i4)
c
 2000 format(1a8)
        read(dummy,2000) anum8
c
c        construct the file name on which the Gnuplot instructions
c        are to be written
c
        file1(1)=gn(1)
        file1(2)=gn(2)
        do 2200 i=1,6
        file1(i+2)=anum1(i)
 2200 continue
c
        file88=file8
c        
        return
        end
