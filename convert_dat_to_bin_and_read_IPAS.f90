      program main

        integer :: i,j,k
        integer, parameter :: ii = 50, jj = 28, kk = 3
        real :: phi(ii),r(jj)
        real :: an(ii,jj),cn(ii,jj),a(ii,jj),d(ii,jj)
        real var1(ii,jj+1),var2(ii,jj+1),var3(ii,jj+1),var4(ii,jj+1)

        r(1) = 1
        do j=2,jj
           if(r(j-1).lt.10)then
              r(j) = r(j-1) + 1
           else if(r(j-1).ge.10.and.r(j-1).lt.100)then
              r(j) = r(j-1) + 10
           else if(r(j-1).ge.100)then
              r(j) = r(j-1) + 100
           end if
        end do

        open(9,file='a_avg_lookup_new.txt')
        open(10,file='a_n_lookup_dd.txt')
        open(11,file='c_n_lookup_dd.txt')
        open(12,file='dd_lookup.txt')
        do i=1,ii
           read(9,*)  (var1(i,j),j=1,jj+1)
           read(10,*) (var2(i,j),j=1,jj+1)
           read(11,*) (var3(i,j),j=1,jj+1)
           read(12,*) (var4(i,j),j=1,jj+1)
        end do
        close(9)
        close(10)
        close(11)
        close(12)

        do i=1,ii
           phi(i) = var1(i,1)
           do j=1,jj
              a(i,j) = var1(i,j+1)  
              an(i,j) = var2(i,j+1) 
              cn(i,j) = var3(i,j+1)
              d(i,j) = var4(i,j+1)
           enddo
        enddo
        
        OPEN(8,file='ACOLLdd.bin',form='unformatted')
        WRITE(8) (phi(i),i=1,ii)
        WRITE(8) (r(j),j=1,jj)
        WRITE(8) ((a(i,j),i=1,ii),j=1,jj)
        WRITE(8) ((an(i,j),i=1,ii),j=1,jj)
        WRITE(8) ((cn(i,j),i=1,ii),j=1,jj)
        WRITE(8) ((d(i,j),i=1,ii),j=1,jj)
        CLOSE(8)

        phi=0
        r=0
        a=0
        an=0
        cn=0
        d=0

        OPEN(8,file='ACOLLdd.bin',form='unformatted')
        READ(8) (phi(i),i=1,ii)
        READ(8) (r(j),j=1,jj)
        READ(8) ((a(i,j),i=1,ii),j=1,jj)
        READ(8) ((an(i,j),i=1,ii),j=1,jj)
        READ(8) ((cn(i,j),i=1,ii),j=1,jj)
        READ(8) ((d(i,j),i=1,ii),j=1,jj)
        CLOSE(8)

        
        do i=1,ii
           do j=1,jj   
              print*,i,j,phi(i),r(j),a(i,j),an(i,j),cn(i,j),d(i,j)
           enddo
        enddo
        
      end program main

