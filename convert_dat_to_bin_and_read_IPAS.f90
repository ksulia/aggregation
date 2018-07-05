      program main

        integer :: i,j,k
        integer, parameter :: ii = 50, jj = 28, kk = 3
        real :: phi(ii),r(jj)
        real :: an(ii,jj),cn(ii,jj),a(ii,jj)
        real var(ii,jj+1)

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
        do i=1,ii
           read(9,*) (var(i,j),j=1,jj+1)
        end do
        close(9)

        do i=1,ii
           phi(i) = var(i,1)
           do j=1,jj
              a(i,j) = var(i,j+1)    
           enddo
        enddo
        
        OPEN(8,file='ACOLL.bin',form='unformatted')
        WRITE(8) (phi(i),i=1,ii)
        WRITE(8) (r(j),j=1,jj)
        WRITE(8) ((a(i,j),i=1,ii),j=1,jj)
        CLOSE(8)

        phi=0
        r=0
        a=0

        do i=1,ii
           do j=1,jj   
              print*,i,j,phi(i),r(j),a(i,j)
           enddo
        enddo


        OPEN(8,file='ACOLL.bin',form='unformatted')
        READ(8) (phi(i),i=1,ii)
        READ(8) (r(j),j=1,jj)
        READ(8) ((a(i,j),i=1,ii),j=1,jj)
        CLOSE(8)

        
        do i=1,ii
           do j=1,jj   
              print*,i,j,phi(i),r(j),a(i,j)
           enddo
        enddo
        
      end program main

