      program main

        integer :: i,j,k,l,m
        integer, parameter :: ii = 5, jj = 6, kk = 6, ll = 8, mm = 9
        real :: ni(ii),an(jj),cn(kk),nu(ll),rho(mm)
        real :: coll(ii,jj,kk,ll,mm),ncoll(ii,jj,kk,ll,mm)
        real :: n(ii),a(jj),c(kk),nn(ll),r(mm)
        real :: co(ii,jj,kk,ll,mm),nco(ii,jj,kk,ll,mm)

        !!!if COLLV2.dat is used and dates back to Oct 4 2018,
        !!must multiply ni by 10^6 coll and ncoll by 10^12 to correct old error
        
        open(9,file='COLLV3.dat') 
        do i=1,ii
           do j=1,jj
              do k=1,kk
                 do l=1,ll
                    do m=1,mm
                       read(9,*) ni(i),an(j),cn(k),nu(l),rho(m),coll(i,j,k,l,m),ncoll(i,j,k,l,m)
                    enddo
                 enddo
              enddo
           enddo
        enddo
        close(9)
        
        open(8,file='COLLV3.bin',form='unformatted')

        WRITE(8) (ni(i),i=1,ii)
        WRITE(8) (an(j),j=1,jj)
        WRITE(8) (cn(k),k=1,kk)
        WRITE(8) (nu(l),l=1,ll)
        WRITE(8) (rho(m),m=1,mm)
        WRITE(8) (((((coll(i,j,k,l,m),i=1,ii),j=1,jj),k=1,kk),l=1,ll),m=1,mm)
        WRITE(8) (((((ncoll(i,j,k,l,m),i=1,ii),j=1,jj),k=1,kk),l=1,ll),m=1,mm)

        !if COLLV2.dat is dated Oct 4, 2018, uncomment below and comment above.

        !WRITE(8) (ni(i)*10**6,i=1,ii)
        !WRITE(8) (an(j),j=1,jj)
        !WRITE(8) (cn(k),k=1,kk)
        !WRITE(8) (nu(l),l=1,ll)
        !WRITE(8) (rho(m),m=1,mm)
        !WRITE(8) (((((coll(i,j,k,l,m)*10**12,i=1,ii),j=1,jj),k=1,kk),l=1,ll),m=1,mm)
        !WRITE(8) (((((ncoll(i,j,k,l,m)*10**12,i=1,ii),j=1,jj),k=1,kk),l=1,ll),m=1,mm)

        CLOSE(8)


        open(8,file='COLLV3.bin',form='unformatted')
        READ(8) (n(i),i=1,ii)
        READ(8) (a(j),j=1,jj)
        READ(8) (c(k),k=1,kk)
        READ(8) (nn(l),l=1,ll)
        READ(8) (r(m),m=1,mm)
        READ(8) (((((co(i,j,k,l,m),i=1,ii),j=1,jj),k=1,kk),l=1,ll),m=1,mm)
        READ(8) (((((nco(i,j,k,l,m),i=1,ii),j=1,jj),k=1,kk),l=1,ll),m=1,mm)
        CLOSE(8)

        
        do i=1,ii
           do j=1,jj
              do k=1,kk
                 do l=1,ll
                    do m=1,mm
                       write(*,*) n(i),a(j),c(k),nn(l),r(m),co(i,j,k,l,m),nco(i,j,k,l,m)
                    enddo
                 enddo
              enddo
           enddo
        enddo
        
        
        
        
      end program main

