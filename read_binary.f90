      program main

        integer :: i,j,k,l,m
        integer, parameter :: ii = 5, jj = 4, kk = 4, ll = 8, mm = 9
        !integer, parameter :: ii = 1, jj = 1, kk = 1, ll = 1, mm = 9
        real :: ni(ii),an(jj),cn(kk),nu(ll),rho(mm)
        real :: coll(ii,jj,kk,ll,mm),ncoll(ii,jj,kk,ll,mm)

        open(8,file='COLLV2.bin',form='unformatted')
        
        READ(8) (ni(i),i=1,ii)
        READ(8) (an(j),j=1,jj)
        READ(8) (cn(k),k=1,kk)
        READ(8) (nu(l),l=1,ll)
        READ(8) (rho(m),m=1,mm)
        READ(8) (((((coll(i,j,k,l,m),i=1,ii),j=1,jj),k=1,kk),l=1,ll),m=1,mm)
        READ(8) (((((ncoll(i,j,k,l,m),i=1,ii),j=1,jj),k=1,kk),l=1,ll),m=1,mm)

        do i=1,ii
           do j=1,jj
              do k=1,kk
                 do l=1,ll
                    do m=1,mm
                       write(*,*) ni(i),an(j),cn(k),nu(l),rho(m),coll(i,j,k,l,m),ncoll(i,j,k,l,m)
                    enddo
                 enddo
              enddo
           enddo
        enddo
        
        
        
        
      end program main

