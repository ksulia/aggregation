PROGRAM COLLECTION_LOOKUP_GENERATOR_LARGE

  IMPLICIT NONE


     integer :: i, j, k, l, m, n, o, p, q, c
     integer, parameter :: E = 1.0, nsum = 80
     real, parameter :: drr = 10.e-6
     real r(nsum), dr(nsum)
     integer, parameter :: ii=5, jj=6, kk=jj, ll=8, mm=9, nn=nsum, oo=nsum, pp=nsum, qq=nsum
     real gammln, pi, gam, maxx, maxy, area, vol1, vol2, n1, n2, v1, v2
     real ncheck(nn,oo), mrate(qq,pp,oo,nn), nrate(qq,pp,oo,nn)
     real coll, ncoll, ni(ii), an(jj), cn(kk), nu(ll), rho(mm), a1(nn), c1(oo), a2(pp), c2(qq) 
     real n3(pp), v(oo,nn)
     real start, finish

     integer ncid,vid,nid,aid,cid,nuid,rhoid
     character(7) :: fname

     pi = 3.14159265359

     r(1) = 0.5
     r(2) = 1.0
     r(3) = 1.5
     dr(1) = 0.5e-6
     dr(2) = (r(2) - r(1))*1.e-6
     dr(3) = (r(3) - r(2))*1.e-6
     c = 0
     DO n = 4, nn, 1
        c = c + 1
        r(n) = r(n-1) + 5.*10.**floor((c-18.)/18.)
        dr(n) = (r(n) - r(n-1))*1.e-6
        !print*,r(n),dr(n)
     END DO
     !r = [0.5,1.0,1.5-->10,15,20,25-->100,150,200,250-->1000,---->100000] microns
     DO n = 1, nn, 1
         DO o = 1, oo, 1
             maxx = max(r(n)*1.e6,r(o)*1.e6)
             !!Heymsfield, Bansemer, and Schmitt 2013
             if(2.*maxx .lt. 41) v(o,n) = 0.0028*(2.*maxx)**2.
             if(2.*maxx .ge. 41 .and. 2.*maxx .lt. 839) v(o,n) = 0.0791*(2.*maxx)**1.101
             if(2.*maxx .ge. 839) v(o,n) = 62.29*(2.*maxx)**0.1098
             v(o,n) = v(o,n)/100.
         END DO
     END DO

     OPEN(1,FILE="COLLV4.dat")
     DO i = 1,ii,1
        ni(i) = 10.**(i-1)   *1000. !L-1 --> m-3 [1,10,100,1000,10000 /L]

        DO j = 1,jj,1
           an(j) = 10.**(j-1) *1.e-6 !microns --> m 

           DO k = 1,kk,1
              cn(k) = 10.**(k-1)  *1.e-6

              DO l = 1,ll,1
                 nu(l) = real(l)
                 gam = exp(gammln(nu(l)))
                 !$OMP PARALLEL DO
                 DO m = 1,mm,1
                    rho(m) = real(m)*100.

                    !print*,ni(i),an(j),cn(k),nu(l),rho(m)
                    coll = 0
                    ncoll = 0

                    !!These inner loops should have a very fine resolution to make sure we are 
                    !!integrating over "all possible" sizes
                    call cpu_time(start)
                    DO n = 1,nn,1                  
                       a1(n) = r(n)  *1.e-6

                       DO o = 1,oo,1
                         c1(o) = r(o)  *1.e-6

                         maxx = max(a1(n),c1(o))
                         v1 = v(o,n)
                         vol1 = a1(n)*a1(n)*c1(o)
                         n1 = ni(i)/gam*gam * (a1(n)/an(j))**(nu(l)-1.)*(exp(-a1(n)/an(j))/an(j)) *&
                                                   (c1(o)/cn(k))**(nu(l)-1.)*(exp(-c1(o)/cn(k))/cn(k))

                         DO p = 1,pp,1
                            a2(p) = r(p)  *1.e-6

                            DO q = 1,qq,1
                               c2(q) = r(q)  *1.e-6
        
                               maxy = max(a2(p),c2(q))    
                               area = (maxx + maxy)*(maxx + maxy)
                               
                               v2 = v(q,p)
                               vol2 = a2(p)*a2(p)*c2(q)

                               n2 = ni(i)/gam*gam * (a2(p)/an(j))**(nu(l)-1.)*(exp(-a2(p)/an(j))/an(j)) *&
                                                   (c2(q)/cn(k))**(nu(l)-1.)*(exp(-c2(q)/cn(k))/cn(k))


                               mrate(q,p,o,n) = pi*rho(m)*(vol1+vol2)*area*abs(v1-v2)*n1*n2*E*dr(n)*dr(o)*dr(q)*dr(p)
                               nrate(q,p,o,n) = pi*area*abs(v1-v2)*n1*n2*E*dr(n)*dr(o)*dr(q)*dr(p)
                             END DO!c2(q)

                          END DO!a2(p)

                       END DO!c1(o)

                    END DO!a1(n)

                    coll = sum(sum(sum(sum(mrate,dim=4),dim=3),dim=2),dim=1)  !kg/m3/s
                    ncoll = sum(sum(sum(sum(nrate,dim=4),dim=3),dim=2),dim=1) !#/m3/s
                    call cpu_time(finish)
                    WRITE(*,*) finish-start," s  ", ni(i), an(j), cn(k), nu(l), rho(m), coll,ncoll
                    WRITE(1,*) ni(i), an(j), cn(k), nu(l), rho(m), coll, ncoll
                 END DO!rho(m)
                !$OMP END PARALLEL DO

              END DO!nu(l)

           END DO!cn(k)

        END DO!an(j)

     END DO!ni(i)
     CLOSE(1)
  
!     coll = coll*pi  !!divide by air density in model!!


END PROGRAM COLLECTION_LOOKUP_GENERATOR_LARGE



    !_____________________________________________________________________
REAL FUNCTION GAMMLN(XX)
      !---------------------------------------------------------------------
      INTEGER J
      REAL*8 COF(6),STP,HALF,ONE,FPF,X,TMP,SER
      REAL XX
      DATA COF,STP/76.18009173D0,-86.50532033D0,24.01409822D0,&
           -1.231739516D0,.120858003D-2,-.536382D-5,2.50662827465D0/
      DATA HALF,ONE,FPF/0.5D0,1.0D0,5.5D0/
      X=XX-ONE
      TMP=X+FPF
      TMP=(X+HALF)*LOG(TMP)-TMP
      SER=ONE
      DO J=1,6
         X=X+ONE
         SER=SER+COF(J)/X
      ENDDO
      GAMMLN=TMP+LOG(STP*SER)
      RETURN
END FUNCTION GAMMLN

