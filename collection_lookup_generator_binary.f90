PROGRAM COLLECTION_LOOKUP_GENERATOR_BINARY

  IMPLICIT NONE


     integer :: i, j, k, l, m, n, o, p, q
     integer, parameter :: E = 1.0, nsum = 100
     !real, parameter :: dr = 10.e-6
     integer, parameter :: ii=5, jj=4, kk=jj, ll=8, mm=9, nn=nsum, oo=nsum, pp=nsum, qq=nsum
     !integer, parameter :: ii=1, jj=1, kk=jj, ll=1, mm=9, nn=nsum, oo=nsum, pp=nsum, qq=nsum
     real dr(nsum)
     real gammln, pi, gam, maxx, maxy, area, vol1, vol2, n1, n2, v1, v2
     real ncheck(nn,oo), mrate(nn,oo,pp,qq), nrate(nn,oo,pp,qq)
     real coll(ii,jj,kk,ll,mm), ncoll(ii,jj,kk,ll,mm), ni(ii), an(jj), cn(kk), nu(ll), rho(mm), a1(nn), c1(oo), a2(pp), c2(qq) 
     real n3(pp)


     pi = 3.14159265359

     OPEN(8,file='COLLV2.bin',form='unformatted')
     OPEN(9,file='COLLV2.dat')
     
     DO n=1,nn,1
     	IF(n.le.10)then
     		dr(n) = 1.e-6
     	ELSE IF(n.gt.10.and.n.le.19)THEN
     		dr(n) = 10.e-6
     	ELSE 
     		dr(n) = 100.e-6
     	END IF
     END DO
     
     DO i = 1,ii,1
        ni(i) = 10.**(i-1)   /1000. !L-1 --> m-3 [1,10,100,1000,10000 L-1]

        DO j = 1,jj,1
           an(j) = 10.**(j-1) *1.e-6 !microns --> m

           DO k = 1,kk,1
              cn(k) = 10.**(k-1)  *1.e-6

              DO l = 1,ll,1
                 nu(l) = real(l)

                 DO m = 1,mm,1
                    rho(m) = real(m)*100.

                    coll(i,j,k,l,m) = 0
                    ncoll(i,j,k,l,m) = 0

                    !!These inner loops should have a very fine resolution to make sure we are 
                    !!integrating over "all possible" sizes

                    DO n = 1,nn,1
                       IF(n.eq.1)THEN 
                       	  a1(n) = 1.e-6         
                       ELSE 
                          a1(n) = a1(n-1)+dr(n)
                       END IF

                       DO o = 1,oo,1
                         IF(o.eq.1)THEN
                           c1(o) = 1.e-6         
                         ELSE 
                           c1(o) = c1(o-1)+dr(o)
                         END IF

                         DO p = 1,pp,1
                         	IF(p.eq.1)THEN
                         	  a2(p) = 1.e-6         
                            ELSE 
                              a2(p) = a2(p-1)+dr(p)
                            END IF

                            DO q = 1,qq,1
                         	   IF(q.eq.1)THEN
                         	     c2(q) = 1.e-6         
                               ELSE 
                                 c2(q) = c2(q-1)+dr(q)
                               END IF
        
                               gam = exp(gammln(nu(l)))  
                               maxx = max(a1(n),c1(o))
                               maxy = max(a2(p),c2(q))    
                               area = (maxx + maxy)**2

                               vol1 = a1(n)**2*c1(o)
                               vol2 = a2(p)**2*c2(q)

                               n1 = ni(i)/gam**2 * (a1(n)/an(j))**(nu(l)-1.)*(exp(-a1(n)/an(j))/an(j)) *&
                                                   (c1(o)/cn(k))**(nu(l)-1.)*(exp(-c1(o)/cn(k))/cn(k))
                               n2 = ni(i)/gam**2 * (a2(p)/an(j))**(nu(l)-1.)*(exp(-a2(p)/an(j))/an(j)) *&
                                                   (c2(q)/cn(k))**(nu(l)-1.)*(exp(-c2(q)/cn(k))/cn(k))
!                                 
                               !!Heymsfield, Bansemer, and Schmitt 2013
                               if(2.*maxx .lt. 41) v1 = 0.0028*(2.*maxx)**2.
                               if(2.*maxx .ge. 41 .and. 2.*maxx .lt. 839) v1 = 0.0791*(2.*maxx)**1.101
                               if(2.*maxx .ge. 839) v1 = 62.29*(2.*maxx)**0.1098
                               v1 = v1/100.

                               if(2.*maxy .lt. 41) v2 = 0.0028*(2.*maxy)**2.
                               if(2.*maxy .ge. 41 .and. 2.*maxy .lt. 839) v2 = 0.0791*(2.*maxy)**1.101
                               if(2.*maxy .ge. 839) v2 = 62.29*(2.*maxy)**0.1098
                               v2 = v2/100.

                               mrate(n,o,p,q) = pi*rho(m)*(vol1+vol2)*area*abs(v1-v2)*n1*n2*E
                               nrate(n,o,p,q) = pi*area*abs(v1-v2)*n1*n2*E
                               ncheck(p,q) = n2
                               
                               coll(i,j,k,l,m) = coll(i,j,k,l,m) + mrate(n,o,p,q)*dr(n)*dr(o)*dr(p)*dr(q)
                               ncoll(i,j,k,l,m) = ncoll(i,j,k,l,m) + nrate(n,o,p,q)*dr(n)*dr(o)*dr(p)*dr(q)
                             END DO!c2(q)

                          END DO!a2(p)

                       END DO!c1(o)
                    
                    END DO!a1(n)
                    
                    !coll(i,j,k,l,m) = sum(sum(sum(sum(mrate,dim=4)*dr,dim=3)*dr,dim=2)*dr,dim=1)*dr
                    !ncoll(i,j,k,l,m) = sum(sum(sum(sum(nrate,dim=4)*dr,dim=3)*dr,dim=2)*dr,dim=1)*dr
                    
                    !write(*,*) ni(i),an(j),cn(k),nu(l),rho(m), coll(i,j,k,l,m), ncoll(i,j,k,l,m)
                    write(9,*) ni(i),an(j),cn(k),nu(l),rho(m), coll(i,j,k,l,m), ncoll(i,j,k,l,m)
                    
                 END DO!rho(m)

              END DO!nu(l)

           END DO!cn(k)

        END DO!an(j)

     END DO!ni(i)


     WRITE(8) (ni(i),i=1,ii)
     WRITE(8) (an(j),j=1,jj)
     WRITE(8) (cn(k),k=1,kk)
     WRITE(8) (nu(l),l=1,ll)
     WRITE(8) (rho(m),m=1,mm)
     WRITE(8) (((((coll(i,j,k,l,m),i=1,ii),j=1,jj),k=1,kk),l=1,ll),m=1,mm)
     WRITE(8) (((((ncoll(i,j,k,l,m),i=1,ii),j=1,jj),k=1,kk),l=1,ll),m=1,mm)
     CLOSE(8)
     CLOSE(9)

END PROGRAM COLLECTION_LOOKUP_GENERATOR_BINARY



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

