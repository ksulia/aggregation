PROGRAM COLLECTION_LOOKUP_GENERATOR

  USE NETCDF_PLUGIN

  IMPLICIT NONE


     integer :: i, j, k, l, m, n, o, p, q
     integer, parameter :: E = 1.0, nsum = 100
     real, parameter :: dr = 10.e-6
     integer, parameter :: ii=5, jj=4, kk=jj, ll=8, mm=9, nn=nsum, oo=nsum, pp=nsum, qq=nsum
     real gammln, pi, gam, maxx, maxy, area, vol1, vol2, n1, n2, v1, v2
     real ncheck(nn,oo), mrate(nn,oo,pp,qq), nrate(nn,oo,pp,qq)
     real coll(ii,jj,kk,ll,mm), ncoll(ii,jj,kk,ll,mm), ni(ii), an(jj), cn(kk), nu(ll), rho(mm), a1(nn), c1(oo), a2(pp), c2(qq) 
     real n3(pp)

     integer ncid,vid,nid,aid,cid,nuid,rhoid
     character(7) :: fname

     pi = 3.14159265359

     OPEN(1,FILE="COLL.dat")
     DO i = 1,ii,1
        ni(i) = 10.**(i-1)   *1000. !L-1 --> m-3

        DO j = 1,jj,1
           an(j) = 10.**(j-1) *1.e-6 !microns --> m

           DO k = 1,kk,1
              cn(k) = 10.**(k-1)  *1.e-6

              DO l = 1,ll,1
                 nu(l) = real(l)

                 DO m = 1,mm,1
                    rho(m) = real(m)*100.

                    print*,ni(i),an(j),cn(k),nu(l),rho(m)
                    coll(i,j,k,l,m) = 0
                    ncoll(i,j,k,l,m) = 0

                    !!These inner loops should have a very fine resolution to make sure we are 
                    !!integrating over "all possible" sizes

                    DO n = 1,nn,1                  
                       a1(n) = n*dr!!10.**(n-1)  *1.e-6  

                       DO o = 1,oo,1
                         c1(o) = o*dr!!10.**(o-1)   *1.e-6

                         DO p = 1,pp,1
                            a2(p) = p*dr!!10.**(p-1)    *1.e-6

                            DO q = 1,qq,1
                               c2(q) = q*dr!!10.**(q-1)     *1.e-6
        
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
!                               n3(q) = ni(i)/gam * (c2(q)/cn(k))**(nu(l)-1.)*(exp(-c2(q)/cn(k))/cn(k))!(a2(p)/an(j))**(nu(l)-1.)*(exp(-a2(p)/an(j))/an(j))          
                               !!Heymsfield, Bessimer, and Schmitt 2013
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
                             END DO!c2(q)

                          END DO!a2(p)

!			print*, ni(i), sum(sum(ncheck,dim=2)*dr,dim=1)*dr

                       END DO!c1(o)

                    END DO!a1(n)

                    coll(i,j,k,l,m) = sum(sum(sum(sum(mrate,dim=4)*dr,dim=3)*dr,dim=2)*dr,dim=1)*dr  !kg/m3/s
                    ncoll(i,j,k,l,m) = sum(sum(sum(sum(nrate,dim=4)*dr,dim=3)*dr,dim=2)*dr,dim=1)*dr !#/m3/s
                    WRITE(1,*) ni(i), an(j), cn(k), nu(l), rho(m), coll(i,j,k,l,m),ncoll(i,j,k,l,m)
                 END DO!rho(m)

              END DO!nu(l)

           END DO!cn(k)

        END DO!an(j)

     END DO!ni(i)
  CLOSE(1)
  print*,maxval(coll),maxval(ncoll)
  !OUTPUT USING NETCDF ROUTINES
  WRITE(FNAME,"(A7)") "COLL.NC"
  CALL CREATENETCDF(NCID,TRIM(FNAME))
  CALL DEFINEDIMS(NCID,NID,II,"II")  !number concentration
  CALL DEFINEDIMS(NCID,AID,JJ,"JJ")  !characteristic a-axis length
  CALL DEFINEDIMS(NCID,CID,KK,"KK")  !characteristic c-axis length
  CALL DEFINEDIMS(NCID,NUID,LL,"LL") !gamma shape parameter
  CALL DEFINEDIMS(NCID,RHOID,MM,"MM")!icey density
  CALL DEFINEVAR1D(NCID,"NI",NID,VID)
  CALL DEFINEVAR1D(NCID,"AN",AID,VID)
  CALL DEFINEVAR1D(NCID,"CN",CID,VID)
  CALL DEFINEVAR1D(NCID,"NU",NUID,VID)
  CALL DEFINEVAR1D(NCID,"RHOI",RHOID,VID)
  CALL DEFINEVAR5D(NCID,"COLL",NID,AID,CID,NUID,RHOID,VID)
  CALL DEFINEVAR5D(NCID,"NCOLL",NID,AID,CID,NUID,RHOID,VID)
  CALL ENDDEF(NCID)
  CALL CLOSENETCDF(NCID)

  CALL OPENNETCDF(NCID,TRIM(FNAME))
  CALL WRITEVAR1D(NCID,NI,II,"NI")
  CALL WRITEVAR1D(NCID,AN,JJ,"AN")
  CALL WRITEVAR1D(NCID,CN,KK,"CN")
  CALL WRITEVAR1D(NCID,NU,LL,"NU")
  CALL WRITEVAR1D(NCID,RHO,MM,"RHOI")
  CALL WRITEVAR5D(NCID,COLL,II,JJ,KK,LL,MM,"COLL")
  CALL WRITEVAR5D(NCID,NCOLL,II,JJ,KK,LL,MM,"NCOLL")
  CALL CLOSENETCDF(NCID)
!     coll = coll*pi  !!divide by air density in model!!


END PROGRAM COLLECTION_LOOKUP_GENERATOR



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

