      PROGRAM TEST_LOOP

      !USE NETCDF_PLUGIN

      IMPLICIT NONE


      integer :: i, j, k, l, m, n, o, p, q, c
      integer, parameter :: E = 1.0, nsum = 92 
!      real, parameter :: dr = 10.e-6
      real r(nsum), dr(nsum)
      integer, parameter :: ii=5, jj=4, kk=jj, ll=8, mm=9, nn=nsum, oo=nsum, pp=nsum, qq=nsum
      real gammln, pi, gam, maxx, maxy, area, vol1, vol2, n1, n2, v1, v2
      real ncheck(nn,oo), mrate(nn,oo,pp,qq), nrate(nn,oo,pp,qq)
      real coll(ii,jj,kk,ll,mm), ncoll(ii,jj,kk,ll,mm), ni(ii), an(jj), cn(kk), nu(ll), rho(mm), a1(nn), c1(oo), a2(pp), c2(qq) 
      real n3(pp)

      integer ncid,vid,nid,aid,cid,nuid,rhoid
      character(7) :: fname

      pi = 3.14159265359

      r(1) = 0.5
      r(2) = 1.0
      r(3) = 1.5
      dr(1) = 0.5
      dr(2) = r(2) - r(1)
      dr(3) = r(3) - r(2)
      c = 0
      DO n = 4, nn, 1
        c = c + 1
        r(n) = r(n-1) + 5.*10.**floor((c-18.)/18.)
        dr(n) = r(n) - r(n-1)
        print*,floor((c-18.)/18.),5.*10**floor((c-18.)/18.),r(n), dr(n)
      END DO
      

      END PROGRAM TEST_LOOP
