      program main


      REAL :: coll_ni(5), coll_an(4), coll_cn(4), coll_nu(8), &
                        coll_rho(9)
      REAL :: coll(5,4,4,8,9), ncoll(5,4,4,8,9)
      REAL :: coll_phi(50), coll_r(27)
      REAL :: acoll(50,27), ancoll(50,27), cncoll(50,27), &
                ddcoll(50,27)
      INTEGER :: ii, jj, kk, ll, mm, nn, oo, iii, jjj, kkk, &
                lll, mmm, nnn, ooo

      iii = 5
      jjj = 4
      kkk = 4
      lll = 8
      mmm = 9
      nnn = 50
      ooo = 27

      OPEN(1,FILE=&
        "COLLV2.bin",form=&
        'unformatted')!!Lookup table for aggregation mass and number
      READ(1) (coll_ni(ii),ii=1,iii) !ni = 1, 10, 100, 1000, 10000 L-1
      READ(1) (coll_an(jj),jj=1,jjj) !an = 1, 10, 100, 1000, 10000 um
      READ(1) (coll_cn(kk),kk=1,kkk) !cn = 1, 10, 100, 1000, 10000 um
      READ(1) (coll_nu(ll),ll=1,lll) !nu = 1, 2, 3, 4, 5, 6, 7, 8
      READ(1) (coll_rho(mm),mm=1,mmm)!rho = 100, 200, 300, 400, 500, 600, 700, 800, 900 kg/m3
      READ(1) (((((coll(ii,jj,kk,ll,mm),ii=1,iii),jj=1,jjj),kk=1,kkk),&
                        ll=1,lll),mm=1,mmm)
      READ(1) (((((ncoll(ii,jj,kk,ll,mm),ii=1,iii),jj=1,jjj),kk=1,kkk),&
                        ll=1,lll),mm=1,mmm)
      CLOSE(1)

!      OPEN(1,FILE=&
!        "ACOLLdd.bin",form=&
!        'unformatted')!!Lookup table for aggregation a, an, and cn
!      READ(1) (coll_phi(nn),nn=1,nnn) !phi = 0.01 --> 100.0 logarithmically spaced
!      READ(1) (coll_r(oo),oo=1,ooo)   !r = 1 -> 10, 20 -> 100, 200 -> 1000 microns
!      READ(1) ((acoll(nn,oo),nn=1,nnn),oo=1,ooo)  !a_avg
!      READ(1) ((ancoll(nn,oo),nn=1,nnn),oo=1,ooo) !an
!      READ(1) ((cncoll(nn,oo),nn=1,nnn),oo=1,ooo) !cn
!      READ(1) ((ddcoll(nn,oo),nn=1,nnn),oo=1,ooo) !delta density
!      CLOSE(1)

        
      end program main

