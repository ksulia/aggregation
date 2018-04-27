MODULE MODULE_MP_SULIAHARRINGTON
    
      IMPLICIT NONE
      
      REAL, PARAMETER :: PI = 3.1415926535897932384626434
      REAL, PARAMETER :: SQRTPI = 0.9189385332046727417803297
      
      REAL, PRIVATE ::  rhoi    !BULK DENSITY OF CLOUD ICE
      REAL, PRIVATE ::  nu      !DISTRIBUTION SHAPE FACTOR, ICE
      REAL, PRIVATE ::  nus     !DISTRIBUTION SHAPE FACTOR, SNOW
      REAL, PRIVATE ::  nuc     !ICE NUCLEATION CONCENTRAION (#/L)
      REAL, PRIVATE ::  rd      !GAS CONSTANT OF DRY AIR
      REAL, PRIVATE ::  cp      !SPECIFIC HEAT FOR DRY AIR (CONST P)
      REAL, PRIVATE ::  ao      !INITIAL CHARACTERISTIC A-AXIS LENGTH
      REAL, PRIVATE ::  co      !INITIAL CHARACTERISTIC C-AXIS LENGTH
      REAL, PRIVATE ::  li0     !INITIAL SEMI-MAJOR AXIS LENGTH
      REAL, PRIVATE ::  mi0     !INITIAL PARTICLE MASS
      REAL, PRIVATE ::  gammnu  !GAMMA DIST WITH SHAPE, NU
      REAL, PRIVATE ::  qsmall  !THRES MIN MIXING RATIO VALUE
      REAL, PRIVATE ::  rv      !GAS CONSTANT OF WATER VAPOR
      REAL, PRIVATE ::  bif     !FALL SPEED EXPONENTIAL COEFF
      REAL, PRIVATE ::  aif     !FALL SPEED LEADING COEFF
      REAL, PRIVATE ::  g       !ACCELERATION OF GRAVITY
      REAL, PRIVATE ::  cpw

!     parameters for snow microphysics
      REAL, PRIVATE ::  DCS     !THRESHOLD AGGREGATION RADIUS
      REAL, PRIVATE ::  BS      !SNOW FALLSPEED PARAMETER
      REAL, PRIVATE ::  RHOSN   !SNOW/AGGREGATE DENSITY
      REAL, PRIVATE ::  EII     !ICE AGGREGATION EFFICIENCY
      REAL, PRIVATE ::  DS      !SNOW SIZE DIST PARAM
      REAL, PRIVATE ::  CS
      REAL, PRIVATE ::  f1s
      REAL, PRIVATE ::  f2s
      REAL, PRIVATE ::  lammins, lammaxs

!     parameters for rain microphysics
      REAL, PRIVATE :: f1r      ! VENTILATION PARAMETER FOR RAIN
      REAL, PRIVATE :: f2r      ! " "
      REAL, PRIVATE :: ar       ! FALL SPEED PARAMETER FOR RAIN
      REAL, PRIVATE :: ars       ! FALL SPEED PARAMETER FOR SNOW
      REAL, PRIVATE :: rhosu    ! APPROXIMATE AIR DENSITY NEAR 850MB
      REAL, PRIVATE :: rhow     ! DENSITY OF LIQUID WATER
      REAL, PRIVATE :: br       ! FALLSPEED PARAMETER FOR RAIN
      REAL, PRIVATE :: bc       ! FALLSPEED PARAMETER FOR CLOUD WATER
      REAL, PRIVATE :: ecr
      REAL, PRIVATE :: aimm     ! FREEZING PARAMETER FOR RAIN
      REAL, PRIVATE :: bimm     ! FREEZING PARAMETER FOR RAIN
      REAL, PRIVATE :: lamminr, lammaxr

      REAL, PRIVATE :: cons1    !
      REAL, PRIVATE :: cons3    ! PARAMETER FOR MASS-WEIGHTED MEAN SNOW FALL SPEED
      REAL, PRIVATE :: cons4    ! PARAMETER FOR MASS-WEIGHTED MEAN RAIN FALL SPEED
      REAL, PRIVATE :: cons5    ! PARAMETER FOR NUMBER-WEIGHTED MEAN SNOW FALL SPEED
      REAL, PRIVATE :: cons6
      REAL, PRIVATE :: cons9    ! EVAPORATION PARAMETER FOR RAIN
      REAL, PRIVATE :: cons10
      REAL, PRIVATE :: cons15
      REAL, PRIVATE :: cons20   ! FREEZING PARAMETER FOR RAIN
      REAL, PRIVATE :: cons26   ! FREEZING PARAMETER FOR CLOUD WATER
      REAL, PRIVATE :: cons29
      REAL, PRIVATE :: cons34   ! " "
      REAL, PRIVATE :: cons35
      REAL, PRIVATE :: cons39   ! FREEZING PARAMETER FOR CLOUD WATER
      REAL, PRIVATE :: cons40   ! FREEZING PARAMETER FOR CLOUD WATER
      REAL, PRIVATE :: cons41
      REAL, PRIVATE :: FUDGE


      REAL, PRIVATE :: coll_ni(5), coll_an(4), coll_cn(4), coll_nu(8), coll_rho(9)
      REAL, PRIVATE :: coll(5,4,4,8,9), ncoll(5,4,4,8,9)
      INTEGER, PRIVATE :: ii, jj, kk, ll, mm, iii, jjj, kkk, lll, mmm
      
      
    CONTAINS
      
    SUBROUTINE SULIAHARRINGTON_INIT
        !************************************************************************************
        
      IMPLICIT NONE
      
      rhoi = 920.
      !     rhoi = 500.
      nu = 4.
      nus = 4.
      nuc = 1.
      rd = 287.15
      cp = 1005.
      cpw = 4187.
      ao = 1.e-6
      co = ao
      li0 = ao
      mi0 = 4./3.*pi*rhoi*(li0)**3
      gammnu = exp(gammln(nu))
      qsmall = 1.e-18
      rv = 461.5
      
      bif = 0.99
      aif = 0.81*(1000.)**bif
      
      !     add constants for snow microphysics !kjs 02/10/2015
      DCS = 125.e-6
      BS = 0.41
      DS = 3.
      rhosn = 100.
      CS = rhosn*pi/6.
      EII = 0.3

!     add constants for rain microphysics
      ar = 841.99667
      ars = 11.72
      rhosu = 85000./(287.15*273.15)
      rhow = 997.
      f1r = 0.78
      f2r = 0.308
      f1s = 0.86
      f2s = 0.28
      g = 9.806
      br = 0.8
      bc = 2.0
      ecr = 1.
      bimm = 100.
      aimm = 0.66
      lammins = 1./2000.e-6
      lammaxs = 1./10.e-6
      lamminr = 1./2800.e-6
      lammaxr = 1./20.e-6

      cons1 = gamma(1.+DS)*CS
      cons3 = gamma(4.+BS)/6.
      cons4 = gamma(4.+ br)/6.
      cons5 = gamma(1.+BS)
      cons6 = gamma(1.+BR)
      cons9 = gamma(5./2.+ br/2.)
      cons10 = 20.*pi*pi*rhow*bimm
      cons15 = -1108.*EII*pi**((1.-BS)/3.)*rhosn**((-2.-BS)/3.)/(4.*720)
      cons20 = 20.*pi*pi*rhow*bimm
      cons26 = pi/6.*rhow
      cons29 = 4./3.*pi*rhow*(25.e-6)**3
      cons34 = 5./2.+ br/2.
      cons35 = 5./2.+BS/2.
      cons39 = pi*pi/36.*rhow*bimm
      cons40 = pi/6.*bimm
      cons41 = pi*pi*ecr*rhow

      FUDGE = 0.9999

      iii = 1!5
      jjj = 1!4
      kkk = 3!4
      lll = 8
      mmm = 9

      OPEN(1,FILE="COLL2.bin",form='unformatted')!!Lookup table for aggregation mass and number
      READ(1) (coll_ni(ii),ii=1,iii) !ni = 1, 10, 100, 1000, 10000 L-1
      READ(1) (coll_an(jj),jj=1,jjj) !an = 1, 10, 100, 1000, 10000 um
      READ(1) (coll_cn(kk),kk=1,kkk) !cn = 1, 10, 100, 1000, 10000 um
      READ(1) (coll_nu(ll),ll=1,lll) !nu = 1, 2, 3, 4, 5, 6, 7, 8
      READ(1) (coll_rho(mm),mm=1,mmm)!rho = 100, 200, 300, 400, 500, 600, 700, 800, 900 kg/m3
      READ(1) (((((coll(ii,jj,kk,ll,mm),ii=1,iii),jj=1,jjj),kk=1,kkk),ll=1,lll),mm=1,mmm)
      READ(1) (((((ncoll(ii,jj,kk,ll,mm),ii=1,iii),jj=1,jjj),kk=1,kkk),ll=1,lll),mm=1,mmm)
      CLOSE(1)

    END SUBROUTINE SULIAHARRINGTON_INIT


    SUBROUTINE mp_suliaharrington(           &
         itimestep,th,                       &
         qv,qc,qr,qi,qs,                     &
         nc,nr,ni,ns,                        &
         ai,ci,as,cs,                        &
         rho,pii,p,dt,dz,ht,w,               &
         icedep,icesub,rainevap,snowevap,    &
         snowmelt,snowdep,snowsub,snowaccr,  &
         cloudcond,cloudevap,icemelt,        &
         icenuc,rainfrz,cloudfrz,            &
         phi,rhoice,relh,cplx,phis,rhos,cplxs&
         ,ids, ide, jds, jde, kds, kde       &
         ,ims, ime, jms, jme, kms, kme       &
         ,its, ite, jts, jte, kts, kte	     &
         )
      
      integer :: itimestep, i, j, k
      INTEGER, INTENT(IN) :: IDS,IDE,JDS,JDE,KDS,KDE,    &
           IMS,IME,JMS,JME,KMS,KME,ITS,ITE,JTS,JTE,KTS,KTE
      REAL, DIMENSION(IMS:IME, KMS:KME, JMS:JME), INTENT(INOUT):: QV,TH 
      REAL, DIMENSION(IMS:IME, KMS:KME, JMS:JME) :: NC,QC,NR,QR,NI,QI,NS,QS        
      REAL, DIMENSION(IMS:IME, KMS:KME, JMS:JME) :: ai,ci,as,cs,rho,pii,p,dz,ht,w
      
      REAL, DIMENSION(IMS:IME, KMS:KME, JMS:JME) :: &
           icedep,icesub,rainevap,snowevap,snowmelt,snowdep,snowsub,    &
           snowaccr,cloudcond,cloudevap,icemelt,icenuc,rainfrz,cloudfrz,&
           phi,rhoice,relh,cplx,phis,rhos,cplxs
      
      REAL, DIMENSION(its:ite, kts:kte) ::                            &
           QV2D, QC2D, QR2D, QS2D, QI2D,                              &
           NC2D, NR2D, NS2D, NI2D,                                    &
           AI2D, CI2D, AS2D, CS2D,                                    &
           T2D, RHO2D, P2D,                                           &
           ICEDEP2D, ICESUB2D, RAINEVAP2D, SNOWEVAP2D, SNOWMELT2D,    &
           SNOWDEP2D, SNOWSUB2D, SNOWACCR2D, CLOUDCOND2D, CLOUDEVAP2D,&
           ICEMELT2D, ICENUC2D, RAINFRZ2D, CLOUDFRZ2D,                &
           PHI2D, RHOICE2D, RELH2D, CPLX2D, PHIS2D, RHOS2D, CPLXS2D
      REAL, DIMENSION(its:ite, kts:kte, jts:jte) :: T
      REAL, DIMENSION(kts:kte) :: DZQ
      REAL :: dt

      DO i=its,ite
         DO j=jts,jte
            DO k=kts,kte

               T(i,k,j) = TH(i,k,j) * PII(i,k,j)

            END DO
         END DO
      END DO

      DO j=jts,jte              !(east-west)
         DO i=its,ite           !(north-south)
            DO k=kts,kte        !(vertical)
               
               QV2D(i,k) = QV(i,k,j) !"QV"    "mixing ratio"                           "kg kg-1"
               QC2D(i,k) = QC(i,k,j) !"QC"    "cloud water mixing ratio"               "kg kg-1"
               QR2D(i,k) = QR(i,k,j) !"QR"    "rain water mixing ratio"                "kg kg-1"
               QI2D(i,k) = QI(i,k,j) !"QI"    "cloud ice mixing ratio"                 "kg kg-1"
               QS2D(i,k) = QS(i,k,j) !"QS"    "snow mixing ratio"                      "kg kg-1"

               NC2D(i,k) = NC(i,k,j) !"cloud num mixing ratio"                          "# kg-1"
               NR2D(i,k) = NR(i,k,j) !"rain num mixing ratio"                           "# kg-1"
               NI2D(i,k) = NI(i,k,j) !"ice num mixing ratio"                            "# kg-1"
               NS2D(i,k) = NS(i,k,j) !"snow num mixing ratio"                           "# kg-1"
               AI2D(i,k) = AI(i,k,j) !"ice a-axis length-weighted volume mixing ratio"  "m3 kg-1"
               CI2D(i,k) = CI(i,k,j) !"ice c-axis length-weighted volume mixing ratio"  "m3 kg-1"
               AS2D(i,k) = AS(i,k,j) !"snow a-axis length-weighted volume mixing ratio" "m3 kg-1"
               CS2D(i,k) = CS(i,k,j) !"snow c-axis length-weighted volume mixing ratio" "m3 kg-1"
               T2D(i,k) = T(i,k,j)   !"temperature" "K"

               ICEDEP2D(i,k) = ICEDEP(i,k,j)        !"ICE DEPOSITIONAL RATE MIXING RATIO" "kg kg-1s-1"
               ICESUB2D(i,k) = ICESUB(i,k,j)        !"ICE SUBLIMATIONAL RATE MIXING RATIO" "kg kg-1s-1"
               RAINEVAP2D(i,k) = RAINEVAP(i,k,j)    !"RAIN EVAPORATIONAL RATE MIXING RATIO" "kg kg-1s-1"
               SNOWEVAP2D(i,k) = SNOWEVAP(i,k,j)    !"SNOW EVAP RATE MIX RATIO" "kg kg-1s-1"
               SNOWMELT2D(i,k) = SNOWMELT(i,k,j)    !"CLOUD WATER MIXING RATIO FROM A CU SCHEME" "kg kg-1s-1"
               SNOWDEP2D(i,k) = SNOWDEP(i,k,j)      !"SNOW DEP RATE MIXING RATIO" "kg kg-1s-1"
               SNOWSUB2D(i,k) = SNOWSUB(i,k,j)      !"SNOW SUBLIMATION RATE MIX RATIO" "kg kg-1s-1"
               SNOWACCR2D(i,k) = SNOWACCR(i,k,j)    !"ACCRET OF SNOW BY RAIN" "kg kg-1s-1"
               CLOUDCOND2D(i,k) = CLOUDCOND(i,k,j)  !"CLOUD DROPLET CONDENSATION RATE MIX RATIO" "kg kg-1s-1"
               CLOUDEVAP2D(i,k) = CLOUDEVAP(i,k,j)  !"CLOUD DROPLET EVAP RATE MIX RATIO" "kg kg-1s-1"
               ICEMELT2D(i,k) = ICEMELT(i,k,j)      !"ICE MELTING RATE" "kg kg-1s-1"
               ICENUC2D(i,k) = ICENUC(i,k,j)        !"ICENUC" "ICE HETERO NUC RATE (FREEZING)" "kg kg-1s-1"
               RAINFRZ2D(i,k) = RAINFRZ(i,k,j)      !"RAINFRZ" "RAIN HOMOGENOUS FREEZING RATE" "kg kg-1s-1"
               CLOUDFRZ2D(i,k) = CLOUDFRZ(i,k,j)    !"CLOUD DROPLET FREEZING RATE" "kg kg-1s-1"
               P2D(i,k) = P(i,k,j)                  !"PRESSURE" "Pa"
               RHO2D(i,k) = RHO(i,k,j)              !"AIR DENSITY "kg m-3"
               DZQ(k) = DZ(i,k,j)                   !"LAYER DEPTH" "m"


            END DO
         END DO

         CALL SULIAHARRINGTON_MICRO(T2D, DT, QV2D, QC2D, QR2D, QS2D, &
              QI2D, NC2D, NR2D, NS2D, NI2D, AI2D, CI2D, AS2D, CS2D, P2D,  &
              RHO2D, DZQ,  &
              IMS, IME, JMS, JME, KMS, KME, ITS, ITE, JTS, JTE, KTS, KTE, &
              ITIMESTEP, ICEDEP2D, ICESUB2D, RAINEVAP2D, SNOWEVAP2D,      &
              SNOWMELT2D, SNOWDEP2D, SNOWSUB2D, SNOWACCR2D, CLOUDCOND2D,  &
              CLOUDEVAP2D, ICEMELT2D, ICENUC2D, RAINFRZ2D, CLOUDFRZ2D,    &
              PHI2D, RHOICE2D, RELH2D, CPLX2D,PHIS2D,RHOS2D,CPLXS2D)

         DO i=its,ite
            DO k=kts,kte

               QV(i,k,j) = QV2D(i,k)
               QC(i,k,j) = QC2D(i,k)
               QI(i,k,j) = QI2D(i,k)
               QR(i,k,j) = QR2D(i,k)
               QS(i,k,j) = QS2D(i,k)

               NC(i,k,j) = NC2D(i,k)
               NR(i,k,j) = NR2D(i,k)
               NI(i,k,j) = NI2D(i,k)
               NS(i,k,j) = NS2D(i,k)
               AI(i,k,j) = AI2D(i,k)
               CI(i,k,j) = CI2D(i,k)
               AS(i,k,j) = AS2D(i,k)
               CS(i,k,j) = CS2D(i,k)
               T(i,k,j) = T2D(i,k)
               TH(i,k,j) = T(i,k,j)/PII(i,k,j)

               ICEDEP(i,k,j) = ICEDEP2D(i,k)
               ICESUB(i,k,j) = ICESUB2D(i,k)
               RAINEVAP(i,k,j) = RAINEVAP2D(i,k)
               SNOWEVAP(i,k,j) = SNOWEVAP2D(i,k)
               SNOWMELT(i,k,j) = SNOWMELT2D(i,k)
               SNOWDEP(i,k,j) = SNOWDEP2D(i,k)
               SNOWSUB(i,k,j) = SNOWSUB2D(i,k)
               SNOWACCR(i,k,j) = SNOWACCR2D(i,k)
               CLOUDCOND(i,k,j) = CLOUDCOND2D(i,k)
               CLOUDEVAP(i,k,j) = CLOUDEVAP2D(i,k)
               ICEMELT(i,k,j) = ICEMELT2D(i,k)
               ICENUC(i,k,j) = ICENUC2D(i,k)
               RAINFRZ(i,k,j) = RAINFRZ2D(i,k)
               CLOUDFRZ(i,k,j) = CLOUDFRZ2D(i,k)

               RHOICE(i,k,j) = RHOICE2D(i,k)
               RELH(i,k,j) = RELH2D(i,k)
               PHI(i,k,j) = PHI2D(i,k)
               CPLX(i,k,j) = CPLX2D(i,k)

               RHOS(i,k,j) = RHOS2D(i,k)
               PHIS(i,k,j) = PHIS2D(i,k)
               CPLXS(i,k,j) = CPLXS2D(i,k)
            END DO              !end k loop
         END DO
      END DO
      
      
    END SUBROUTINE mp_suliaharrington
    
    !*****************************************************************************
    SUBROUTINE SULIAHARRINGTON_MICRO(t, dt, qv, qc, qr, qs, qi, nc,            &
      nr, ns, ni, ai, ci, as, cs, p, rho, dzq, ims, ime, jms, jme, kms, kme,   &
      its, ite, jts, jte, kts, kte, itimestep,icedep, icesub, rainevap,        &
      snowevap, snowmelt, snowdep, snowsub, snowaccr, cloudcond,               &
      cloudevap, icemelt, icenuc, rainfrz, cloudfrz, phi, rhoice, relh,        &
      cplx, phis, rhos, cplxs)
!****************************************************************************

      IMPLICIT NONE

      INTEGER, INTENT(IN) ::  &
      ims,ime,jms,jme,kms,kme,&
      its,ite,jts,jte,kts,kte,&
      itimestep

      REAL, DIMENSION(its:ite,kts:kte),INTENT(INOUT) :: &
      qv, qc, qr, qs, qi, nc, nr, ns, ni, ai, ci, as, cs, t
 
      REAL, DIMENSION(its:ite,kts:kte), INTENT(OUT) :: &
      icedep, icesub, rainevap, snowevap, snowmelt, snowdep, snowsub,&
      snowaccr, cloudcond, cloudevap, icemelt, icenuc, rainfrz, cloudfrz

      REAL, DIMENSION(its:ite,kts:kte), INTENT(OUT) :: &
           rhoice,phi,relh,cplx,&
           rhos,phis,cplxs
      
      REAL, DIMENSION(its:ite,kts:kte), INTENT(IN) :: p, rho
      REAL, DIMENSION(kts:kte), INTENT(IN) :: dzq
      REAL, INTENT(IN) :: dt


      INTEGER i, k, iflag, nstep, n, iaspect, homofreeze, &
      ipart, resupply, sphrflag, SEDON, RAINON,       &
      DSTRCHECKS, ICE_CALCS, EVOLVE_ON, snowflag, ice_start_time,    &
      processes, LTRUE, CNUM, redden, agg_flag

!     process rates: calculates changes due to nucleation,deposition,cond/evap
      REAL mnuccd               !change in qi freezing aerosol
      REAL nnuccd               !change in ni freezing aerosol
      REAL prd, ard, crd        !change in qi, ai, and ci, deposition of ice
      REAL pcc                  !condendation/evaporation cloud droplets
      REAL prc                  !autoconversion of droplets
      REAL pra                  !accretion of droplets by rain
      REAL pre                  !evaporation of rain
      REAL nprci, prci          !change in ni, qi autoconversion cloud ice by snow
      REAL psmlt                !change in q melting snow to rain
      REAL evpms                !change in q melting snow evaporating (snow -->melt-->evap)
      REAL pracs                !change in q rain-snow collection
      REAL prai                 !change in q accretion of cloud ice by snow
      REAL npra                 !change in n due to droplet accretion by rain
      REAL nprai                !change in n accretion of cloud ice
      REAL nprc                 !change in nc autoconversion droplets
      REAL nprc1                !change in nr autoconversion droplets
      REAL nsmlts               !change in n melting snow
      REAL nsmltr               !change in n melting snow to rain
      REAL prds, eprds          !change in q deposition/sublimation snow
      REAL nsubr, nsubs         !loss of nr,ns during evap,sub
      REAL agg, nagg, nsagg, nragg !change in q and n aggregation
      REAL nnew, prd1

!     ice characteristics
      REAL rhobar

      REAL dum, dum1, dum2, temp, theta, celsius, press, rhoa, cpm, qsdum
      REAL evs, evi, qvs, qvi, qvv, xxls, xxlv, xxlf, dqsidt, abi, epss
      REAL ani, cni, rni, anf, cnf, rnf, phii, phif
      REAL deltastr, alphstr, alphv, betam
      REAL iwci, iwcf, vi, nidum, voltmp
      REAL weight, igr1, igr2, igr, losstot
      REAL vtbarb, vtbarbm, vtbarblen, vtdumm, vtdumn
      REAL alpham, alpha_v, beta_v, alpha_a, beta_a, lmean, cap
      REAL qt_adv, qt_sed, qi_frac, qiold, niold
      REAL*8 sui, sup, qvqvs, qvqvsi

      !snow characteristics (aggregate microphysics)
      REAL ans, cns, rns, rhobars, deltastrs, vs, swci, swcf, phisf, rnsf, &
           vtbarbs, vtbarbms, vtbarblens, ards, ansf, crds, cnsf, phiis

      REAL lammin, lammax, lamc, pgam

      REAL, DIMENSION(its:ite,kts:kte) :: vtrmi1, vtrni1, vtrli1, &
      effi1, vtrmc

!     sedimenation terms
      REAL, DIMENSION(kts:kte) ::                                      &
      fc, fr, fs, fi, fnc, fnr, fns, fni, fci, fai,                    &
      dumc, dumr, dums, dumi, dumnc, dumnr, dumns, dumni, dumai, dumci,&
      falloutc, falloutr, fallouts, fallouti, falloutnc,               &
      falloutnr, falloutns, falloutni, falloutai, falloutci,           &
      qcsten, qrsten, qssten, qisten, ncsten, nrsten, nssten, nisten,  &
      aisten, cisten, arn, acn, asn,                                   &
      qiloss2,qcloss,qiloss,qrloss,qcloss2,qrloss2
      
      REAL, DIMENSION(its:ite,kts:kte) :: qsum
      REAL ums, umr, uns, unr
      
      REAL, DIMENSION(kts:kte) :: thetal,dthl

      REAL falltndc, falltndr, falltnds, falltndi, falltndnc, &
      falltndnr, falltndns, falltndni, falltndai, falltndci,  &
      rgvm, precrt, snowrt

      REAL ncloud, cdist1
      REAL  qsum2, dthlmax
      INTEGER  level, level2, bottomflag,topflag

!     parameters for rain microphysics

      REAL n0rr, dv, lamr, mu, kap, sc, ab, epsr, dqsdt, ratio, mnuccr
      REAL nnuccr, lams, n0s

!     parameters for cloud water processes

      REAL mnuccc, nnuccc

      REAL alphap,alphamsp,betamp

      REAL gammaFindSumTrip(60)
      DATA gammaFindSumTrip/  0.910547, 0.81807, 0.6874, 0.60127, 1.59767&
      ,2.32423,  2.08818,  1.61921, 1.15865, 0.863071, 0.617586, 0.453917&
      ,0.351975, 0.28794, 0.269298, 0.28794, 0.333623, 0.418883,  0.56992&
      ,0.796458, 1.14325,  1.64103, 1.90138,  1.82653, 1.61921,   1.47436&
      ,1.32463,  1.25556,  1.22239,   1.206,  1.11522,  1.10751,  1.10738&
      ,1.11484,  1.12234,  1.12221, 1.14529,  1.16884,  1.20104,  1.22573&
      ,1.25094,  1.27666,  1.31183,  1.3388,  1.35704,  1.37553,  1.38479&
      ,1.39411,  1.40349,  1.41294, 1.42245,  1.43202,  1.44166,  1.45137&
      ,1.46114,  1.47097,  1.48087, 1.50105,  1.50087,  1.51098/

      resupply=0
      ipart=1                   !0 for RAMS dendrites and needles (Walko et al., 1995)
!     1 for RAMS plates and columns (Walko et al., 1995)
!     2 for Mitchell's mass relations (hex plates, columns)
!     3 for Mitchell's mass relations (stellars, columns)
      !     4 for Wood's (2007) mass-size relations
      iaspect    = 0            !set constant aspect ratio (sensitivity study)
      sphrflag   = 0            !all ice assumed spheres
      redden = 0
      homofreeze = 1            !homogeneous freezing
      snowflag   = 2            !snow calculations, 0 = all snow off, 1 = only old snow, 2 = only new snow
      SEDON      = 1            !sedimentation
      EVOLVE_ON  = 1            !depositional growth
      RAINON     = 1            !rain processes
      ICE_CALCS  = 1            !all ice calculations
      ice_start_time = 0.!60.*60.*4.0 !time to begin ice nucleation & homogeneous freezing
      LTRUE = 0
      CNUM = 0                  !0 for constant droplet number
      processes = 1

      IF(sphrflag  .eq. 1) rhoi = 920.
      If(redden .eq. 1) rhoi = 500.
  
      alphap = 0.0
      alphamsp = (4./3.)*pi*rhoi
      betamp = 3.0
      
      precrt=0.0
      snowrt=0.0
      qt_adv=0.0
      qt_sed=0.0

      DO i = its,ite
         DO k = kts,kte

            !ICE
            IF(qi(i,k).gt.qsmall.and.ni(i,k).gt.qsmall)THEN
               ni(i,k) = max(ni(i,k),qsmall)
               ai(i,k) = max(ai(i,k),qsmall)
               ci(i,k) = max(ci(i,k),qsmall)

               ani = ((ai(i,k)**2)/(ci(i,k)*nu*ni(i,k)))**(1./3.)
               cni = ((ci(i,k)**2)/(ai(i,k)*nu*ni(i,k)))**(1./3.)
               ai(i,k) = nu*ni(i,k)*ani
               ci(i,k) = nu*ni(i,k)*cni
            ELSE               
               qi(i,k)=0.
               ni(i,k)=0.
               ai(i,k)=0.
               ci(i,k)=0.
            END IF
            
            !CLOUD WATER
            IF(qc(i,k).lt.qsmall)THEN
               qc(i,k) = 0.
               nc(i,k) = 0.
            END IF

            !RAIN WATER
            IF(qr(i,k).lt.qsmall)THEN
               qr(i,k) = 0.
               nr(i,k) = 0.
            END IF

            !SNOW (AGGREGATES)

            IF(snowflag .eq. 2.and.qs(i,k).gt.qsmall.and.ns(i,k).gt.qsmall)THEN
               ns(i,k) = max(ns(i,k),qsmall)
               as(i,k) = max(as(i,k),qsmall)
               cs(i,k) = max(cs(i,k),qsmall)
               
               ans = ((as(i,k)**2)/(cs(i,k)*nus*ns(i,k)))**(1./3.)
               cns = ((cs(i,k)**2)/(as(i,k)*nus*ns(i,k)))**(1./3.)
               as(i,k) = nus*ns(i,k)*ans
               cs(i,k) = nus*ns(i,k)*cns
            ELSE
               qs(i,k) = 0.
               ns(i,k) = 0.
               as(i,k) = 0.
               cs(i,k) = 0.
            END IF

            
!     set droplet concentration (units of kg-1)
            IF(CNUM.eq.0)THEN
               nc(i,k)=200.e6/rho(i,k)
            END IF

            !     initialize process rates
            mnuccd = 0.
            nnuccd = 0.
            mnuccr = 0.
            nnuccr = 0.
            mnuccc = 0.
            nnuccc = 0.
            prd = 0.
            ard = 0.
            crd = 0.
            pcc = 0.
            prc = 0.
            pra = 0.
            prds = 0.
            eprds = 0.
            nsubr = 0.
            nsubs = 0.
            pre = 0.
            npra = 0.
            nprci = 0.
            nprai = 0.
            nprc = 0.
            nprc1 = 0.
            prci = 0.
            agg = 0.
            nagg = 0.
            nsagg = 0.
            nragg = 0.
            psmlt = 0.
            pracs = 0.
            evpms = 0.
            nsmlts = 0.
            nsmltr = 0.
            qrsten(k)  = 0.0
            qcsten(k)  = 0.0
            qssten(k) = 0.0
            qisten(k) = 0.0
            ncsten(k) = 0.0
            nrsten(k) = 0.0
            nssten(k) = 0.0
            nisten(k) = 0.0
            aisten(k) = 0.0
            cisten(k) = 0.0
            qsum(i,k) = 0.0
            qiloss(k) = 0.0
            qrloss(k) = 0.0
            qcloss(k) = 0.0
            qiloss2(k)=0.0
            qcloss2(k)=0.0
            qrloss2(k)=0.0
            phi(i,k)=1.0

            thetal(k)=0.0
            dthl(k)=0.0

!     initialize ice characteristics in case first time ice nucleation

            deltastr = 1.
            rhobar = 920.
            If(redden .eq. 1) rhoi = 500.
            temp = t(i,k)
            celsius = temp-273.15
            press = p(i,k)
            rhoa = rho(i,k)
            theta = temp*(100000./press)**(rd/cp)

!     set parameters that vary in space and time
!     add fix for low pressure
            evs = min(FUDGE*press,polysvp1(temp,0)) !Pa
            evi = min(FUDGE*press,polysvp1(temp,1)) !pa
            qvs = 0.622*evs/(press-evs)
            qvi = 0.622*evi/(press-evi)
            qvqvs = qv(i,k)/qvs
            qvqvsi = qv(i,k)/qvi
            qvv = qvs/qvi - 1. 
            sup = qv(i,k)/qvs-1.
            relh(i,k) = qv(i,k)/qvs*100.
            sui = qv(i,k)/qvi-1.
            xxls = 3.15e6-2370.*temp+0.3337e6
            xxlv = 3.1484E6-2370.*temp
            xxlf = xxls-xxlv
            dqsdt = xxlv*qvs/(rv*temp**2)
            dqsidt = xxls*qvi/(rv*temp**2)
            cpm = cp*(1.+0.887*qv(i,k))
            abi = 1.+dqsidt*xxls/cpm
!     diffusivity of water vapor
            dv = 8.794E-5*temp**1.81/press
!     viscocity of air
            mu = (1.496E-6*temp**1.5/(temp+120.))
            kap = 1.414E3*mu

!     schmidt
            sc = mu/(rho(i,k)*dv)

!     parameter for rain evaporation
            ab = 1.+dqsdt*xxlv/cpm

!     "a" parameter for rain fallspeed
            arn(k) = (rhosu/rho(i,k))**0.54*ar
!     "a" parameter for cloud droplet fallspeed, based on Stokes Law
            acn(k) = g*rhow/(18.*mu)
!     "a" parameter for snow fall speed
            asn(k) = (rhosu/rho(i,k))**0.54*ars

            IF(celsius .le. -1. .and. celsius .ge. -59.)THEN 
               weight = (ABS(real(int(celsius))) + 1.0) - ABS(celsius)
               igr1 = gammaFindSumTrip(int(celsius)*(-1))
               igr2 = gammaFindSumTrip((int(celsius)*(-1))+1)
               igr = weight*igr1 + (1.-weight)*igr2    
            ELSE IF(celsius .lt. -60.)THEN
               igr = 1.5
            ELSE
               igr = 1.0
            END IF

            
            IF(qc(i,k).lt.qsmall.and.qi(i,k).lt.qsmall.and.&
               qs(i,k).lt.qsmall.and.qr(i,k).lt.qsmall)THEN
               IF((temp.lt.273.15.and.qvqvsi.lt.FUDGE).or.&
                  (temp.ge.273.15.and.qvqvs.lt.FUDGE))THEN
                  Print*,'Time = ', itimestep, 'there is no condensate'
                  GOTO 200
               END IF
            END IF
!     RAIN PROCESSES ----------------------------------------------------------------------
!     calculate rain size distribution parameters
            nr(i,k) = max(0.,nr(i,k))
            nc(i,k) = max(0.,nc(i,k))
            ns(i,k) = max(0.,ns(i,k))
            IF(RAINON.eq.1) THEN
               IF(qr(i,k).ge.qsmall)THEN
                  lamr = (pi*rhow*nr(i,k)/qr(i,k))**0.3333
                  n0rr = nr(i,k)*lamr
                  
                  IF(lamr .lt. lamminr)THEN
                     lamr = lamminr
                     n0rr = lamr**4*qr(i,k)/(pi*rhow)
                     nr(i,k) = n0rr/lamr
                  ELSE IF(lamr.gt.lammaxr)THEN
                     lamr = lammaxr
                     n0rr = lamr**4*qr(i,k)/(pi*rhow)
                     nr(i,k) = n0rr/lamr
                  END IF
               END IF
            END IF
            !     CLOUD DROPLETS ---------------------------------------------------------------------
!     Martin et al. (1994) formula for pgam
            
            IF(qc(i,k).gt.qsmall)THEN
               dum = press/(287.15*temp)
               pgam = 0.0005714*(nc(i,k)/1.e6*dum)+0.2714
               pgam = 1./(pgam**2)-1.
               pgam = max(pgam,2.)
               pgam = min(pgam,10.)
               
               lamc = (cons26*nc(i,k)*GAMMA(pgam+4.)/&
               (qc(i,k)*GAMMA(pgam+1.)))**0.3333
               
               lammin = (pgam+1.)/60.e-6
               lammax = (pgam+1.)/1.e-6
               
               IF(lamc.lt.lammin)THEN
                  lamc = lammin
                  nc(i,k) = exp(3.*log(lamc)+log(qc(i,k))+&
                  log(GAMMA(pgam+1.))-log(GAMMA(pgam+4.)))/cons26
               ELSE IF(lamc.gt.lammax)THEN
                  lamc = lammax
                  nc(i,k) = exp(3.*log(lamc)+log(qc(i,k))+&
                  log(GAMMA(pgam+1.))-log(GAMMA(pgam+4.)))/cons26
               END IF
            END IF              !qc >=qsmall
            !     SNOW---------------------------------------------------------------------------------
            IF(snowflag .eq. 1)THEN
               IF(qs(i,k).ge.qsmall)THEN
                  lams = (cons1*ns(i,k)/qs(i,k))**(1./DS)
                  n0s = ns(i,k)*lams
                  
                  IF(lams.lt.lammins)THEN
                     lams = lammins
                     n0s = lams**4*qs(i,k)/cons1
                     ns(i,k) = n0s/lams
                  ELSE IF(lams.gt.lammaxs)THEN
                     lams = lammaxs
                     n0s = lams**4*qs(i,k)/cons1
                     ns(i,k) = n0s/lams
                  END IF
                  qsdum = 2.*pi*n0s*rho(i,k)*(f1s/(lams*lams)+&
                       f2s*SQRT(asn(k)*rho(i,k)/mu)*sc**0.3333*cons10/&
                       (lams**cons35))
               END IF           !qs >= qsmall
            END IF!snowflag
                 
            IF(temp.ge.273.15)THEN
               !     autoconversion of cloud liquid water to rain
               IF(RAINON.eq.1)THEN
                  IF(qc(i,k).ge.1.E-6)THEN
                     !     formula from khairoutdinov and kogan 2000, mwr
                     prc=1350.*qc(i,k)**2.47*(nc(i,k)/1.e6*&
                          rho(i,k))**(-1.79)
                     
                     nprc1 = prc/cons29
                     nprc = prc/(qc(i,k)/nc(i,k))
                     
                     nprc = min(nprc,nc(i,k)/dt)
                     nprc1 = min(nprc1,nprc)
                  END IF
                  !     accretion of cloud water by rain
                  !     continuous collection eq with grav. coll. kernel, droplet fall spd neglected
                  
                  IF(qr(i,k).ge.1.E-8 .and. qc(i,k).ge.1.E-8)THEN
                     dum=(qc(i,k)*qr(i,k))
                     pra = 67.*dum**1.15
                     npra = pra/qc(i,k)*nc(i,k)
                  END IF
                  
                  !     self-collection of rain drops
                  IF(qr(i,k).ge.1.e-8)THEN
                     dum1 = 300.e-6
                     IF(1./lamr.lt.dum1)THEN
                        dum = 1.
                     ELSE IF(1./lamr.ge.dum1)THEN
                        dum = 2.-exp(2300.*(1./lamr-dum1))
                     END IF
                     nragg = -5.78*dum*nr(i,k)*qr(i,k)*rho(i,k)
                  END IF
                  !     evaporation of rain (RUTLEDGE AND HOBBS 1983)
                  
                  IF(qr(i,k).ge.qsmall)THEN
                     epsr = 2.*pi*n0rr*rho(i,k)*dv*(f1r/(lamr*lamr)+f2r*&
                          SQRT(arn(k)*rho(i,k)/mu)*sc**0.33333*cons9/&
                          (lamr**cons34))
                  ELSE
                     epsr = 0.
                  END IF
                  !     no condensation onto rain, only evap allowed
                  
                  IF (qv(i,k).lt.qvs) THEN
                     pre = min((epsr*(qv(i,k)-qvs)/ab),0.)
                  ELSE
                     pre = 0.
                  END IF
                  !     collection of snow by rain above freezing
                  !     formula from Ikawa and Saito (1991)
                  IF(snowflag.eq.1)THEN
                     IF(qr(i,k).ge.1.e-8.and.qs(i,k).ge.1.e-8)THEN
                        ums = asn(k)*cons3/(lams**bs)
                        umr = arn(k)*cons4/lamr**br
                        uns = asn(k)*cons5/lams**bs
                        unr = arn(k)*cons6/lamr**br
                        !     set realistic limits on fallspeeds
                        
                        dum = (rhosu/rho(i,k))**0.54
                        ums = min(ums,1.2*dum)
                        uns = min(uns,1.2*dum)
                        umr = min(umr,9.1*dum)
                        unr = min(unr,9.1*dum)
                        
                        pracs = cons41*(SQRT((1.2*umr-0.95*ums)**2 + 0.08*&
                             ums*umr)*rho(i,k)*n0rr*n0s/lamr**3* &
                             (5./(lamr**3*lams)+2./(lamr*lamr*lams*lams)+0.5/&
                             (lamr*lams**3)))
                     END IF        !qr & qs >=1.e-8
                     
                     !     melting of snow
                     !     snow may persist above freezing, from Rutledge and Hobbs (1984)
                     !     if supersat, snow melts to form rain
                     IF(qs(i,k).ge.1.e-8)THEN
                        dum = -cpw/xxlf*(temp-273.15)*pracs
                        psmlt = qsdum*kap*(273.15-temp)/xxlf+dum
                        IF(qvqvs.lt.1.)THEN
                           epss = qsdum*dv
                           evpms = (qv(i,k)-qvs)*epss/ab
                           evpms = max(evpms,psmlt)
                           psmlt = psmlt-evpms
                        END IF
                     END IF        !qs>=1.e-8
                     pracs = 0.
                  END IF!snowflag
               END IF!RAINON
               !     check cloud water conservation, update qc and qr with process rates
               
               dum = (prc+pra)*dt
               IF(dum.gt.qc(i,k).and.qc(i,k).ge.qsmall)THEN
                  ratio = qc(i,k)/dum
                  prc = prc*ratio
                  pra = pra*ratio
               END IF
               
!     conservation of qr, NOTE: pre is a negative number

               dum = (-pre-pra-prc+psmlt)*dt
               IF(dum.gt.qr(i,k).and.qr(i,k).ge.qsmall) THEN
                  ratio = (qr(i,k)/dt+pra+prc+pracs-psmlt)/(-pre)
                  pre = pre*ratio
               END IF
               
               !     conservation of qs
               IF(snowflag .eq. 1)THEN
                  dum = (pracs-evpms-psmlt)*dt !melting, evap, & accretion of snow
                  IF(dum.gt.qs(i,k).and.qs(i,k).ge.qsmall)THEN
                     ratio = qs(i,k)/dum
                     pracs = pracs*ratio
                     psmlt = psmlt*ratio
                     evpms = evpms*ratio
                  END IF
               END IF!snowflag
            
            ELSE !temp
               IF (RAINON.eq.1)THEN
                  !     autoconversion of cloud liquid water to rain
                  IF(qc(i,k).ge.1.E-6)THEN
                     !     formula from khairoutdinov and kogan 2000, mwr
                     prc=1350.*qc(i,k)**2.47*(nc(i,k)/1.e6*&
                          rho(i,k))**(-1.79)
                     
                     nprc1 = prc/cons29
                     nprc = prc/(qc(i,k)/nc(i,k))
                     
                     nprc = min(nprc,nc(i,k)/dt)
                     nprc1 = min(nprc1,nprc)
                  END IF
                  
                  !     accretion of cloud water by rain
                  !     continuous collection eq with grav. coll. kernel, droplet fall spd neglected
                  
                  IF(qr(i,k).ge.1.E-8 .and. qc(i,k).ge.1.E-8)THEN
                     dum=(qc(i,k)*qr(i,k))
                     pra = 67.*dum**1.15
                     npra = pra/qc(i,k)*nc(i,k)
                  END IF
                  
                  !     self-collection of rain drops
                  IF(qr(i,k).ge.1.e-8)THEN
                     dum1 = 300.e-6
                     IF(1./lamr.lt.dum1)THEN
                        dum = 1.
                     ELSE IF(1./lamr.ge.dum1)THEN
                        dum = 2.-exp(2300.*(1./lamr-dum1))
                     END IF
                     nragg = -5.78*dum*nr(i,k)*qr(i,k)*rho(i,k)
                  END IF
!     evaporation of rain (RUTLEDGE AND HOBBS 1983)
                  
                  IF(qr(i,k).ge.qsmall)THEN
                     epsr = 2.*pi*n0rr*rho(i,k)*dv*(f1r/(lamr*lamr)+f2r*&
                          SQRT(arn(k)*rho(i,k)/mu)*sc**0.33333*cons9/&
                          (lamr**cons34))
                  ELSE
                     epsr = 0.
                  END IF
               
!     no condensation onto rain, only evap allowed
               
                  IF (qv(i,k).lt.qvs) THEN
                     pre = min((epsr*(qv(i,k)-qvs)/ab),0.)
                  ELSE
                     pre = 0.
                  END IF
               !               end if
               END IF!RAINON
               IF(snowflag.eq.1)THEN
                  !     deposition of qs
                  !               IF(qs(i,k).ge.qsmall)THEN
                  !                  epss = qsdum*dv
                  !                  prds = epss*(qv(i,k)-qvi)/abi
                  !               END IF      
                  !               dum = (qv(i,k)-qvi)/dt
                  !               IF((dum.gt.0..and.prds.gt.dum*FUDGE).or.&
                  !                  (dum.lt.0..and.prds.lt.dum*FUDGE))THEN
                  !                  prds = FUDGE*dum
                  !               END IF
                  !               IF(prds.lt.0.)THEN
                  !                  eprds = prds
                  !                  prds = 0.
                  !               END IF
               
                  !     aggregation of qs
                  IF(qs(i,k).ge.1.e-8)THEN
                     nsagg = cons15*asn(k)*(qs(i,k)*rho(i,k))**&
                          ((2.+BS)/3.)*(ns(i,k)*rho(i,k))**((4.-BS)/3.)&
                          /rho(i,k)
                  END IF
               END IF!snowflag
               
!     check cloud water conservation, update qc and qr with process rates
            
               dum = (prc+pra)*dt
               IF(dum.gt.qc(i,k).and.qc(i,k).ge.qsmall)THEN
                  ratio = qc(i,k)/dum
                  prc = prc*ratio
                  pra = pra*ratio
               END IF
               
!     conservation of qr, NOTE: pre is a negative number
               
               dum = (-pre-pra-prc)*dt
               IF(dum.gt.qr(i,k).and.qr(i,k).ge.qsmall) THEN
                  ratio = (qr(i,k)/dt+pra+prc)/(-pre)
                  pre = pre*ratio                  
               END IF

               IF(pre.lt.0.)THEN
                  dum = pre*dt/qr(i,k)
                  dum = max(-1.,dum)
                  nsubr = dum*nr(i,k)/dt
               END IF
               
!     conservation of qs
               dum = (-prds-eprds)*dt !deposition and sublimation of snow
               IF(dum.gt.qs(i,k).and.qs(i,k).ge.qsmall)THEN
                  ratio = (qs(i,k)/dt+prds)/(-eprds)
                  eprds = eprds*ratio
               END IF

               IF(eprds.lt.0.)THEN
                  dum = eprds*dt/qs(i,k)
                  dum = max(-1.,dum)
                  nsubs = dum*ns(i,k)/dt
               END IF

            END IF              !temp<273.15
            IF(qi(i,k).gt.qsmall.and.ni(i,k).gt.qsmall)THEN
               
!     set minimum values for ni,ai,ci, otherwise
               
               ni(i,k) = max(ni(i,k),qsmall)
               ai(i,k) = max(ai(i,k),qsmall)
               ci(i,k) = max(ci(i,k),qsmall)
               
!     get characteristic ci,ai (i.e. cni,ani) from ci,ai
               
               ani = ai(i,k)/(nu*ni(i,k))
               cni = ci(i,k)/(nu*ni(i,k))

               
               !check that deltastr, rhobar, and rni are within reasonable bounds
               CALL ICE_CHECKS(1,ni(i,k),qi(i,k),ani,cni,rni,deltastr,rhobar,&
                    iaspect,sphrflag,redden,betam,alphstr,alphv)

               ci(i,k)=nu*ni(i,k)*cni
               ai(i,k)=nu*ni(i,k)*ani

               
!     get iwc to calculate iwc tendency
!     by substracting final and initial values
               vi = 4./3.*pi*rni**3.*exp(gammln(nu+deltastr+2.))/gammnu 
               iwci = ni(i,k)*rhobar*vi*rho(i,k)
               
!     calculate the particle inherent growth ratio (IGR) based on temp
!     igr=1 for ice particles pre-diagnosed as spheres


               IF(iaspect .eq. 1) igr = .27
               IF(sphrflag .eq. 1) igr=1.0
               If(redden .eq. 1) rhoi = 500.
!     calculate number concentration from number mixing ratio

               IF(EVOLVE_ON .eq. 1) THEN
                  
                  nidum = ni(i,k)*rho(i,k)
                  
                  CALL EVOLVE(1,ani,nidum,sui,sup,qvv,temp,press,igr,dt,iwcf,&
                  cnf,iwci,phii,phif,cni,rni,rnf,anf,deltastr,mu,&
                  rhobar,vtbarb,vtbarbm,alphstr,vtbarblen,rhoa,i,k,iaspect&
                  ,ipart,sphrflag,redden,itimestep) 
                  
                  betam = deltastr + 2.0
                  alphstr = co/ao**(deltastr)
                  alphv = 4./3.*pi*alphstr
                  
                  phi(i,k)=phif
                  
!     get deposition/sublimation process rates for qi, ai, and ci
!     deposition rate for ice [=] kg/kg/s (mixing ratio rate)
                  
                  prd=(iwcf-iwci)/rho(i,k)/dt              
                  ard=(anf-ani)*nu*ni(i,k)/dt
                  crd=(cnf-cni)*nu*ni(i,k)/dt 
               END IF           !EVOLVE_ON
            END IF              !ice is present
            IF(snowflag .eq.2)THEN
               IF(qs(i,k).gt.qsmall.and.ns(i,k).gt.qsmall)THEN
               
                  ns(i,k) = max(ns(i,k),qsmall)
                  as(i,k) = max(as(i,k),qsmall)
                  cs(i,k) = max(cs(i,k),qsmall)
                  
                  ans = as(i,k)/(nus*ns(i,k))
                  cns = cs(i,k)/(nus*ns(i,k))
       
                  !check that deltastr, rhobar, and rni are within reasonable bounds
                  CALL ICE_CHECKS(2,ns(i,k),qs(i,k),ans,cns,rns,deltastrs,rhobars,&
                       iaspect,sphrflag,redden,betam,alphstr,alphv)
                  
                  cs(i,k)=nus*ns(i,k)*cns
                  as(i,k)=nus*ns(i,k)*ans
                  
                  
                  !     get iwc to calculate iwc tendency
                  !     by substracting final and initial values
                  vs = 4./3.*pi*rns**3.*exp(gammln(nus+deltastrs+2.))/exp(gammln(nus)) 
                  swci = ns(i,k)*rhobars*vs*rho(i,k)
                  
                  IF(iaspect .eq. 1) igr = .27
                  IF(sphrflag .eq. 1) igr=1.0
                  If(redden .eq. 1) rhobars = 100.
                  !     calculate number concentration from number mixing ratio
                  
                  IF(EVOLVE_ON .eq. 1) THEN
                     
                     nidum = ns(i,k)*rho(i,k)
                     
                     CALL EVOLVE(2,ans,nidum,sui,sup,qvv,temp,press,igr,dt,swcf,&
                          cnsf,swci,phiis,phisf,cns,rns,rnsf,ans,deltastrs,mu,&
                          rhobars,vtbarbs,vtbarbms,alphstr,vtbarblens,rhoa,i,k,iaspect&
                          ,ipart,sphrflag,redden,itimestep) 
                     
                     betam = deltastr + 2.0
                     alphstr = co/ao**(deltastr)
                     alphv = 4./3.*pi*alphstr
                     
                     phis(i,k)=phisf
                     
                     !     get deposition/sublimation process rates for qi, ai, and ci
                     !     deposition rate for ice [=] kg/kg/s (mixing ratio rate)
                     
                     prds=(swcf-swci)/rho(i,k)/dt              
                     ards=(ansf-ans)*nus*ns(i,k)/dt
                     crds=(cnsf-cns)*nus*ns(i,k)/dt 
                  END IF           !EVOLVE_ON
               END IF              !ice is present
            END IF!snowflag
!     get ice nucleation 
            IF(temp.lt.268.15 .and.sui.ge.0.05)THEN
               IF(REAL(itimestep)*dt.gt.ice_start_time)THEN 
              
                  dum = nuc*1000./rho(i,k) !1/L*1000/rho ~ 1000/kg
                  IF(ni(i,k).lt.dum) nnuccd=(dum-ni(i,k))/dt
  
                  mnuccd = nnuccd*mi0
               END IF!ice start time
            END IF!temp and sui
            
!     make sure doesn't push into subsat or supersat

            iflag=0
            prd1=prd+prds
            IF(prd1.lt.0.0.and.sui.ge.0.0)THEN
               IF(prd1.lt.-FUDGE*sui*qvi/abi/dt)THEN
                  prd = -FUDGE*sui*qvi/abi/dt*(prd/prd1)
                  prds = -FUDGE*sui*qvi/abi/dt*(prds/prd1)
                  iflag = 1
               END IF
            END IF 
            IF(prd1.lt.0.0.and.sui.lt.0.0)THEN
               IF(prd1.lt.FUDGE*sui*qvi/abi/dt)THEN
                  prd = FUDGE*sui*qvi/abi/dt*(prd/prd1)
                  prds = FUDGE*sui*qvi/abi/dt*(prds/prd1)
                  iflag = 1
               END IF
            END IF
            IF((prd1.gt.0..and.prd1.gt.FUDGE*sui*qvi/abi/dt))THEN
               prd=FUDGE*sui*qvi/abi/dt*(prd/prd1)
               prds=FUDGE*sui*qvi/abi/dt*(prds/prd1)
               iflag=1
            END IF

          
!     conservation checks
            prd=max(prd,-qi(i,k)/dt)
!     add growth tendencies to get updated variables
            qi(i,k)=qi(i,k)+(prd)*dt
            ai(i,k)=ai(i,k)+(ard)*dt
            ci(i,k)=ci(i,k)+(crd)*dt

            IF(prd .gt. 0.0)THEN
               icedep(i,k) = prd
            ELSE
               icesub(i,k) = prd
            ENDIF


            !     if sublimation, reduce crystal number,
            !     NOTE: this should not impact rhobar, since
            !     rhobar contains terms with qi/ni and ratio
            !     of qi/ni is assumed constant during loss of ni
            
            IF(qi(i,k).ge.qsmall.and.ni(i,k).gt.qsmall)THEN
               
               IF(prd.lt.0.)THEN
                  ni(i,k)=ni(i,k)+prd*ni(i,k)/qi(i,k)*dt
               END IF
               
               !     set minimum ni to avoid taking root of a negative number
               ni(i,k)=max(ni(i,k),qsmall)
               
               !     if iflag=1, then recalculate ai and ci assuming same deltastr
               !     and same rhobar
               !     this is needed for consistency between qi, ai, ci, etc.
               !     since growth rate prd must be scaled back
               
               IF(iflag.eq.1)THEN
                  alphstr=co/ao**(deltastr)
                  alphv=4./3.*pi*alphstr
                  betam=2.+deltastr   
                  
                  ani=((qi(i,k)*gammnu)/(rhobar*ni(i,k)*alphv*&
                       exp(gammln(nu+betam))))**(1./betam)
                  
                  cni=co*(ani/ao)**deltastr
               END IF           ! iflag = 1

               CALL R_CHECK(1,qi(i,k),ni(i,k),cni,ani,rni,rhobar,deltastr,&
                    betam,alphstr,alphv)
               
               ci(i,k)=nu*ni(i,k)*cni
               ai(i,k)=nu*ni(i,k)*ani
            END IF              ! q > qsmall


            if(snowflag.eq.2)then
               prds=max(prds,-qs(i,k)/dt)
               
               qs(i,k)=qs(i,k)+(prds)*dt
               as(i,k)=as(i,k)+(ards)*dt
               cs(i,k)=cs(i,k)+(crds)*dt
               
               IF(prds .gt. 0.0)THEN
                  snowdep(i,k) = prds
               ELSE
                  snowsub(i,k) = prds
               ENDIF

               !     if sublimation, reduce snow number,
               !     NOTE: this should not impact rhobars, since
               !     rhobars contains terms with qs/ns and ratio
               !     of qs/ns is assumed constant during loss of ns
               
               IF(qs(i,k).ge.qsmall.and.ns(i,k).gt.qsmall)THEN
               
                  IF(prds.lt.0.)THEN
                     ns(i,k)=ns(i,k)+prds*ns(i,k)/qs(i,k)*dt
                  END IF

                  !     set minimum ni to avoid taking root of a negative number
                  ns(i,k)=max(ns(i,k),qsmall)
                  
                  !     if iflag=1, then recalculate as and cs assuming same deltastrs
                  !     and same rhobars
                  !     this is needed for consistency between qs, as, cs, etc.
                  !     since growth rate prd must be scaled back
               
                  IF(iflag.eq.1)THEN
                     alphstr=co/ao**(deltastrs)
                     alphv=4./3.*pi*alphstr
                     betam=2.+deltastrs   
                     
                     ans=((qs(i,k)*exp(gammln(nus+betam)))/(rhobars*ns(i,k)*alphv*&
                          exp(gammln(nus+betam))))**(1./betam)
                     
                     cns=co*(ans/ao)**deltastrs
                  END IF           ! iflag = 1
                  
                  CALL R_CHECK(2,qs(i,k),ns(i,k),cns,ans,rns,rhobars,deltastrs,&
                       betam,alphstr,alphv)
               
                  cs(i,k)=nus*ns(i,k)*cns
                  as(i,k)=nus*ns(i,k)*ans
               END IF              ! q > qsmall
            end if
            

            !     calculate simplified aggregation for snow category 
            !     kjs 02/2015
            
            !     First, calulcate the aggregate-available ice,
            !     which is the amount of ice that would autoconvert to snow
            !     in traditional schemes. This is for r_i > 125 um.
            
            betam=2.+deltastr
            alphstr=co/ao**(deltastr)
            alphv=4./3.*pi*alphstr
            
            IF(snowflag.eq.1)THEN
               IF(qi(i,k).ge.1.e-8.and.qv(i,k)/qvi.gt.1.)THEN
                  nprci = 4./DCS/rhobar*(qv(i,k)-qvi)*rho(i,k)*ni(i,k)&
                       /rni*exp(-DCS/rni)*dv/abi ! #/kg/s
                  prci = 4./3.*pi*rhobar*DCS**3*nprci !kg/kg/s
                  !     Second, if the aggregate-available ice mass mixing ratio is
                  !     > 1.e-8 kg/kg, then calculate total aggregate number and mass
                  !     and add it to the snow category.  The aggeagates are assumed 
                  !     to be spheres and have a density of rhosn=100 kg/m3.
                  agg = 0.
                  nagg = 0.
                  IF(prci*dt.gt.1.e-8)THEN
                     dum2 = cons15*asn(k)
                     nagg = dum2*(rho(i,k)*prci*dt)**((2.+BS)/3.)*&
                          (nprci*dt*rho(i,k))**((4.-BS)/3.)/rho(i,k) !#/kg/s
                     agg = 4./3.*pi*rhosn*DCS**3*nagg !kg/kg/s
                     
                  END IF
                  !     Third, update qs, ns, qi, ni, ai, and ci due to aggregation. 
                  !     Assume rhobar remains the same.
                  !     Note: agg & nagg are negative numbers
                  
                  IF(nagg .ne. 0.0)THEN 
                     qs(i,k) = qs(i,k) - agg*dt
                     ns(i,k) = ns(i,k) - nagg*dt
                     
                     qi(i,k) = qi(i,k) + agg*dt
                     ni(i,k) = ni(i,k) + nagg*dt
                  END IF
                  
               END IF
            
               IF(qi(i,k).gt.qsmall.and.ni(i,k).gt.qsmall)THEN
                  ani = ((qi(i,k)*gammnu)/(rhobar*ni(i,k)*alphv*&
                       exp(gammln(nu+betam))))**(1./betam)
                  cni=co*(ani/ao)**deltastr
                  ci(i,k)=nu*ni(i,k)*cni
                  ai(i,k)=nu*ni(i,k)*ani  
               END IF
               ns(i,k)=max(ns(i,k),qsmall)
               IF(qs(i,k).lt.qsmall.or.ns(i,k).lt.qsmall)THEN
                  qs(i,k) = 0.0
                  ns(i,k) = 0.0
               END IF
            END IF!snowflag
            
!..................................................................
!     now update qi, ni, ai, and ci due to ice nucleation
!     for simplicity, assume that rhobar and deltastr are 
!     constant during nucleation (probably not the best assumption,
!     should look into future modifications....)
            
            qiold = qi(i,k)
            qi(i,k)=qi(i,k)+mnuccd*dt
            ni(i,k)=ni(i,k)+nnuccd*dt

!     set minimum ni to avoid division by zero
            ni(i,k)=max(ni(i,k),qsmall)

            IF(qi(i,k).ge.qsmall)THEN

               qi_frac = min(max(qiold/qi(i,k),0.0),1.0)
               rhobar = rhoi*(1.-qi_frac)+rhobar*(qi_frac)

               ani=((qi(i,k)*gammnu)/(rhobar*ni(i,k)*alphv*&
                   exp(gammln(nu+betam))))**(1./betam)
               cni=co*(ani/ao)**deltastr
               ci(i,k)=nu*ni(i,k)*cni
               ai(i,k)=nu*ni(i,k)*ani

               ai(i,k)=max(ai(i,k),qsmall)
               ci(i,k)=max(ci(i,k),qsmall)

!     final check on limits for deltastr, rhobar, and rn
               
!     get cni and ani from ci and ai

               ani=ai(i,k)/(nu*ni(i,k))
               cni=ci(i,k)/(nu*ni(i,k))

               CALL DSTR_CHECK(1,ni(i,k),ani,cni,deltastr,iaspect,sphrflag)
               CALL RHO_CHECK(1,deltastr,qi(i,k),ni(i,k),ani,cni,sphrflag,redden,&
                    betam,alphstr,alphv,rhobar)
               ci(i,k)=nu*ni(i,k)*cni
               ai(i,k)=nu*ni(i,k)*ani
         

            END IF ! qi > qsmall


!     add tendencies to temp, water vapor
!     if qi < qsmall, then zero out ni, ai, ci
            IF(qi(i,k).lt.qsmall)THEN
               qi(i,k)=0.
               ai(i,k)=0.
               ni(i,k)=0.
               ci(i,k)=0.
               rhobar=920.
            ELSE
               qv(i,k)=qv(i,k)-(prd+mnuccd)*dt
               temp=temp+((prd+mnuccd)*xxls/cpm*dt)
               icenuc(i,k) = mnuccd
            END IF

            IF(qs(i,k).lt.qsmall)THEN
               qs(i,k)=0.
               ns(i,k)=0.
            END IF

!     get fallspeed parameters
!     for now, use formulation from LH74, side planes
!     express v-D relationship in terms of 'a' length for simplicity
            
!     note: no air density correction factor applied yet!
            
            ni(i,k) = max(ni(i,k),qsmall)
            IF(qi(i,k).gt.qsmall.and.ni(i,k).gt.qsmall) THEN
               
!     updated ani, cni are calculated right before final checks above
               ni(i,k) = max(ni(i,k),qsmall)
               ani=ai(i,k)/(nu*ni(i,k))
               cni=ci(i,k)/(nu*ni(i,k))

               dum=exp(gammln(nu+deltastr+2.+bif))/&
               exp(gammln(nu+deltastr+2.))
               
               vtdumm=aif*dum*ani**bif
               dum=exp(gammln(nu+bif))/gammnu
               vtdumn=aif*dum*ani**bif
                   
!     limit fallspeed to 5 m/s

               
               vtrmi1(i,k)=max(min(vtbarbm,5.),0.)
               vtrni1(i,k)=max(min(vtbarb,5.),0.)
               vtrli1(i,k)=max(min(vtbarblen,5.),0.)

            ELSE
               
               vtrli1(i,k)=0.
               vtrmi1(i,k)=0.
               vtrni1(i,k)=0.
               
            END IF
!     get ice effective radius

!     for now just set to 10 micron

            effi1(i,k)=10.
            
!     get phi
            
            IF(qi(i,k).ge.qsmall.and.ani.gt.qsmall)THEN
               phii = cni/ani*exp(gammln(nu+deltastr-1.))/gammnu
            ELSE
               phii=1.
               ani=0.
               cni=0.
               rhobar=920.
            END IF
            IF(iaspect .eq. 1) phii = 0.27
            IF(sphrflag .eq. 1) rhobar =920.
            If(redden .eq. 1) rhoi = 500.
            !rhobar=500.
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     add tendencies to cloud, rain, and snow

            qc(i,k)=qc(i,k)+(-pra-prc)*dt
            qr(i,k)=qr(i,k)+(pra+prc+pre+pracs-psmlt)*dt
            nr(i,k)=nr(i,k)+(nprc1+nragg+nsubr)*dt
            IF(CNUM.ne.0)THEN
               nc(i,k)=nc(i,k)+(-npra-nprc)*dt
            END IF
            IF(snowflag .eq. 1)THEN
               qs(i,k) = qs(i,k)+(-pracs+psmlt+evpms+prds+eprds)*dt
               ns(i,k) = ns(i,k)+(nsagg+nsubs+nsmlts)*dt
            END IF
            qv(i,k)=qv(i,k)+(-pre-evpms-prds-eprds)*dt
            temp=temp+(xxlv*pre+(evpms+prds+eprds)*xxls+&
                   (psmlt-pracs)*xxlf)*dt/cpm

            rainevap(i,k) = -1.*pre
            if(snowflag.eq.1)then
               snowevap(i,k) = -1.*evpms
               snowmelt(i,k) = -1.*psmlt
               snowdep(i,k) = prds
               snowsub(i,k) = eprds
               snowaccr(i,k) = pracs
            end if
  
!     add liquid, do saturation adjustment to get updated qc
!     get updated thermodynamics (updated after ice microphysics)
            
!     get temp from 0theta (in K)
            
            qvs=0.622*polysvp1(temp,0)/(press-polysvp1(temp,0))
            dum=qv(i,k)-qvs
            pcc=dum/(1.+xxlv**2*qvs/(cpm*rv*temp**2))/dt
            IF(pcc*dt+qc(i,k).lt.0.)THEN
               pcc=-qc(i,k)/dt
            END IF
            !condrate(i,k)=pcc
!     update variables due to condensation/evaporation

            temp=temp+pcc*xxlv/cpm*dt
            qv(i,k)=qv(i,k)-pcc*dt
            qc(i,k)=qc(i,k)+pcc*dt


            !print*,'here',i,k,pcc
            IF(pcc .gt. 0.0)THEN
               cloudcond(i,k) = pcc
            ELSE
               cloudevap(i,k) = pcc
            END IF
!     for now, assume droplet fallspeed is zero           
            vtrmc(i,k)=0.

            !temp3(i,k)=temp
            t(i,k) = temp
           

!     add melting, melt all ice within one time step above freezing
            IF(t(i,k).ge.273.15)THEN
               IF(qi(i,k).ge.qsmall)THEN
                  qr(i,k)=qr(i,k)+qi(i,k)
                  nr(i,k)=nr(i,k)+ni(i,k)
                  t(i,k)=t(i,k)-qi(i,k)*xxlf/cpm
                  icemelt(i,k) = qi(i,k)/dt
                  qi(i,k)=0.
                  ni(i,k)=0.
                  ai(i,k)=0.
                  ci(i,k)=0.
               END IF

               IF(snowflag.eq.2)then
                  IF(qs(i,k).ge.qsmall)THEN
                     qr(i,k)=qr(i,k)+qs(i,k)
                     nr(i,k)=nr(i,k)+ns(i,k)
                     t(i,k)=t(i,k)-qs(i,k)*xxlf/cpm
                     snowmelt(i,k) = qs(i,k)/dt
                     qs(i,k)=0.
                     ns(i,k)=0.
                  END IF
               END IF
            END IF
!     homogeneous freezing, freeze andd cloud and rain water within one time-step below -40
            IF(homofreeze .eq. 1.and.t(i,k).le.233.15.and.&
               real(itimestep)*dt.gt.ice_start_time)THEN
               IF(qr(i,k).ge.qsmall)THEN

                  qi(i,k)=qi(i,k)+qr(i,k)
                  ni(i,k)=ni(i,k)+nr(i,k)
                  t(i,k)=t(i,k)+qr(i,k)*(xxls-xxlv)/cpm
!                  nnew = qr(i,k)/((4./3.)*pi*rhoi*(200.e-6)**3)
                  ai(i,k) = ai(i,k) + nu*nr(i,k)*(qr(i,k)*3./(4.*pi*&
                  rhoi*nr(i,k))*gammnu/exp(gammln(nu+3)))**(1./3.)
                  ci(i,k) = ci(i,k) + nu*nr(i,k)*(qr(i,k)*3./(4.*pi*&
                  rhoi*nr(i,k))*gammnu/exp(gammln(nu+3)))**(1./3.)
                  qr(i,k)=0.
                  nr(i,k)=0.
                  rainfrz(i,k) = qr(i,k)/dt 
               END IF
               IF(qc(i,k).ge.qsmall)THEN

                  qi(i,k)=qi(i,k)+qc(i,k)
                  t(i,k)=t(i,k)+qc(i,k)*(xxls-xxlv)/cpm
                  nnew = qc(i,k)/((4./3.)*pi*rhoi*(20.e-6)**3)
                  ai(i,k) = ai(i,k) + nu*nnew*(qc(i,k)*3./(4.*pi*&
                  rhoi*nnew)*gammnu/exp(gammln(nu+3)))**(1./3.)
                  ci(i,k) = ci(i,k) + nu*nnew*(qc(i,k)*3./(4.*pi*&
                  rhoi*nnew)*gammnu/exp(gammln(nu+3)))**(1./3.)
                  ni(i,k) = ni(i,k) + nnew
                  qc(i,k)=0.
                  cloudfrz(i,k) = qc(i,k)/dt
               END IF
            END IF
 
            LTRUE = 1
            
200         CONTINUE

            IF(qi(i,k) .gt. qsmall)THEN
               ani = ai(i,k)/(ni(i,k)*nu)
               cni = ci(i,k)/(ni(i,k)*nu)

               ai(i,k) = ani**2*cni*nu*ni(i,k)
               ci(i,k) = cni**2*ani*nu*ni(i,k)
            END IF
            rhoice(i,k)=rhobar
         END DO
               IF(LTRUE.eq.0)GOTO 400
!     SEDIMENTATION --------------------------------------------------------------
        IF(SEDON .eq. 1)THEN
  
         nstep = 1
         
         DO k = kts,kte
           
            IF(qi(i,k).ge.qsmall)THEN
               fi(k) = vtrmi1(i,k)
               fni(k) = vtrni1(i,k)
               fai(k) = vtrmi1(i,k)
               fci(k) = vtrmi1(i,k)
            ELSE
               fi(k) = 0.0
               fni(k) = 0.0
               fai(k) = 0.0
               fci(k) = 0.0
            END IF
            !sedm_i(i,k) = fi(k)

!     calculate rain fallspeed
            nr(i,k) = max(nr(i,k),0.0)
            IF(qr(i,k).ge.qsmall)THEN!.and.nr(i,k).gt.0.0)THEN
               lamr = (pi*rhow*nr(i,k)/qr(i,k))**0.33333
               lamr = max(lamr,lamminr)
               lamr = min(lamr,lammaxr)
               fr(k) = arn(k)*cons4/lamr**br
               fnr(k) = arn(k)*cons6/lamr**br
               fr(k) = min(fr(k),9.1*(rhosu/rho(i,k))**0.54)
               fnr(k) = min(fnr(k),9.1*(rhosu/rho(i,k))**0.54)
            ELSE
               fr(k) = 0.
               fnr(k) = 0.
            END IF
            !sedm_l(i,k) = fr(k)
!     hm add cloud water sedimentation
            nc(i,k) = max(nc(i,k),0.)
            IF(qc(i,k).ge.qsmall)THEN!.and.nc(i,k).gt.0.0)THEN
               DUM = p(i,k)/(287.15*t(i,k))
               PGAM = 0.0005714*(nc(i,k)/1.E6*DUM)+0.2741
               PGAM = 1./(PGAM**2)-1.
               PGAM = MAX(PGAM,2.)
               PGAM = MIN(PGAM,10.)

               lamc = (CONS26*nc(i,k)*GAMMA(PGAM+4.)/&
               (qc(i,k)*GAMMA(PGAM+1.)))**0.33333
               LAMMIN = (PGAM+1.)/60.E-6
               LAMMAX = (PGAM+1.)/1.E-6
               lamc=MAX(lamc,LAMMIN)
               lamc=MIN(lamc,LAMMAX)

               fc(k) = acn(k)*GAMMA(4.+BC+PGAM)/(lamc**BC*GAMMA(PGAM+4.))
               fnc(k) = acn(k)*gamma(1.+BC+pgam)/(lamc**BC*gamma(pgam+1.))
             ELSE
               fc(k) = 0.
               fnc(k) = 0.
             END IF
            
             !     calculate snow/aggregate sedimentation
             IF(snowflag.eq.1)THEN
                IF(qs(i,k).ge.qsmall)THEN
                   dum = (cons1*ns(i,k)/qs(i,k))**(1./DS)
                   fs(k) = asn(k)*cons3/dum**BS
                   fns(k) = asn(k)*cons5/dum**BS
                   fs(k) = min(fs(k),1.2*(rhosu/rho(i,k))**0.54)
                   fns(k) = min(fns(k),1.2*(rhosu/rho(i,k))**0.54)
                ELSE
                   fs(k) = 0.
                   fns(k) = 0.
                END IF
             END IF!snowflag


!     calculate number of split time steps

            rgvm = max(fi(k),fni(k),fai(k),fci(k),fr(k),fc(k),fs(k),fns(k))
            nstep = max(int(rgvm*dt/dzq(k)+1.),nstep)
            
            
!     multiply variables by rho for right units
            dumi(k) = qi(i,k)*rho(i,k) !kgm^-3
            dumni(k) = ni(i,k)*rho(i,k) !#m^-3
            dumai(k) = ai(i,k)*rho(i,k) !unitless (ai is volume now)
            dumci(k) = ci(i,k)*rho(i,k) !unitless (ci is volume now)
            dumr(k) = qr(i,k)*rho(i,k) !kgm^-3
            dumc(k) = qc(i,k)*rho(i,k) !kgm^-3
            dumnr(k) = nr(i,k)*rho(i,k) !kgm^-3
            dumnc(k) = nc(i,k)*rho(i,k) !kgm^-3
            dums(k) = qs(i,k)*rho(i,k)
            dumns(k) = ns(i,k)*rho(i,k)
         END DO

         DO n = 1, nstep
            DO k = kts,kte
               fallouti(k) = fi(k)*dumi(k) !kgm^-2 s^-1
               falloutni(k) = fni(k)*dumni(k) !#m^-2 s^-1
               falloutai(k) = fai(k)*dumai(k) !m s^-1
               falloutci(k) = fci(k)*dumci(k) !m s^-1
               falloutr(k) = fr(k)*dumr(k) !kgm^-2 s^-1
               falloutc(k) = fc(k)*dumc(k) !kgm^-2 s^-1
               falloutnr(k) = fnr(k)*dumnr(k) !kgm^-2 s^-1
               falloutnc(k) = fnc(k)*dumnc(k) !kgm^-2 s^-1
               fallouts(k) = fs(k)*dums(k)
               falloutns(k) = fns(k)*dumns(k)
            END DO
            
!     top of model
            
            k = kte
            falltndi = fallouti(k)/dzq(k) !kgm^-3 s^-1
            falltndni = falloutni(k)/dzq(k) !#m^-3 s^-1
            falltndai = falloutai(k)/dzq(k) !s^-1
            falltndci = falloutci(k)/dzq(k) !s^-1
            falltndr = falloutr(k)/dzq(k) !kgm^-3 s^-1
            falltndc = falloutc(k)/dzq(k) !kgm^-3 s^-1
            falltndnr = falloutnr(k)/dzq(k) !kgm^-3 s^-1
            falltndnc = falloutnc(k)/dzq(k) !kgm^-3 s^-1
            falltnds = fallouts(k)/dzq(k)
            falltndns = falloutns(k)/dzq(k)
            
!     sedimentation tendencies should be renewed every time step
!     so be sure to initialize to zero at beginning of code (kjs)
            
            qisten(k) = qisten(k) - falltndi/nstep/rho(i,k) !s^-1
            nisten(k) = nisten(k) - falltndni/nstep/rho(i,k) !#kg^-1 s^-1
            aisten(k) = aisten(k) - falltndai/nstep/rho(i,k) !m^3 kg^-1 s^-1
            cisten(k) = cisten(k) - falltndci/nstep/rho(i,k) !m^3 kg^-1 s^-1
            qrsten(k) = qrsten(k) - falltndr/nstep/rho(i,k) !s^-1
            qcsten(k) = qcsten(k) - falltndc/nstep/rho(i,k) !s^-1
            nrsten(k) = nrsten(k) - falltndnr/nstep/rho(i,k) !s^-1
            ncsten(k) = ncsten(k) - falltndnc/nstep/rho(i,k) !s^-1
            qssten(k) = qssten(k) - falltnds/nstep/rho(i,k)
            nssten(k) = nssten(k) - falltndns/nstep/rho(i,k)
            
            dumi(k) = dumi(k) - falltndi*dt/nstep             
            dumni(k) = dumni(k) - falltndni*dt/nstep
            dumai(k) = dumai(k) - falltndai*dt/nstep
            dumci(k) = dumci(k) - falltndci*dt/nstep
            dumr(k) = dumr(k) - falltndr*dt/nstep
            dumc(k) = dumc(k) - falltndc*dt/nstep
            dumnr(k) = dumnr(k) - falltndnr*dt/nstep
            dumnc(k) = dumnc(k) - falltndnc*dt/nstep
            dums(k) = dums(k) - falltnds*dt/nstep
            dumns(k) = dumns(k) - falltndns*dt/nstep
            
            
            DO k = kte-1,kts,-1
               
               falltndi = (fallouti(k+1) - fallouti(k))/dzq(k)
               falltndni = (falloutni(k+1) - falloutni(k))/dzq(k)
               falltndai = (falloutai(k+1) - falloutai(k))/dzq(k)
               falltndci = (falloutci(k+1) - falloutci(k))/dzq(k)
               falltndr = (falloutr(k+1) - falloutr(k))/dzq(k)
               falltndc = (falloutc(k+1) - falloutc(k))/dzq(k)
               falltndnr = (falloutnr(k+1) - falloutnr(k))/dzq(k)
               falltndnc = (falloutnc(k+1) - falloutnc(k))/dzq(k)
               falltnds = (fallouts(k+1) - fallouts(k))/dzq(k)
               falltndns = (falloutns(k+1) - falloutns(k))/dzq(k)
               
               qisten(k) = qisten(k) + falltndi/nstep/rho(i,k)
               nisten(k) = nisten(k) + falltndni/nstep/rho(i,k)
               aisten(k) = aisten(k) + falltndai/nstep/rho(i,k)
               cisten(k) = cisten(k) + falltndci/nstep/rho(i,k)
               qrsten(k) = qrsten(k) + falltndr/nstep/rho(i,k)
               qcsten(k) = qcsten(k) + falltndc/nstep/rho(i,k)
               nrsten(k) = nrsten(k) + falltndnr/nstep/rho(i,k)
               ncsten(k) = ncsten(k) + falltndnc/nstep/rho(i,k)
               qssten(k) = qssten(k) + falltnds/nstep/rho(i,k)
               nssten(k) = nssten(k) + falltndns/nstep/rho(i,k)

               dumi(k) = dumi(k) + falltndi*dt/nstep
               dumni(k) = dumni(k) + falltndni*dt/nstep
               dumai(k) = dumai(k) + falltndai*dt/nstep
               dumci(k) = dumci(k) + falltndci*dt/nstep
               dumr(k) = dumr(k) + falltndr*dt/nstep
               dumc(k) = dumc(k) + falltndc*dt/nstep
               dumnr(k) = dumnr(k) + falltndnr*dt/nstep
               dumnc(k) = dumnc(k) + falltndnc*dt/nstep
               dums(k) = dums(k) + falltnds*dt/nstep
               dumns(k) = dumns(k) + falltndns*dt/nstep

               !falloutL(i,k) = falloutr(k)+falloutc(k) 
               !falloutICE(i,k) = fallouti(k)
  
               qiloss(k)=qiloss(k)+fallouti(k)
               qrloss(k)=qrloss(k)+falloutr(k)+falloutc(k)
            END DO
            
!     get precipitation and snowfall accumulation during time step
!     NOTE: factor of 1000 converts m to mm, but division by density cancels this
            
            precrt = precrt + (fallouti(kts)+falloutr(kts)+&
            falloutc(kts)+fallouts(kts))*dt/nstep !kgm^-2
            snowrt = snowrt + (fallouti(kts)+fallouts(kts))*dt/nstep

         END DO                 !end nstep loop
         END IF  !SEDON       

!     add on sedimenation tendencies for mixing ratio to rest of tendencies
         DO k=kts,kte
!     add new tendencies to mixing ratios

            qi(i,k) = qi(i,k) + qisten(k)*dt !kg kg^-1
            ni(i,k) = ni(i,k) + nisten(k)*dt !#kg^-1
            ai(i,k) = ai(i,k) + aisten(k)*dt !m kg^-1
            ci(i,k) = ci(i,k) + cisten(k)*dt !m kg^-1
            qr(i,k) = qr(i,k) + qrsten(k)*dt !kg kg^-1
            qc(i,k) = qc(i,k) + qcsten(k)*dt !kg kg^-1
            nr(i,k) = nr(i,k) + nrsten(k)*dt !kg kg^-1
            nc(i,k) = nc(i,k) + ncsten(k)*dt !kg kg^-1
            qs(i,k) = qs(i,k) + qssten(k)*dt
            ns(i,k) = ns(i,k) + nssten(k)*dt

            IF(qi(i,k) .lt. qsmall)THEN
               qi(i,k) = 0.
               ni(i,k) = 0.
               ai(i,k) = 0.
               ci(i,k) = 0.
               qiloss(k) = 0.
            END IF
           
            IF(qr(i,k) .lt. qsmall)THEN
               qr(i,k) = 0.
               nr(i,k) = 0.
               qrloss(k) = 0.
            END IF
            IF(qc(i,k) .lt. qsmall)THEN
               qc(i,k) = 0.
               nc(i,k) = 0.
               qcloss(k) = 0.
            END IF
            IF(qs(i,k).lt.qsmall)THEN
               qs(i,k) = 0.0
               ns(i,k) = 0.0
            END IF
            
!     recalculate cloud, rain, and snow distributions
            IF(qr(i,k).ge.qsmall)THEN
               lamr = (pi*rhow*nr(i,k)/qr(i,k))**0.33333
               
               IF(lamr.lt.lamminr)THEN
                  lamr = lamminr
                  n0rr = lamr**4*qr(i,k)/(pi*rhow)
                  nr(i,k) = n0rr/lamr
               ELSE IF(lamr.gt.lammaxr)THEN
                  lamr = lammaxr
                  n0rr = lamr**4*qr(i,k)/(pi*rhow)
                  nr(i,k) = n0rr/lamr
               END IF
            END IF

            IF(qc(i,k) .ge. qsmall)THEN
               dum = p(i,k)/(287.15*t(i,k))
               pgam = 0.0005714*(nc(i,k)/1.e6*dum)+0.2714
               pgam = 1./(pgam**2)-1.
               pgam = max(pgam,2.)
               pgam = min(pgam,10.)

               lamc = (cons26*nc(i,k)*gamma(pgam+4.)/&
               (qc(i,k)*gamma(pgam+1.)))**0.3333

               lammin = (pgam+1.)/60.e-6
               lammax = (pgam+1.)/1.e-6

               IF(lamc.lt.lammin)THEN
                  lamc = lammin
                  nc(i,k) = exp(3.*log(lamc)+log(qc(i,k))+&
                  log(gamma(pgam+1.))-log(gamma(pgam+4.)))/cons26
               ELSE IF(lamc.gt.lammax)THEN
                  lamc = lammax
                  nc(i,k) = exp(3.*log(lamc)+log(qc(i,k))+&
                  log(gamma(pgam+1.))-log(gamma(pgam+4.)))/cons26
               END IF
            END IF
            
            IF(snowflag.eq.1)THEN
               IF(qs(i,k).ge.qsmall)THEN
                  lams = (cons1*ns(i,k)/qs(i,k))**(1./DS)
                  
                  IF(lams.lt.lammins)THEN
                     lams = lammins
                     n0s = lams**4*qs(i,k)/cons1
                     ns(i,k) = n0s/lams
                  ELSE IF(lams.gt.lammaxs)THEN
                     lams = lammaxs
                     n0s = lams**4*qs(i,k)/cons1
                     ns(i,k) = n0s/lams
                  END IF
               END IF
            END IF!snowflag

         END DO

 400     CONTINUE
      END DO     
               !end i loop
      DO k=kts,kte
         DO i=its,ite
            !temp4(i,k)=t(i,k)
            !qtot2(i,k)=(qv(i,k)+qc(i,k)+qi(i,k)+qr(i,k)+qs(i,k))!*rho(i,k)
            !qt_sed=qt_sed+qtot2(i,k)*rho(i,k)*dzq(k)
            thetal(k)=theta/t(i,k)*(1.-(xxlv/cp*(qc(i,k)+qi(i,k))))
         END DO
         IF(k.gt.kts) dthl(k)=thetal(k)-thetal(k-1)
      END DO



      IF(resupply .EQ. 1)THEN   
        level2=0
        topflag=1
        dthlmax=-10.
        DO k=kte-4,kts,-1
           IF(topflag.eq.1)THEN
              IF(dthlmax.lt.dthl(k))THEN
                 dthlmax=dthl(k)
              ELSE  
                 level2=k    !mixed-layer top
                 topflag=0
              END IF
           END IF
        END DO
        
        level=0
        bottomflag=1
        dthlmax=-10.
        qsum2=0.0
        DO k=kts,kte
           IF(bottomflag .eq. 1)THEN
              IF(dthlmax.lt.dthl(k))THEN
                 dthlmax=dthl(k)
              ELSE
                 level=k        !mixed-layer bottom
                 bottomflag=0
              END IF
           END IF
!           IF(k.ne.0.and.k.eq.level)THEN
!              qsum2=-qiloss(k)-qrloss(k)
!              qsum2=qsum2*dt
!           END IF
        END DO
        qsum2=(qt_adv-qt_sed)

        IF(level.lt.level2)THEN        
           DO k=level,level2
              DO i=its,ite
                 qv(i,k)=qv(i,k)+qsum2/rho(i,k)/dzq(k)&
                 /float(ite-its+1)/float(level2-level)
                 t(i,k)=t(i,k)-qsum2*xxlv/cpm/float(level2-level)/ &
                 float(ite-its+1)
              END DO
           END DO
         END IF
      END IF


      
!      print*,'Time = ',itimestep*dt, ' seconds'
      RETURN
      
    END SUBROUTINE SULIAHARRINGTON_MICRO


!************************************************************************
! SUBROUTINE EVOLVE
!   Purpose: This routine solves the equations for the growth of 
!            non-spherical ice following the Chen and Lamb (1994)
!            method as per the document of Harrington and Sulia (2010)
!
!   Variables: (NOTE: All lengths are characteristic lengths!!!)
!           ani = initial a-axis length at time t
!           nu = shape of the distribution
!           ni = number density of ice crystals
!           sui = ice supersaturation
!           temp = temperaure
!           press = pressure
!           igr = inherent growth ratio
!           deltt = time-step 
!           iwc = final ice water content at time t+deltt
!           cf = final c-axis length at time t+deltt
!           iwci = initial ice water content at time t
!           phii = initial aspect ratio at time t
!           phif = final aspect ratio at time t+deltt
!           cni = initial c-axis length at time t
!           rni = initial equivalent volume radius at time t
!           rnf = final equivalent volume radius at time t+deltt
!           anf = final a-axis length at time t+deltt
!           deltastr = historical factor relating the a and c axes
!
!------------------------------------------------------------------------

      SUBROUTINE EVOLVE(iflag,ani,ni,sui,sup,qvv,temp,press,igr,deltt,iwc,cf,&
          iwci,phii,phif,cni,rni,rnf,anf,deltastr,mu,rhoavg,vtbarb,&
          vtbarbm,alphstr,vtbarblen,rhoa,i,k,iaspect,&
          ipart,sphrflag,redden,itimestep)
!***********************************************************************

      IMPLICIT NONE

      INTEGER ipart,ivent,i,k,iaspect,redden,iflag
      INTEGER sphrflag,itimestep,MICRO_ON
      REAL ani,ni,temp,press,qvv, nn
      REAL igr,gi,capgam,dmdt,deltt,qe
      REAL dlndt,mu
      REAL l,lmean,fs,gammnu1,gammnubet,anf,betam
      REAL iwc,cf,alphstr,deltastr,iwci,phii,phif
      REAL cni,gammnu2delt,gammnu1delt,r,rf,vf
      REAL gamrats,alphan,rni,rnf,alphanr,rhoavg,rhodep,rhosub
      REAL drho,igr_den,vi,betavol,Vmin,videp,fallcheck
      REAL fv,dvs,kts,grav,rhoa,am,bm,bv1,bv2
      REAL bt1,bt2,bx,gv,gt,fh,ntherm,xvent,nre,nsch,npr,xm
      REAL bl,ba,aa,al,etaa,celsius,wghtv1
      REAL igrvent,wghtv2,wghtv3,vtbarb,vtbarbm,vtbarblen,vt
      REAL cv, bv, bm1
      REAL alpham,alpha_a,beta_a,alpha_v,beta_v,cap
      REAL alphams, alpha_as, alphap, alphamsp, betamp
      REAL*8 sui,sup,xn

      MICRO_ON = 1
      IF(MICRO_ON .EQ. 1)THEN
      bv = 0.5
      cv = 12.*2.**bv
      bm1 = 3.0

      ivent=1
      grav = 9.81
      celsius = temp-273.15
      fv = 1.0
      fh = 1.0
      gi = findgtp(temp,press,fv,fh)
      dvs = 0.211*(temp/273.15)**1.94 * &
          (1013.25*100./press)*1.0/100.0**2.0 ! vapor diffusivity
      kts = 2.3823e-2 +7.1177e-5*(temp-273.15) ! thermal conductivity
      
      l=ani
      drho = ((polysvp1(temp,0)-polysvp1(temp,1))/(Rv*temp))*1000.0
      IGR_DEN = IGR

      betam = 2. + deltastr

      if(iflag.eq.1) nn = nu
      if(iflag.eq.2) nn = nus
      
      gammnu=exp(gammln(nn))
      gammnu1=exp(gammln(nn+1.0))
      gammnubet=exp(gammln(nn+betam))

      capgam = capacitance_gamma(l,deltastr,nn,alphstr)
      lmean = l*nn

      cap = 0.
      alpham = 0.
      alphams = 0.
      alpha_as = 0.
      alphap = 0.
      alphamsp = 0.
      betamp = 0.
      alpha_a = 0.
      beta_a = 0.
      alpha_v = 0.
      beta_v = 0.

      dmdt = 4.*pi*gi*sui*ni*capgam

        !     get alpham and betam from deltastr
      dlndt = dmdt/(ni*lenconv(ani,deltastr,nn,alpham,&
      betam,alphstr,rhoavg))
      gammnubet = exp(gammln(nn+betam))
      gammnu2delt = exp(gammln(nn+2.+deltastr))
      gammnu1delt = exp(gammln(nn+deltastr-1.))
      gamrats = exp(gammln(nn+betam/3.))/gammnu
      
      phii = cni/ani*gammnu1delt/gammnu ! diagnose initial aspect ratio
      IF(iaspect .eq. 1) phii = 0.27
      
      r = alphstr**(1./3.) * (l)**(betam/3.)
!      fs = capgam/(r*gamrats)
      fs = capgam/(rni*gamrats)
      alphan = cni*ani**(-igr)
      alphanr = ani/rni**(3./(2.+igr))

      if(phii.lt.1.0)then
         bl = 1.
         al = 2.
         aa = pi
         ba = 2.
         qe = (rhoavg/920.)**(2./3.)
      else if(phii .gt. 1.0)then
         al = 2.*alphstr**(1./3.)
         bl = (deltastr + 2.)/3.
         aa = pi*alphstr
         ba = deltastr + 1.
         qe = 1.0
      else if(phii .eq. 1.0)then
         bl = 1.
         al = 2.
         aa = pi
         ba = 2.
         qe = 1.0
      endif
      
      etaa = (1.718 + 0.0049*celsius - 1.2e-5*celsius**2)*1.e-4 &
      *100./1000.               ! dynamic viscisoty
      nsch = etaa/(rhoa*dvs)    ! Schmidt number
      npr = etaa/(rhoa*kts)     ! Prandlt number

!      xn =  2./rhoavg*(rhoavg-rhoa)*grav*rhoa/etaa**2 & 
!     *(4./3.*pi*rhoavg)*alphstr*al**2/(aa*qe) * qe**(3./4.)
      xn =  2./rhoavg*(rhoavg-rhoa)*grav*rhoa/etaa**2 & 
      *(4./3.*pi*rhoavg)*10.**(alog10(alphstr)+2.* &
      alog10(al)-alog10(aa))/qe * qe**(3./4.)
      
      bx = deltastr+2.+2.*bl-ba
      

      wghtv1 = exp(gammln(nn+deltastr+2.+2.*bl-ba))/gammnu
      xm = xn*ani**bx * wghtv1  ! average Best number
      
      if(xm.le.10.)then
         am = 0.04394
         bm = 0.970
      else if(xm.gt.10.0.and.xm.le.585.)then
         am = 0.06049
         bm = 0.831
      else if(xm.gt.585.and.xm.le.1.56e5)then
         am = 0.2072
         bm = 0.638
      else if(xm.gt.1.56e5.and.xm.lt.1.e8)then
         am = 1.0865
         bm = 0.499
      else
!print*,' out of range!!!'
         am = 1.0865
         bm = 0.499

      endif
      
      Nre = am*xm**bm
      
      wghtv2 = exp(gammln(nn+bx*bm-1.))/gammnu
      vtbarb = etaa/rhoa*0.5 * am*(xn)**bm*ani**(bx*bm-1.) * wghtv2
      wghtv3 = exp(gammln(nn+bx*bm-1.+2.+deltastr))&
      /exp(gammln(nn+2.+deltastr))
      vtbarbm = etaa/rhoa*0.5 * am*(xn)**bm*ani**(bx*bm-1.) * wghtv3

 
! add length
      if(phii.lt.1.0)then
         wghtv3 = exp(gammln(nn+bx*bm-1.+1)) &
         /exp(gammln(nn+1.))
         vtbarblen = etaa/rhoa*0.5 * am*(xn)**bm*ani**(bx*bm-1.) &
         * wghtv3
      else if(phii.gt.1.0)then
         wghtv3 = exp(gammln(nn+bx*bm-1.+deltastr))&
         /exp(gammln(nn+deltastr))
         vtbarblen = etaa/rhoa*0.5 * am*(xn)**bm*ani**(bx*bm-1.) &
         * wghtv3
      else if (phii.eq.1.0)then
         wghtv3 = exp(gammln(nn+bx*bm-1.+1))&
         /exp(gammln(nn+1.))
         vtbarblen = etaa/rhoa*0.5 * am*(xn)**bm*ani**(bx*bm-1.) &
         * wghtv3
      endif
      
      

         
      xvent = nsch**(1./3.)*SQRT(nre)
      ntherm = SQRT(nre)*npr**(1./3.)

      if(xvent.le.1.0)then
         bv1 = 1.0
         bv2 = 0.14
         gv = 2.
      else
         bv1 = 0.86
         bv2 = 0.28
         gv = 1.
      endif

      if(ntherm.lt.1.4)then
         bt1 = 1.0
         bt2 = 0.108
         gt = 2.0
         !print*,'ntherm = ',ntherm
      else
         bt1 = 0.78
         bt2 = 0.308
         gt = 1.0
      endif

      if(ivent.eq.1)then
         fv = bv1 + bv2*xvent**gv
!         fva = bv1 + bv2*xvent**gv*(ani/rni)**(gv/2.)
!         fva = bv1 + bv2*xvent**gv*(alphstr**(-1./3.)*
!     +        ln**(1.-(deltastr+2.)/3.) * 
!     +        exp(gammln(nn + 1. -(deltastr+2.)/3.))/gammnu)**(gv/2.)
!         fvc = bv1 + bv2*xvent**gv*(cni/rni)**(gv/2.)
!         fvc = bv1 + bv2*xvent**gv*(alphstr**(2./3.)*
!     +        ln**(deltastr-(deltastr+2.)/3.) * 
!     +        exp(gammln(nn + deltastr -(deltastr+2.)/3.))/gammnu)
!     +        **(gv/2.)
         fh = bt1 + bt2*ntherm**gt

         igrvent = igr
      else
         fv = 1.0
         fh = 1.0
         igrvent = igr
      endif

      gi = findgtp(temp,press,fv,fh)
      

!      rhodep = 920.*exp((-3.0*max(drho-0.05,0.0))/IGR_DEN) ! get density of ice deposited
 
!      rhodep = 920.                 
                                              ! during growth
      vt = mu/rhoa*0.5*Nre/ani
      fallcheck = sqrt(dvs*pi*2*cni/(vt*nn))
      rhosub = rhoi*((1.-sui/qvv)+igr*(sui/qvv))
      IF(igr .lt. 1.0)THEN
         IF(sup .ge. 0.0)THEN
            IF(ani .gt. fallcheck)THEN
               rhodep = igr*rhoi
            ELSE
               rhodep = rhoi
            END IF
         ELSE
            IF(ani .gt. fallcheck)THEN
               rhodep = rhosub
            ELSE
               rhodep = rhoi
            END IF
         END IF
      ELSE IF(igr.gt.1.0)THEN 
         IF(sup .ge. 0.0)THEN
            rhodep = rhoi/igr
         ELSE
            rhodep = rhosub
         END IF
      ELSE
         rhodep = 920.
      END IF
      if(sui.lt.0.0)then        ! OR get density of ice removed 
         rhodep = rhoavg        ! during sublimation using polynomial
         videp = 4./3.*pi*rni**3 & 
         *exp(gammln(nn+deltastr+2.))/gammnu
         Vmin = 4./3.*pi*ao**2*co
         if(Vmin.lt.videp)then
            betavol = alog(rhoi/rhoavg)*1./(alog(Vmin/videp))
            rhodep = rhoavg*(1.+betavol)
         else
            rhodep = rhoavg 
         endif
      endif


      IF(sphrflag.eq.1) rhodep=920.
      If(redden .eq. 1) rhodep = 500.
      !rhodep=500.
!     rf = ((rni*gamrats)**2 + 2.*gi*sui*fs/rhodep &
!           *gammnu/gammnubet*gamrats**3*deltt)**(0.5) ! rmean after deltat

      rf = (max(real((rni*gamrats)**2 + 2.*gi*sui*fs/rhodep &
           *gammnu/gammnubet*gamrats**3*deltt),1.e-20))**(0.5) ! rmean after deltat

      rnf = rf/gamrats           ! convert to rn (characteristic r)
      anf = alphanr*rnf**(3./(2.+igr)) ! characteristic a-axis after deltat
      vi = 4./3.*pi*rni**3 * exp(gammln(nn+deltastr+2.))/gammnu ! mean initial volume
      vf = 4./3.*pi*rnf**3 * exp(gammln(nn+deltastr+2.))/gammnu ! mean volume after deltat
      rhoavg = rhoavg*(vi/vf) + rhodep*(1.-vi/vf) ! new average ice density
      rhoavg = max(min(rhoavg,920.),50.)

      IF(sphrflag.eq.1)rhoavg=920.
      !rhoavg=500.
      iwc = ni*rhoavg*vf        ! IWC after deltat

      if(igr.ne.1.0.and.(log(anf)-log(co)).gt.0.01 .and.&
      (3.*log(rnf)-2.*log(anf)-log(ao)).gt.0.001)THEN
         deltastr = (3.*log(rnf)-2.*log(anf)-log(co)) &
         /(log(anf)-log(ao))    ! diagnose new deltastar
      else
         deltastr = 1.
      endif
      IF(iaspect .eq. 1) deltastr = 0.8
      deltastr = min(max(deltastr,0.55),1.5)
  
      if(deltastr.ge.1.0)then   ! if columns, get c from c-a relation
         cf = co*(anf/ao)**deltastr
      endif

      phif = phii*(rnf**3/r**3)**((igr-1.)/(igr+2.)) ! if plates, get c from aspect ratio after deltat
      IF(iaspect .eq.1) phif = 0.27

!      rnf = rf
      if(deltastr.lt.1.0)then
         cf = phif*anf*gammnu/exp(gammln(nn+deltastr-1.))
      endif

     
      END IF !MICRO_ON

 
      return
    END SUBROUTINE EVOLVE

    !********************************************************
    !THIS SUBROUTINE COMPUTES DELTASTR, ICE DENSITY, AND ICE
    !EQUIVALENT VOL RADIUS AND SUBSEQUENTLY PERFORMS LIMIT
    !CHECKS.
    !********************************************************
    SUBROUTINE ICE_CHECKS(iflag,ni,qi,ani,cni,rni,deltastr,rhobar,iaspect,&
         sphrflag,redden,betam,alphstr,alphv)
                 
      IMPLICIT NONE
      INTEGER iaspect, sphrflag, redden,iflag
      REAL qi, ni, ani, cni, rni, deltastr, rhobar
      REAL betam,alphstr,alphv
      
      CALL DSTR_CHECK(iflag,ni,ani,cni,deltastr,iaspect,sphrflag)
      CALL RHO_CHECK(iflag,deltastr,qi,ni,ani,cni,sphrflag,redden,betam,alphstr,alphv,rhobar)
      CALL R_CHECK(iflag,qi,ni,cni,ani,rni,rhobar,deltastr,betam,alphstr,alphv)
                 

    END SUBROUTINE ICE_CHECKS

    SUBROUTINE DSTR_CHECK(iflag,ni,ani,cni,deltastr,iaspect,sphrflag)
      
      IMPLICIT NONE
      INTEGER :: iaspect,sphrflag,iflag
      REAL ::  ani, cni, deltastr, voltmp, ai, ci, ni, n, gn

      !     get deltastr from cni and ani
      !     deltastr = 1 for ice particles pre-diagnosed as spheres

      IF(iflag .eq. 1)THEN
         n = nu
      ELSE IF(iflag .eq. 2)THEN
         n = nus
      END IF
      gn = exp(gammln(n))
      
      IF((log(ani)-log(ao)).gt. 0.01 &
           .and.(log(cni)-log(co)).gt.0.001)THEN
         deltastr = (log(cni)-log(co))/(log(ani)-log(ao))
      ELSE
         deltastr = 1.
      ENDIF
      
      
      IF(iaspect .eq. 1) deltastr = 0.8
      IF(sphrflag .eq. 1) deltastr = 1.0
      
               
      if(deltastr.lt.0.55) then
         voltmp=(4./3.)*pi*co/ao**(deltastr)*ani**(2.+deltastr)* &
              (exp(gammln(n+deltastr+2.)))/gn
         deltastr=0.55
         ani=((3.*voltmp*gn)/ &
              (4.*pi*co/ao**(deltastr)*(exp(gammln(n+deltastr+2.)))))** &
              (1./(2.+deltastr))
         ai=max((n*ni*ani),1.e-20)
         ani=ai/(n*ni)
      else if (deltastr.gt.1.5) then
         voltmp=(4./3.)*pi*co/ao**(deltastr)*ani**(2.+deltastr)* &
              (exp(gammln(n+deltastr+2.)))/gn
         deltastr=1.5
         ani=((3.*voltmp*gn)/ &
              (4.*pi*co/ao**(deltastr)*(exp(gammln(n+deltastr+2.)))))** &
              (1./(2.+deltastr))
         ai=max((n*ni*ani),1.e-20)
         ani=ai/(n*ni)
         cni=co*(ani/ao)**deltastr
         ci=max((n*ni*cni),1.e-20)
         cni=ci/(n*ni)
      endif
      
    END SUBROUTINE DSTR_CHECK

    SUBROUTINE RHO_CHECK(iflag,deltastr,qi,ni,ani,cni,sphrflag,redden,betam,alphstr,alphv,rhobar)
      IMPLICIT NONE
      INTEGER sphrflag, redden, iflag
      REAL deltastr, qi, ni, ani, cni, n, gn
      REAL alphstr, alphv, betam
      REAL rhobar

      IF(iflag .eq. 1)THEN
         n = nu
      ELSE IF(iflag .eq. 2)THEN
         n = nus
      END IF
      gn = exp(gammln(n))
      
      betam = 2.+deltastr
      alphstr = co/ao**(deltastr)
      alphv = (4./3.)*pi*alphstr
      
      !     get avg ice density 
               
      rhobar = qi*gn/(ni*alphv* &
           ani**betam*exp(gammln(n+betam)))
               
      IF(sphrflag.eq.1) rhobar = 920.
      If(redden .eq. 1) rhoi = 500.

      IF(rhobar.gt.920.)THEN 
                      
         rhobar=920.
         ani=((qi*gn)/(rhobar*ni*alphv*&
              exp(gammln(n+betam))))**(1./betam)
         cni=co*(ani/ao)**deltastr
                  
      ELSE IF(rhobar.lt.50.)THEN
                      
         rhobar=50.
         ani=((qi*gn)/(rhobar*ni*alphv*&
              exp(gammln(n+betam))))**(1./betam)
         cni=co*(ani/ao)**deltastr
                      
      END IF

              

      
    END SUBROUTINE RHO_CHECK

    !     get rni (characteristic equivalent volume ice radius)
    SUBROUTINE R_CHECK(iflag,qi,ni,cni,ani,rni,rhobar,deltastr,betam,alphstr,alphv)
      IMPLICIT NONE
      INTEGER iflag
      REAL qi, ni, cni, ani, rhobar, deltastr, n, gn
      REAL rni
      REAL alphstr, alphv, betam

      IF(iflag .eq. 1)THEN
         n = nu
      ELSE IF(iflag .eq. 2)THEN
         n = nus
      END IF
      gn = exp(gammln(n))
      
      betam = 2.+deltastr
      alphstr = co/ao**(deltastr)
      alphv = (4./3.)*pi*alphstr

      rni = (qi*3./(ni*rhobar*4.*pi*&
           (exp(gammln(n+deltastr+2.))/gn)))**(1./3.)
      
      
      !     make sure rni is within reasonable bounds,
      
      IF(rni.lt.2.e-6)THEN
         
         rni=2.e-6
         ni=3.*qi*gn/(4.*pi*rhobar*rni**3.*&
              (exp(gammln(n+deltastr+2.))))
         ani=((qi*gn)/(rhobar*ni*alphv* &
              exp(gammln(n+betam))))**(1./betam) 
         cni=co*(ani/ao)**deltastr
         
      ELSE IF(rni.gt.2.e-3)THEN
         
         rni=2.e-3
         ni=3.*qi*gn/(4.*pi*rhobar*rni**3.* &
              (exp(gammln(n+deltastr+2.))))
         ani=((qi*gn)/(rhobar*ni*alphv* &
              exp(gammln(n+betam))))**(1./betam)
         cni=co*(ani/ao)**deltastr
         
      END IF
    END SUBROUTINE R_CHECK


!------------------------------------------------------------------
!     FUNCTION FINDGTP
!     
!     This function calculates Gi(T,P) for ice growth based on the
!     function in Pruppacher and Klett (1997)
!     
!------------------------------------------------------------------
      
      real function findgtp(temp,p0,fv,fh)
      
      implicit none
      real temp,p0,xdvstar,xk,eeq,esi_new,rv,L,gtp1,gtp2,dvstar
      real xkstar,dv,fv,fh
!      real polysvp1
      
      L = FLHS(temp)            !2.833e6               ! latent heat of sublimation
      dv = 0.211*(temp/273.15)**1.94 *(101300.25/p0)*1.0/100.0**2.0 ! vapor diffusivity
      xdvstar = dv*fv              ! vapor diffusivity in cm^2/s
!     xk = xkstar              ! thermal conductivity in W/m/K
      xk = (2.3823e-2 +7.1177e-5*(temp-273.15))*fh ! thermal conductivity
      eeq =  polysvp1(temp,1)           ! equilibrium vp
      rv = 461.5                ! H2O gas constant
!      WRITE(*,*) 'T', temp, 'P', p0, 'Ls', L, 'Kt', xk, 'Dv', dv, 'ei', 
!     +     eeq
      if(xdvstar.eq.0.)then
         findgtp=0.0
      else
         gtp1 = rv*temp/(xdvstar*eeq)
         gtp2 = (L**2.0/(xk*rv*temp**2.0)) - (L/(xk*temp))
         findgtp = 1.0/(gtp1+gtp2)
      endif
      return
      end function findgtp

!------------------------------------------------------------------------------

      real function capacitance_gamma(axis,deltastr,nu,alpha_s)
      
      implicit none
      real l,phi,a,c,axis,deltastr
      real a1,a2,b1,capone,c1,c2,d1,d2,nu
      real li,phicoli,phiplai,bcol,bpla
      real acol,apla,b2,alpha_s,beta_s
      real gammnu,gammad1,gammad2!,gammln

      a1 = 0.58
      a2 = 0.95
      b1 = 0.75

      l = axis

      li = ao/nu
      phicoli = 1.0
      phiplai = 1.0
      bcol = 1./deltastr
      bpla = deltastr
      acol = 1./phicoli * li**(1.-bcol)
      apla = phiplai * li**(1.-bpla)

!    Oblate Spheroid
      if (deltastr.le.1.) then
       a1 = 0.6369427        !0.58/0.9105999
       b1 = 0.0
       a2 = 0.57*a1   ! fit to spheroid
       b2 = 0.95   ! fit to spheroid
       alpha_s = apla
       beta_s = bpla-1.0

!     Prolate Spheroid
      elseif (deltastr.gt.1.) then
       a1 = 0.5714285   ! 0.58/1.015000
       b1 = -1.0
       a2 = 0.75*a1      ! fit to spheroid
       b2 = 0.82-1.0  ! fit to spheroid 
       alpha_s = apla
       beta_s = bpla-1.0

      endif

 
      if (deltastr.eq.1..and.nu.eq.1.) then
       gammnu = exp(gammln(nu))
       gammad1 = exp(gammln(nu+1.0))
       capacitance_gamma = l*gammad1/gammnu
       
      else if(deltastr.le.1.)then
       c1 = a1*alpha_s**b1
       c2 = a2*alpha_s**b2
       d1 = beta_s*b1+1.
       d2 = beta_s*b2+1.

       gammnu = exp(gammln(nu))
       gammad1 = exp(gammln(nu+d1))
       gammad2 = exp(gammln(nu+d2))
       capacitance_gamma = c1*l**d1 * gammad1/gammnu &
          + c2*l**d2 * gammad2/gammnu


      else if(deltastr.gt.1.) then
       c1 = a1*alpha_s**(b1+1.)
       c2 = a2*alpha_s**(b2+1.)
       d1 = beta_s*b1+deltastr
       d2 = beta_s*b2+deltastr

       gammnu = exp(gammln(nu))
       gammad1 = exp(gammln(nu+d1))
       gammad2 = exp(gammln(nu+d2))
       capacitance_gamma = c1*l**d1 * gammad1/gammnu &
          + c2*l**d2 * gammad2/gammnu

      endif

      return
      end function capacitance_gamma

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real function lenconv(axis,deltastr,nu,alpha_m,beta_m,alphstr,&
      rhoi)

      implicit none
      real l,phi,a,c,axis
      real a1,a2,b1,capone,c1,c2,d1,d2
      real deltastr,li,phicoli,phiplai,bcol,bpla
      real acol,apla,b2!,pi
      real gammnu,gammabeta,alpha_m,beta_m
      real alphstr,nu,rhoi

! NOTE: USE SEMI-MAJOR AXIS NOT MAJOR AXIS LENGTH!!

!      pi = 3.14159

      
      l = axis

      li = ao/nu
      phicoli = 1.0
      phiplai = 1.0
      bcol = 1./deltastr
      bpla = deltastr
      acol = 1./phicoli * li**(1.-bcol)
      apla = phiplai * li**(1.-bpla)

      if(deltastr.ne.1.) then
       beta_m = 2.+bpla
       alpha_m = 4./3.*pi*apla*rhoi
       alphstr = apla

      elseif (deltastr.eq.1.)then
       l = axis
       beta_m = 3.0
       alpha_m = 4./3.*pi*rhoi
       alphstr = 1.0
      endif

      gammnu = exp(gammln(nu))
      gammabeta = exp(gammln(nu+beta_m))
      lenconv = alpha_m*beta_m*l**(beta_m-1.)*gammabeta/gammnu

      return
      end function lenconv


!**********************************************************************
!     This function calculates the enthalpy of sublimation for       
!     water as per Bohren and Albrecht (1999, pg 197, eq 5.64). 
!     Arguement is in kelvin and the returned result is in mks 
!     Coded: JYH Feb 2005    
!**********************************************************************
      REAL FUNCTION FLHS (tmp)

      REAL tmp,cpv,ci,t0,xlvs0,delcv,delt,xlvs
        
      cpv = 1850.               ! specific heat of vapor const p
      ci = 2106.                ! specific heat of ice
      t0 = 273.15
      xlvs0 = 2.834e6           ! enthalpy of sublimation at 0 C 
      delcv = cpv - ci
      delt = t0 - tmp
      xlvs = (xlvs0 + delcv*delt) ! B&A equation 5.64
      FLHS = xlvs
  
      RETURN
      END FUNCTION FLHS

    
      subroutine getmasssize(celsiushab,lmean,alpham,betam,&
      alpha_a,beta_a,alpha_v,beta_v,qi,ipart)


      implicit none
      real celsiushab,lmean,alpham,betam,&
      alpha_a,beta_a,alpha_v,beta_v,qi
      integer ipart

      if(lmean*1.e6.lt.125.0)then
         if(ipart.eq.0)then
            alpham = 1.23e-3    ! needles
            betam = 1.8
         else if(ipart.eq.1)then
            alpham = 110.8      ! columns
            betam = 2.91
         endif
         beta_a = 2.0
         alpha_a = 0.684*100.**(beta_a-2.)
      else
         if(ipart.eq.0)then
            alpham = 1.23e-3    ! needles
            betam = 1.8
         else if(ipart.eq.1)then
            alpham =  2.739e-3  ! columns
            betam = 1.74
         endif
         beta_a = 1.5
         alpha_a = 0.0696*100.**(beta_a-2.)
      endif

      if(ipart.eq.2.or.ipart.eq.3)then
         if(lmean*1.e6.le.100.)then
            betam = 2.91
            alpham = 0.167*100.**betam/1000.
            beta_a = 2.0
            alpha_a = 0.684*100.**(beta_a-2.)
         else if(lmean*1.e6.gt.100..and.lmean*1.e6.le.300.)then
            betam = 1.91
            alpham = 0.00166*100.**betam/1000.
            beta_a = 1.5
            alpha_a = 0.0696*100.**(beta_a-2.)
         else if(lmean*1.e6.gt.300.)then
            betam = 1.74
            alpham = 0.000907*100.**betam/1000.
            beta_a = 1.414
            alpha_a = 0.0512*100.**(beta_a-2.)
         endif
      endif

      if(celsiushab.lt.-10.0.and.celsiushab.gt.-20.0)then
         if(ipart.eq.0)then
            alpham = .377e-2   !dendrites
            betam = 2.
            if(lmean*1.e6.le.90.)then
               beta_a = 1.85
               alpha_a = 0.24*100.**(beta_a-2.)
            else
               beta_a = 1.63
               alpha_a = 0.11*100.**(beta_a-2.)
            endif
! Alex's data
!     alpha_m = 2.423e-1 !stellar dendrites
!     beta_m = 2.53
!     alpha_m = 2.327e-2  ! classic dendrites
!     beta_m = 2.29
         else if(ipart.eq.1)then
            alpham = .8854     ! hex plates
            betam =  2.5
            if(lmean*1.e6.le.90.)then
               beta_a = 1.85
               alpha_a = 0.24*100.**(beta_a-2.)
            else
               beta_a = 2.0
               alpha_a = 0.65*100.**(beta_a-2.)
            endif
         endif

         if(ipart.eq.2)then
            betam = 2.45
            alpham = 0.00739*100.**betam/1000. ! hex plates
            if(lmean*1.e6.le.90.)then
               beta_a = 1.85
               alpha_a = 0.24*100.**(beta_a-2.)
            else
               beta_a = 1.63
               alpha_a = 0.11*100.**(beta_a-2.)
            endif
         endif

         if(ipart.eq.3)then
            if(lmean*1.e6.le.90.)then ! dendrites
               betam = 2.42
               alpham = 0.00583*100.**betam/1000.
               beta_a = 1.85
               alpha_a = 0.24*100.**(beta_a-2.)
            else
               betam = 1.67
               alpham = 0.00027*100.**betam/1000.
               beta_a = 1.63
               alpha_a = 0.11*100.**(beta_a-2.)
            endif
         endif
      endif


      if(ipart.eq.4)then
      if(celsiushab.lt.0.0.and.celsiushab.ge.-3)then
         qi = 0.179*2.          !cap for a column
         betam = 3.0            ! columns
         alpham = 0.0450*1000.**betam/(1000.*1000.)
         beta_v = 0.650
         alpha_v = 0.825*1000.**beta_v
         if(lmean*1.e6.le.100.)then
            beta_a = 2.0
            alpha_a = 0.684*100.**(beta_a-2.)
         else if(lmean*1.e6.gt.100..and.lmean*1.e6.le.300.)then
            beta_a = 1.5
            alpha_a = 0.0696*100.**(beta_a-2.)
         else if(lmean*1.e6.gt.300.)then
            beta_a = 1.414
            alpha_a = 0.0512*100.**(beta_a-2.)
         endif
      else if(celsiushab.lt.-3.0.and.celsiushab.ge.-8.)then
         qi = .1803*2.
         betam = 2.01           ! needles
         alpham = 0.0092*1000.**betam/(1000.*1000.)
         beta_v = 0.590
         alpha_v = 0.444*1000.**beta_v
         if(lmean*1.e6.le.100.)then
            beta_a = 2.0
            alpha_a = 0.684*100.**(beta_a-2.)
         else if(lmean*1.e6.gt.100..and.lmean*1.e6.le.300.)then
            beta_a = 1.5
            alpha_a = 0.0696*100.**(beta_a-2.)
         else if(lmean*1.e6.gt.300.)then
            beta_a = 1.414
            alpha_a = 0.0512*100.**(beta_a-2.)
         endif
      else if(celsiushab.lt.-8.0.and.celsiushab.ge.-11.0)then
         qi = 0.429*2.
         betam = 1.9            ! sector, broad-branched
         alpham = 0.037*1000.**betam/(1000.*1000.)
         beta_v = 0.410
         alpha_v = 0.690*1000.**beta_v
         if(lmean*1.e6.le.90.)then
            beta_a = 1.85
            alpha_a = 0.24*100.**(beta_a-2.)
         else
            beta_a = 2.0
            alpha_a = 0.65*100.**(beta_a-2.)
         endif
      else if(celsiushab.lt.-11.0.and.celsiushab.ge.-17.0)then
         qi = .3183*2.
         betam = 2.19           ! dendrites
         alpham = 0.0141*1000.**betam/(1000.*1000.)
         beta_v = 0.217
         alpha_v = 0.376*1000.**beta_v
         if(lmean*1.e6.le.90.)then
            beta_a = 1.85
            alpha_a = 0.24*100.**(beta_a-2.)
         else
            beta_a = 1.63
            alpha_a = 0.11*100.**(beta_a-2.)
         endif
      else if(celsiushab.lt.-17.0.and.celsiushab.ge.-21.0)then
         qi = .3183*2.
         betam = 2.19           ! radiating assemblage of dendrites
         alpham = 0.0141*1000.**betam/(1000.*1000.)
         beta_v = 0.217
         alpha_v = 0.376*1000.**beta_v
         if(lmean*1.e6.le.90.)then
            beta_a = 1.85
            alpha_a = 0.24*100.**(beta_a-2.)
         else
            beta_a = 1.63
            alpha_a = 0.11*100.**(beta_a-2.)
         endif
      else if(celsiushab.lt.-21.0)then
         qi = 0.429*2.
         betam = 1.9            ! cold type
         alpham = 0.037*1000.**betam/(1000.*1000.)
         beta_v = 0.410
         alpha_v = 0.690*1000.**beta_v
         if(lmean*1.e6.le.90.)then
            beta_a = 1.85
            alpha_a = 0.24*100.**(beta_a-2.)
         else
            beta_a = 2.0
            alpha_a = 0.65*100.**(beta_a-2.)
         endif
      endif
      endif
      
      return
      END SUBROUTINE GETMASSSIZE
    !___________________________________________________________________
      REAL FUNCTION POLYSVP1 (T,type)
!--------------------------------------------------------------------

!  Compute saturation vapor pressure by using 
! function from Goff and Gatch (1946)

!  Polysvp returned in units of pa.
!  T is input in units of K.
!  type refers to saturation with respect to liquid (0) or ice (1)

      implicit none

      real dum
      real T,t0
      integer type

      t0=273.16

! ice


! ice

      if (type.eq.1) then
         
!     Goff Gatch equation (good down to -100 C)

            polysvp1 = 10.**(-9.09718*(273.16/t-1.)-3.56654*&
            alog10(273.16/t)+0.876793*(1.-t/273.16)+&
            alog10(6.1071))*100.

        
            
         
      end if
      
      
! liquid
      if (type.eq.0) then

!     Goff Gatch equation, uncertain below -70 C

         polysvp1 = 10.**(-7.90298*(373.16/t-1.)+&
          5.02808*alog10(373.16/t)- &
          1.3816e-7*(10**(11.344*(1.-t/373.16))-1.)+&
          8.1328e-3*(10**(-3.49149*(373.16/t-1.))-1.)+ &
          alog10(1013.246))*100. 

    
      end if
      
      end function polysvp1
      
    
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
    
    REAL FUNCTION GAMMA(X)
      
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!CD    DOUBLE PRECISION FUNCTION DGAMMA(X)
!C----------------------------------------------------------------------
!C
!C THIS ROUTINE CALCULATES THE GAMMA FUNCTION FOR A REAL ARGUMENT X.
!C   COMPUTATION IS BASED ON AN ALGORITHM OUTLINED IN REFERENCE 1.
!C   THE PROGRAM USES RATIONAL FUNCTIONS THAT APPROXIMATE THE GAMMA
!C   FUNCTION TO AT LEAST 20 SIGNIFICANT DECIMAL DIGITS.  COEFFICIENTS
!C   FOR THE APPROXIMATION OVER THE INTERVAL (1,2) ARE UNPULBISHED.
!C   THOSE FOR THE APPROXIMATION FOR X .GE. 12 ARE FROM REFERENCE 2.
!C   THE ACCURACY ACHIEVED DEPENDS ON THE ARITHMETIC SYSTEM, THE
!C   COMPILER, THE INTRINSIC FUNCTIONS, AND PROPER SELECTION OF THE
!C   MACHINE-DEPENDENT CONSTANTS.
!C
!C
!C*******************************************************************
!C*******************************************************************
!C
!C EXPLANATION OF MACHINE-DEPENDENT CONSTANTS
!C
!C BETA   - RADIX FOR THE FLOATING-POINT REPRESENTATION
!C MAXEXP - THE SMALLEST POSITIVE POWER OF BETA THAT OVERFLOWS
!C XBIG   - THE LARGEST ARGUMENT FOR WHICH GAMMA(X) IS REPRESENTABLE
!C          IN THE MACHINE, I.E., THE SOLUTION TO THE EQUATION
!C                  GAMMA(XBIG) = BETA**MAXEXP
!C XINF   - THE LARGEST MACHINE REPRESENTABLE FLOATING-POINT NUMBER;
!C          APPROXIMATELY BETA**MAXEXP
!C EPS    - THE SMALLEST POSITIVE FLOATING-POINT NUMBER SUCH THAT
!C          1.0+EPS .GT. 1.0
!C XMININ - THE SMALLEST POSITIVE FLOATING-POINT NUMBER SUCH THAT
!C          1/XMININ IS MACHINE REPRESENTABLE
!C
!C     APPROXIMATE VALUES FOR SOME IMPORTANT MACHINES ARE:
!C
!C                            BETA       MAXEXP        XBIG
!C
!C CRAY-1         (S.P.)        2         8191        966.961
!C CYBER 180/855
!C   UNDER NOS    (S.P.)        2         1070        177.803
!C IEEE (IBM/XT,
!C   SUN, ETC.)   (S.P.)        2          128        35.040
!C IEEE (IBM/XT,
!C   SUN, ETC.)   (D.P.)        2         1024        171.624
!C IBM 3033       (D.P.)       16           63        57.574
!C VAX D-FORMAT   (D.P.)        2          127        34.844
!C VAX G-FORMAT   (D.P.)        2         1023        171.489
!C
!C                            XINF         EPS        XMININ
!C
!C CRAY-1         (S.P.)   5.45E+2465   7.11E-15    1.84E-2466
!C CYBER 180/855
!C   UNDER NOS    (S.P.)   1.26E+322    3.55E-15    3.14E-294
!C IEEE (IBM/XT,
!C   SUN, ETC.)   (S.P.)   3.40E+38     1.19E-7     1.18E-38
!C IEEE (IBM/XT,
!C   SUN, ETC.)   (D.P.)   1.79D+308    2.22D-16    2.23D-308
!C IBM 3033       (D.P.)   7.23D+75     2.22D-16    1.39D-76
!C VAX D-FORMAT   (D.P.)   1.70D+38     1.39D-17    5.88D-39
!C VAX G-FORMAT   (D.P.)   8.98D+307    1.11D-16    1.12D-308
!C
!C*******************************************************************
!C*******************************************************************
!C
!C ERROR RETURNS
!C
!C  THE PROGRAM RETURNS THE VALUE XINF FOR SINGULARITIES OR
!C     WHEN OVERFLOW WOULD OCCUR.  THE COMPUTATION IS BELIEVED
!C     TO BE FREE OF UNDERFLOW AND OVERFLOW.
!C
!C
!C  INTRINSIC FUNCTIONS REQUIRED ARE:
!C
!C     INT, DBLE, EXP, LOG, REAL, SIN
!C
!C
!C REFERENCES:  AN OVERVIEW OF SOFTWARE DEVELOPMENT FOR SPECIAL
!C              FUNCTIONS   W. J. CODY, LECTURE NOTES IN MATHEMATICS,
!C              506, NUMERICAL ANALYSIS DUNDEE, 1975, G. A. WATSON
!C              (ED.), SPRINGER VERLAG, BERLIN, 1976.
!C
!C              COMPUTER APPROXIMATIONS, HART, ET. AL., WILEY AND
!C              SONS, NEW YORK, 1968.
!C
!C  LATEST MODIFICATION: OCTOBER 12, 1989
!C
!C  AUTHORS: W. J. CODY AND L. STOLTZ
!C           APPLIED MATHEMATICS DIVISION
!C           ARGONNE NATIONAL LABORATORY
!C           ARGONNE, IL 60439
!C
!C----------------------------------------------------------------------
      INTEGER I,N
      LOGICAL PARITY
      REAL C(7),CONV,EPS,FACT,HALF,ONE,P(8),PI,Q(8),RES,SQRTPI,SUM,&
         TWELVE,TWO,X,XBIG,XDEN,XINF,XMININ,XNUM,Y,Y1,YSQ,Z,ZERO
!      DIMENSION C(7),P(8),Q(8)
!C----------------------------------------------------------------------
!C  MATHEMATICAL CONSTANTS
!C----------------------------------------------------------------------
      DATA ONE,HALF,TWELVE,TWO,ZERO/1.0E0,0.5E0,12.0E0,2.0E0,0.0E0/,&
          SQRTPI/0.9189385332046727417803297E0/,&
          PI/3.1415926535897932384626434E0/
!CD    DATA ONE,HALF,TWELVE,TWO,ZERO/1.0D0,0.5D0,12.0D0,2.0D0,0.0D0/,
!CD   1     SQRTPI/0.9189385332046727417803297D0/,
!CD   2     PI/3.1415926535897932384626434D0/
!C----------------------------------------------------------------------
!C  MACHINE DEPENDENT PARAMETERS
!C----------------------------------------------------------------------
      DATA XBIG,XMININ,EPS/35.040E0,1.18E-38,1.19E-7/, XINF/3.4E38/
!CD    DATA XBIG,XMININ,EPS/171.624D0,2.23D-308,2.22D-16/,
!CD   1     XINF/1.79D308/
!C----------------------------------------------------------------------
!C  NUMERATOR AND DENOMINATOR COEFFICIENTS FOR RATIONAL MINIMAX
!C     APPROXIMATION OVER (1,2).
!C----------------------------------------------------------------------
     DATA P/-1.71618513886549492533811E+0,2.47656508055759199108314E+1,&
            -3.79804256470945635097577E+2,6.29331155312818442661052E+2,&
            8.66966202790413211295064E+2,-3.14512729688483675254357E+4,&
            -3.61444134186911729807069E+4,6.64561438202405440627855E+4/
     DATA Q/-3.08402300119738975254353E+1,3.15350626979604161529144E+2,&
            -1.01515636749021914166146E+3,-3.10777167157231109440444E+3&
            ,2.25381184209801510330112E+4,4.75584627752788110767815E+3,&
            -1.34659959864969306392456E+5,-1.15132259675553483497211E+5/
!CD    DATA P/-1.71618513886549492533811D+0,2.47656508055759199108314D+1,
!CD   1       -3.79804256470945635097577D+2,6.29331155312818442661052D+2,
!CD   2       8.66966202790413211295064D+2,-3.14512729688483675254357D+4,
!CD   3       -3.61444134186911729807069D+4,6.64561438202405440627855D+4/
!CD    DATA Q/-3.08402300119738975254353D+1,3.15350626979604161529144D+2,
!CD   1      -1.01515636749021914166146D+3,-3.10777167157231109440444D+3,
!CD   2        2.25381184209801510330112D+4,4.75584627752788110767815D+3,
!CD   3      -1.34659959864969306392456D+5,-1.15132259675553483497211D+5/
!C----------------------------------------------------------------------
!C  COEFFICIENTS FOR MINIMAX APPROXIMATION OVER (12, INF).
!C----------------------------------------------------------------------
      DATA C/-1.910444077728E-03,8.4171387781295E-04,&
           -5.952379913043012E-04,7.93650793500350248E-04,&
           -2.777777777777681622553E-03,8.333333333333333331554247E-02&
           ,5.7083835261E-03/
!C----------------------------------------------------------------------
      CONV(I) = REAL(I)
!CD    CONV(I) = DBLE(I)
      PARITY=.FALSE.
      FACT=ONE
      N=0
      Y=X
      IF(Y.LE.ZERO)THEN
!C----------------------------------------------------------------------
!C  ARGUMENT IS NEGATIVE
!C----------------------------------------------------------------------
        Y=-X
        Y1=AINT(Y)
        RES=Y-Y1
        IF(RES.NE.ZERO)THEN
          IF(Y1.NE.AINT(Y1*HALF)*TWO)PARITY=.TRUE.
          FACT=-PI/SIN(PI*RES)
          Y=Y+ONE
        ELSE
          RES=XINF
          GOTO 900
        ENDIF
      ENDIF
!C----------------------------------------------------------------------
!C  ARGUMENT IS POSITIVE
!C----------------------------------------------------------------------
      IF(Y.LT.EPS)THEN
!C----------------------------------------------------------------------
!C  ARGUMENT .LT. EPS
!C----------------------------------------------------------------------
        IF(Y.GE.XMININ)THEN
          RES=ONE/Y
        ELSE
          RES=XINF
          GOTO 900
        ENDIF
      ELSEIF(Y.LT.TWELVE)THEN
        Y1=Y
        IF(Y.LT.ONE)THEN
!C----------------------------------------------------------------------
!C  0.0 .LT. ARGUMENT .LT. 1.0
!C----------------------------------------------------------------------
          Z=Y
          Y=Y+ONE
        ELSE
!C----------------------------------------------------------------------
!C  1.0 .LT. ARGUMENT .LT. 12.0, REDUCE ARGUMENT IF NECESSARY
!C----------------------------------------------------------------------
          N=INT(Y)-1
          Y=Y-CONV(N)
          Z=Y-ONE
        ENDIF
!C----------------------------------------------------------------------
!C  EVALUATE APPROXIMATION FOR 1.0 .LT. ARGUMENT .LT. 2.0
!C----------------------------------------------------------------------
        XNUM=ZERO
        XDEN=ONE
        DO 260 I=1,8
          XNUM=(XNUM+P(I))*Z
          XDEN=XDEN*Z+Q(I)
  260   CONTINUE
        RES=XNUM/XDEN+ONE
        IF(Y1.LT.Y)THEN
!C----------------------------------------------------------------------
!C  ADJUST RESULT FOR CASE  0.0 .LT. ARGUMENT .LT. 1.0
!C----------------------------------------------------------------------
!C  ADJUST RESULT FOR CASE  2.0 .LT. ARGUMENT .LT. 12.0
!C----------------------------------------------------------------------
          DO 290 I=1,N
            RES=RES*Y
            Y=Y+ONE
  290     CONTINUE
        ENDIF
      ELSE
!C----------------------------------------------------------------------
!C  EVALUATE FOR ARGUMENT .GE. 12.0,
!C----------------------------------------------------------------------
        IF(Y.LE.XBIG)THEN
          YSQ=Y*Y
          SUM=C(7)
          DO 350 I=1,6
            SUM=SUM/YSQ+C(I)
  350     CONTINUE
          SUM=SUM/Y-Y+SQRTPI
          SUM=SUM+(Y-HALF)*LOG(Y)
          RES=EXP(SUM)
        ELSE
          RES=XINF
          GOTO 900
        ENDIF
      ENDIF
!C----------------------------------------------------------------------
!C  FINAL ADJUSTMENTS AND RETURN
!C----------------------------------------------------------------------
      IF(PARITY)RES=-RES
      IF(FACT.NE.ONE)RES=FACT/RES
  900 GAMMA=RES
!CD900 DGAMMA = RES
      RETURN
!C ---------- LAST LINE OF GAMMA ----------
      END function gamma
  
  
end MODULE MODULE_MP_SULIAHARRINGTON
