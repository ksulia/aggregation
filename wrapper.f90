PROGRAM MICROWRAPPER


  USE MODULE_MP_SULIAHARRINGTON
  USE NETCDF_PLUGIN

  IMPLICIT NONE

  !NETCDF 
  CHARACTER(150) :: FNAME
  INTEGER :: NCID,TID,VID

  ! DEFINE VARIABLES
  INTEGER :: ITIMESTEP
  INTEGER, PARAMETER :: NTIMES=1000

  REAL, DIMENSION(1,1,1) :: ICEDEP,ICESUB,RAINEVAP,SNOWEVAP,SNOWMELT
  REAL, DIMENSION(1,1,1) :: SNOWDEP,SNOWSUB,SNOWACCR,CLOUDCOND
  REAL, DIMENSION(1,1,1) :: CLOUDEVAP,ICEMELT,ICENUC,RAINFRZ,CLOUDFRZ
  REAL, DIMENSION(1,1,1) :: PHI,RHOICE,RELH,CPLX,PHIS,RHOS,CPLXS

  REAL, DIMENSION(1,1,1) :: QV_CURR,TH !  &
  REAL, DIMENSION(1,1,1) :: QNC_CURR,QC_CURR,QNR_CURR,QR_CURR,QNI_CURR,QI_CURR,QNS_CURR,QS_CURR
  REAL, DIMENSION(1,1,1) :: AICE_CURR, CICE_CURR, ASNOW_CURR, CSNOW_CURR, RHO, PI_PHY, P,DZ8W,HT,W

  REAL :: DT, DQV, DTEMP

  REAL :: T_INIT, P_INIT

  !OUTPUT VARIABLES
  REAL, DIMENSION(NTIMES+1) :: QVOUT,QCOUT,QROUT,QSOUT,QIOUT,QNCOUT,QNROUT,QNSOUT,QNIOUT,TIMESOUT,&
          RHOICEOUT,PHIOUT,RELHOUT,CPLXOUT,PHISOUT,RHOSOUT,CPLXSOUT
  REAL, DIMENSION(NTIMES+1) :: THOUT,POUT,AICEOUT,CICEOUT,ASNOWOUT,CSNOWOUT,ICEDEPOUT,ICESUBOUT,RAINEVAPOUT,SNOWEVAPOUT,SNOWMELTOUT
  REAL, DIMENSION(NTIMES+1) :: SNOWDEPOUT,SNOWSUBOUT,SNOWACCROUT,CLOUDCONDOUT,CLOUDEVAPOUT,ICEMELTOUT,&
          ICENUCOUT,RAINFRZOUT,CLOUDFRZOUT
  !INITIALIZE VARIABLES
  WRITE(FNAME,"(A9)") "OUTPUT.NC"

  DQV = 1.E-4
  DTEMP = 0.27
  DT = 1.
  !NTIMES = 1000.
  T_INIT = 260. !K
  P_INIT = 800. !HPA
  QV_CURR=1.E-3 !KG/KG
  QC_CURR=0.
  QR_CURR=0.
  QS_CURR=0.
  QI_CURR=0.
  QNC_CURR=0.
  QNR_CURR=0.
  QNS_CURR=0.
  QNI_CURR=0.
  AICE_CURR =0.
  CICE_CURR=0.
  ASNOW_CURR =0.
  CSNOW_CURR=0.
  DZ8W = 100. !M PROBABLY IRRELEVANT ONCE SEDIMENATION IS SHUT OFF
  HT= 1000. !M
  W = 0. ! COULD CHANGE TO MIMIC PARCEL MODEL BUT WOULD HAVE TO INCLUDE EXPANSION TERM FOR THETA AND P TENDENCY

  !SET DIAGNOSTICS TO 0
  ICEDEP=0.
  ICESUB=0.
  RAINEVAP=0.
  SNOWEVAP=0.
  SNOWMELT=0.
  SNOWDEP=0.
  SNOWSUB=0.
  SNOWACCR=0.
  CLOUDCOND=0.
  CLOUDEVAP=0.
  ICEMELT=0.
  ICENUC=0.
  RAINFRZ=0.
  CLOUDFRZ=0.
  PHI=0.
  RHOICE=0.
  RELH=0.
  CPLX = 0.
  PHIS=0.
  RHOS=0.
  CPLXS=0.

  !COMPUTE OTHER VARIABLES
  RHO=P_INIT*100./T_INIT/287.04
  P=P_INIT*100.
  PI_PHY=(P_INIT/1000)**(287.04/1004.) !COULD REPLACE WITH T-DEPENDENT VALUES IF YOU WISH
  TH=T_INIT/PI_PHY

  QVOUT(1)=QV_CURR(1,1,1)
  TIMESOUT(1)=0.
  QCOUT(1)=QC_CURR(1,1,1)
  QROUT(1)=QR_CURR(1,1,1)
  QSOUT(1)=QS_CURR(1,1,1)
  QIOUT(1)=QI_CURR(1,1,1)
  QNCOUT(1)=QNC_CURR(1,1,1)
  QNROUT(1)=QNR_CURR(1,1,1)
  QNSOUT(1)=QNS_CURR(1,1,1)
  QNIOUT(1)=QNI_CURR(1,1,1)

  AICEOUT(1)=AICE_CURR(1,1,1)
  CICEOUT(1)=CICE_CURR(1,1,1)
  ASNOWOUT(1)=ASNOW_CURR(1,1,1)
  CSNOWOUT(1)=CSNOW_CURR(1,1,1)


  THOUT(1)=TH(1,1,1)
  POUT(1)=P(1,1,1)

  ICEDEPOUT(1)=ICEDEP(1,1,1)
  ICESUBOUT(1)=ICESUB(1,1,1)
  RAINEVAPOUT(1)=RAINEVAP(1,1,1)
  SNOWEVAPOUT(1)=SNOWEVAP(1,1,1)
  SNOWMELTOUT(1)=SNOWMELT(1,1,1)
  SNOWDEPOUT(1)=SNOWDEP(1,1,1)
  SNOWACCROUT(1)=SNOWACCR(1,1,1)
  CLOUDCONDOUT(1)=CLOUDCOND(1,1,1)
  CLOUDEVAPOUT(1)=CLOUDEVAP(1,1,1)
  ICEMELTOUT(1)=ICEMELT(1,1,1)
  ICENUCOUT(1)=ICENUC(1,1,1)
  RAINFRZOUT(1)=RAINFRZ(1,1,1)
  CLOUDFRZOUT(1)=CLOUDFRZ(1,1,1)
  PHIOUT(1)=PHI(1,1,1)
  RHOICEOUT(1)=RHOICE(1,1,1)
  RELHOUT(1)=RELH(1,1,1)
  CPLXOUT(1)=CPLX(1,1,1)
  PHISOUT(1)=PHIS(1,1,1)
  RHOSOUT(1)=RHOS(1,1,1)
  CPLXSOUT(1)=CPLXS(1,1,1)
  CALL SULIAHARRINGTON_INIT
  !BEGIN TIME ITERATION

  DO ITIMESTEP = 1,NTIMES
     !    WRITE (*,*) "DOING TIME ",itimestep*dt," seconds of ",ntimes*dt," seconds"
     
!     QV_CURR(1,1,1) = QV_CURR(1,1,1) + DQV*DT
     TH(1,1,1) = TH(1,1,1) - DTEMP*DT!T_INIT*(1.E5/P_INIT/100.)**(287.05/1005.)

     QV_CURR(1,1,1) = QV_CURR(1,1,1) + DQV*DT
!     TH(1,1,1) = T_INIT*(1.E5/P_INIT/100.)**(287.05/1005.)

    CALL MP_SULIAHARRINGTON(                       &
    ITIMESTEP=ITIMESTEP,                &  !*
    TH=TH,                              &  !*
    QV=QV_CURR,                         &  !*
    QC=QC_CURR,                         &  !*
    QR=QR_CURR,                         &
    QI=QI_CURR,                     &  !*
    QS=QS_CURR,                         &
    NC=QNC_CURR,                        &
    NR=QNR_CURR,                        &
    NI=QNI_CURR,                    &  !*
    NS=QNS_CURR,                        &
    AI=AICE_CURR,                     &
    CI=CICE_CURR,                     &
    AS=ASNOW_CURR,                     &
    CS=CSNOW_CURR,                     &
    RHO=RHO,                            &  !*
    PII=PI_PHY,                         &  !*
    P=P,                                &  !*
    DT=DT,                              &  !*
    DZ=DZ8W,                            &  !* !HM
    HT=HT,                              &  !*
    W=W,                                &  !*
    ICEDEP=ICEDEP,                      &
    ICESUB=ICESUB,                      &
    RAINEVAP=RAINEVAP,                  &
    SNOWEVAP=SNOWEVAP,                  &
    SNOWMELT=SNOWMELT,                  &
    SNOWDEP=SNOWDEP,                    &
    SNOWSUB=SNOWSUB,                    &
    SNOWACCR=SNOWACCR,                  &
    CLOUDCOND=CLOUDCOND,                &
    CLOUDEVAP=CLOUDEVAP,                &
    ICEMELT=ICEMELT,                    &
    ICENUC=ICENUC,                      &
    RAINFRZ=RAINFRZ,                    &
    CLOUDFRZ=CLOUDFRZ,                  &
    PHI=PHI,                            &
    RHOICE=RHOICE,                      &
    RELH=RELH,                          &
    CPLX=CPLX,                          &
    PHIS=PHIS,                          &
    RHOS=RHOS,                          &
    CPLXS=CPLXS                         &
    ,IDS=1,IDE=1, JDS=1,JDE=1, KDS=1,KDE=1 &
    ,IMS=1,IME=1, JMS=1,JME=1, KMS=1,KME=1 &
    ,ITS=1,ITE=1, JTS=1,JTE=1, KTS=1,KTE=1 &
    )


    QVOUT(ITIMESTEP+1)=QV_CURR(1,1,1)
    QCOUT(ITIMESTEP+1)=QC_CURR(1,1,1)
    QROUT(ITIMESTEP+1)=QR_CURR(1,1,1)
    QSOUT(ITIMESTEP+1)=QS_CURR(1,1,1)
    QIOUT(ITIMESTEP+1)=QI_CURR(1,1,1)
    QNCOUT(ITIMESTEP+1)=QNC_CURR(1,1,1)
    QNROUT(ITIMESTEP+1)=QNR_CURR(1,1,1)
    QNSOUT(ITIMESTEP+1)=QNS_CURR(1,1,1)
    QNIOUT(ITIMESTEP+1)=QNI_CURR(1,1,1)

    AICEOUT(ITIMESTEP+1)=AICE_CURR(1,1,1)
    CICEOUT(ITIMESTEP+1)=CICE_CURR(1,1,1)
    ASNOWOUT(ITIMESTEP+1)=ASNOW_CURR(1,1,1)
    CSNOWOUT(ITIMESTEP+1)=CSNOW_CURR(1,1,1)

    THOUT(ITIMESTEP+1)=TH(1,1,1)
    POUT(ITIMESTEP+1)=P(1,1,1)

    ICEDEPOUT(ITIMESTEP+1)=ICEDEP(1,1,1)
    ICESUBOUT(ITIMESTEP+1)=ICESUB(1,1,1)
    RAINEVAPOUT(ITIMESTEP+1)=RAINEVAP(1,1,1)
    SNOWEVAPOUT(ITIMESTEP+1)=SNOWEVAP(1,1,1)
    SNOWMELTOUT(ITIMESTEP+1)=SNOWMELT(1,1,1)
    SNOWDEPOUT(ITIMESTEP+1)=SNOWDEP(1,1,1)
    SNOWACCROUT(ITIMESTEP+1)=SNOWACCR(1,1,1)
    CLOUDCONDOUT(ITIMESTEP+1)=CLOUDCOND(1,1,1)
    CLOUDEVAPOUT(ITIMESTEP+1)=CLOUDEVAP(1,1,1)
    ICEMELTOUT(ITIMESTEP+1)=ICEMELT(1,1,1)
    ICENUCOUT(ITIMESTEP+1)=ICENUC(1,1,1)
    RAINFRZOUT(ITIMESTEP+1)=RAINFRZ(1,1,1)
    CLOUDFRZOUT(ITIMESTEP+1)=CLOUDFRZ(1,1,1)
    PHIOUT(ITIMESTEP+1)=PHI(1,1,1)
    RHOICEOUT(ITIMESTEP+1)=RHOICE(1,1,1)
    RELHOUT(ITIMESTEP+1)=RELH(1,1,1)
    CPLXOUT(ITIMESTEP+1)=CPLX(1,1,1)
    PHISOUT(ITIMESTEP+1)=PHIS(1,1,1)
    RHOSOUT(ITIMESTEP+1)=RHOS(1,1,1)
    CPLXSOUT(ITIMESTEP+1)=CPLXS(1,1,1)
    TIMESOUT(ITIMESTEP+1)=DT*ITIMESTEP
  ENDDO

  !OUTPUT USING NETCDF ROUTINES
  CALL CREATENETCDF(NCID,TRIM(FNAME))
  CALL DEFINEDIMS(NCID,TID,NTIMES+1)
  CALL DEFINEVAR1D(NCID,"TIMES",TID,VID)
  CALL DEFINEVAR1D(NCID,"QVAPOR",TID,VID)
  CALL DEFINEVAR1D(NCID,"QCLOUD",TID,VID)
  CALL DEFINEVAR1D(NCID,"QRAIN",TID,VID)
  CALL DEFINEVAR1D(NCID,"QSNOW",TID,VID)
  CALL DEFINEVAR1D(NCID,"QICE",TID,VID)
  CALL DEFINEVAR1D(NCID,"QNCLOUD",TID,VID)
  CALL DEFINEVAR1D(NCID,"QNRAIN",TID,VID)
  CALL DEFINEVAR1D(NCID,"QNSNOW",TID,VID)
  CALL DEFINEVAR1D(NCID,"QNICE",TID,VID)

  CALL DEFINEVAR1D(NCID,"AICE",TID,VID)
  CALL DEFINEVAR1D(NCID,"CICE",TID,VID)
  CALL DEFINEVAR1D(NCID,"ASNOW",TID,VID)
  CALL DEFINEVAR1D(NCID,"CSNOW",TID,VID)
  CALL DEFINEVAR1D(NCID,"THETA",TID,VID)
  CALL DEFINEVAR1D(NCID,"PRESS",TID,VID)

  CALL DEFINEVAR1D(NCID,"ICEDEP",TID,VID)
  CALL DEFINEVAR1D(NCID,"ICESUB",TID,VID)
  CALL DEFINEVAR1D(NCID,"RAINEVAP",TID,VID)
  CALL DEFINEVAR1D(NCID,"SNOWEVAP",TID,VID)
  CALL DEFINEVAR1D(NCID,"SNOWMELT",TID,VID)
  CALL DEFINEVAR1D(NCID,"SNOWDEP",TID,VID)
  CALL DEFINEVAR1D(NCID,"SNOWSUB",TID,VID)
  CALL DEFINEVAR1D(NCID,"SNOWACCR",TID,VID)
  CALL DEFINEVAR1D(NCID,"CLOUDCOND",TID,VID)
  CALL DEFINEVAR1D(NCID,"CLOUDEVAP",TID,VID)
  CALL DEFINEVAR1D(NCID,"ICEMELT",TID,VID)
  CALL DEFINEVAR1D(NCID,"ICENUC",TID,VID)
  CALL DEFINEVAR1D(NCID,"RAINFRZ",TID,VID)
  CALL DEFINEVAR1D(NCID,"CLOUDFRZ",TID,VID)
  CALL DEFINEVAR1D(NCID,"PHI",TID,VID)
  CALL DEFINEVAR1D(NCID,"RHOICE",TID,VID)
  CALL DEFINEVAR1D(NCID,"RELH",TID,VID)
  CALL DEFINEVAR1D(NCID,"CPLX",TID,VID)
  CALL DEFINEVAR1D(NCID,"PHIS",TID,VID)
  CALL DEFINEVAR1D(NCID,"RHOS",TID,VID)
  CALL DEFINEVAR1D(NCID,"CPLXS",TID,VID)

  CALL ENDDEF(NCID)
  CALL CLOSENETCDF(NCID)

  CALL OPENNETCDF(NCID,TRIM(FNAME))
  CALL WRITEVAR1D(NCID,TIMESOUT,NTIMES+1,"TIMES")
  CALL WRITEVAR1D(NCID,QVOUT,NTIMES+1,"QVAPOR")
  CALL WRITEVAR1D(NCID,QCOUT,NTIMES+1,"QCLOUD")
  CALL WRITEVAR1D(NCID,QROUT,NTIMES+1,"QRAIN")
  CALL WRITEVAR1D(NCID,QSOUT,NTIMES+1,"QSNOW")
  CALL WRITEVAR1D(NCID,QIOUT,NTIMES+1,"QICE")

  CALL WRITEVAR1D(NCID,QNCOUT,NTIMES+1,"QNCLOUD")
  CALL WRITEVAR1D(NCID,QNROUT,NTIMES+1,"QNRAIN")
  CALL WRITEVAR1D(NCID,QNSOUT,NTIMES+1,"QNSNOW")
  CALL WRITEVAR1D(NCID,QNIOUT,NTIMES+1,"QNICE")

  CALL WRITEVAR1D(NCID,AICEOUT,NTIMES+1,"AICE")
  CALL WRITEVAR1D(NCID,CICEOUT,NTIMES+1,"CICE")
  CALL WRITEVAR1D(NCID,ASNOWOUT,NTIMES+1,"ASNOW")
  CALL WRITEVAR1D(NCID,CSNOWOUT,NTIMES+1,"CSNOW")
  CALL WRITEVAR1D(NCID,THOUT,NTIMES+1,"THETA")
  CALL WRITEVAR1D(NCID,POUT,NTIMES+1,"PRESS")

  CALL WRITEVAR1D(NCID,ICEDEPOUT,NTIMES+1,"ICEDEP")
  CALL WRITEVAR1D(NCID,ICESUBOUT,NTIMES+1,"ICESUB")
  CALL WRITEVAR1D(NCID,RAINEVAPOUT,NTIMES+1,"RAINEVAP")
  CALL WRITEVAR1D(NCID,SNOWEVAPOUT,NTIMES+1,"SNOWEVAP")
  CALL WRITEVAR1D(NCID,SNOWMELTOUT,NTIMES+1,"SNOWMELT")
  CALL WRITEVAR1D(NCID,SNOWDEPOUT,NTIMES+1,"SNOWDEP")
  CALL WRITEVAR1D(NCID,SNOWSUBOUT,NTIMES+1,"SNOWSUB")
  CALL WRITEVAR1D(NCID,SNOWACCROUT,NTIMES+1,"SNOWACCR")
  CALL WRITEVAR1D(NCID,CLOUDCONDOUT,NTIMES+1,"CLOUDCOND")
  CALL WRITEVAR1D(NCID,CLOUDEVAPOUT,NTIMES+1,"CLOUDEVAP")
  CALL WRITEVAR1D(NCID,ICEMELTOUT,NTIMES+1,"ICEMELT")
  CALL WRITEVAR1D(NCID,ICENUCOUT,NTIMES+1,"ICENUC")
  CALL WRITEVAR1D(NCID,RAINFRZOUT,NTIMES+1,"RAINFRZ")
  CALL WRITEVAR1D(NCID,CLOUDFRZOUT,NTIMES+1,"CLOUDFRZ")
  CALL WRITEVAR1D(NCID,PHIOUT,NTIMES+1,"PHI")
  CALL WRITEVAR1D(NCID,RHOICEOUT,NTIMES+1,"RHOICE")
  CALL WRITEVAR1D(NCID,RELHOUT,NTIMES+1,"RELH")
  CALL WRITEVAR1D(NCID,CPLXOUT,NTIMES+1,"CPLX")
  CALL WRITEVAR1D(NCID,PHISOUT,NTIMES+1,"PHIS")
  CALL WRITEVAR1D(NCID,RHOSOUT,NTIMES+1,"RHOS")
  CALL WRITEVAR1D(NCID,CPLXSOUT,NTIMES+1,"CPLXS")
  CALL CLOSENETCDF(NCID)

END PROGRAM MICROWRAPPER
