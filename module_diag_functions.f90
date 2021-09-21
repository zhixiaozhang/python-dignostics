


MODULE diag_functions

CONTAINS

  
  
  
  
  
  
  
  
  
  
  FUNCTION calc_rh ( p, t, qv ) result ( rh )
    
    IMPLICIT NONE
 
    REAL, INTENT(IN) :: p, t, qv
    REAL :: rh

    
    
    REAL :: evapor, es
  
    

      evapor = p/((0.622/qv)+1.0);
      es = 611.2*exp(17.67*(t-273.15)/(t-29.65));
      rh=100.*evapor/es
      IF (rh .gt. 100.) THEN
        rh=100.
      ENDIF

  END FUNCTION calc_rh



  
  
  
  
  
  
  
  
  
  
  FUNCTION uv_wind ( u, v ) result ( wind_speed )
 
    IMPLICIT NONE
 
    REAL, INTENT(IN) :: u, v
    REAL :: wind_speed

    wind_speed = sqrt( u*u + v*v )

  END FUNCTION uv_wind


  
  
  
  
  
  
  
  

  
  
  FUNCTION Theta ( t, p )
  IMPLICIT NONE

     
     
     REAL, INTENT ( IN ) :: t
     REAL, INTENT ( IN ) :: p
     REAL                :: theta

     REAL :: Rd 
     REAL :: Cp 
     REAL :: p0 
  
     Rd =  287.04
     Cp = 1004.67
     p0 = 1000.00


     
     theta = t * ( (p0/p)**(Rd/Cp) )
  
  END FUNCTION Theta



  
  
  
  
  
  
  
  
  
  
  FUNCTION Thetae ( tK, p, rh, mixr )
  IMPLICIT NONE

     
     
     REAL :: tK        
     REAL :: p         
     REAL :: rh        
     REAL :: mixr      
     REAL :: te        
     REAL :: thetae    
  
     REAL, PARAMETER :: R  = 287.04         
     REAL, PARAMETER :: P0 = 1000.0         
     REAL, PARAMETER :: lv = 2.54*(10**6)   
                                            
     REAL, PARAMETER :: cp = 1004.67        
                                            
     REAL :: tlc                            
  
     
     
     tlc = TLCL ( tK, rh )
  
     
     
     thetae = (tK * (p0/p)**( (R/Cp)*(1.- ( (.28E-3)*mixr*1000.) ) ) )* &
                 exp( (((3.376/tlc)-.00254))*&
                    (mixr*1000.*(1.+(.81E-3)*mixr*1000.)) )
  
  END FUNCTION Thetae



  
  
  
  
  
  
  
  
  
  
  
  
  
  
  FUNCTION The2T ( thetaeK, pres, flag ) result ( tparcel )
  IMPLICIT NONE
  
     
     
     REAL,    INTENT     ( IN ) :: thetaeK
     REAL,    INTENT     ( IN ) :: pres
     LOGICAL, INTENT ( INOUT )  :: flag
     REAL                       :: tparcel
  
     REAL :: thetaK
     REAL :: tovtheta
     REAL :: tcheck
     REAL :: svpr, svpr2
     REAL :: smixr, smixr2
     REAL :: thetae_check, thetae_check2
     REAL :: tguess_2, correction
  
     LOGICAL :: found
     INTEGER :: iter
  
     REAL :: R     
     REAL :: Cp    
     REAL :: kappa 
     REAL :: Lv    
  
     R     = 287.04
     Cp    = 1004.67
     Kappa = R/Cp
     Lv    = 2.500E+6

     
     
     tovtheta = (pres/100000.0)**(r/cp)
     tparcel  = thetaeK/exp(lv*.012/(cp*295.))*tovtheta

     iter = 1
     found = .false.
     flag = .false.

     DO
        IF ( iter > 105 ) EXIT

        tguess_2 = tparcel + REAL ( 1 )

        svpr   = 6.122 * exp ( (17.67*(tparcel-273.15)) / (tparcel-29.66) )
        smixr  = ( 0.622*svpr ) / ( (pres/100.0)-svpr )
        svpr2  = 6.122 * exp ( (17.67*(tguess_2-273.15)) / (tguess_2-29.66) )
        smixr2 = ( 0.622*svpr2 ) / ( (pres/100.0)-svpr2 )

        
        
        
        
        
        
        
        
        
        
        

        
        
        

        
        thetae_check  = Thetae ( tparcel,  pres/100., 100., smixr  )
        thetae_check2 = Thetae ( tguess_2, pres/100., 100., smixr2 )

        
        
        IF ( ABS (thetaeK-thetae_check) < .001) THEN
           found = .true.
           flag  = .true.
           EXIT
        END IF

        
        

        
        correction = ( thetaeK-thetae_check ) / ( thetae_check2-thetae_check )
        tparcel = tparcel + correction

        iter = iter + 1
     END DO

     
     
     
     

  END FUNCTION The2T



  
  
  
  
  
  
  
  
  
  
  FUNCTION VirtualTemperature ( tK, w ) result ( Tv )
  IMPLICIT NONE

     
     real, intent ( in ) :: tK 
     real, intent ( in ) :: w  
     real                :: Tv 

     Tv = tK * ( 1.0 + (w/0.622) ) / ( 1.0 + w )

  END FUNCTION VirtualTemperature




  
  
  
  
  
  
  
  
  
  
  
  
  
  
  FUNCTION SaturationMixingRatio ( tK, p ) result ( ws )

    IMPLICIT NONE

    REAL, INTENT ( IN ) :: tK
    REAL, INTENT ( IN ) :: p
    REAL                :: ws

    REAL :: es

    es = 6.122 * exp ( (17.67*(tK-273.15))/ (tK-29.66) )
    ws = ( 0.622*es ) / ( (p/100.0)-es )

  END FUNCTION SaturationMixingRatio



  
  
  
  
  
  
  

  
  
  
  
  
  FUNCTION TLCL ( tk, rh )
    
    IMPLICIT NONE
 
    REAL, INTENT ( IN ) :: tK   
    REAL, INTENT ( IN ) :: rh   
    REAL                :: tlcl
    
    REAL :: denom, term1, term2

    term1 = 1.0 / ( tK - 55.0 )
    IF ( rh > REAL (0) ) THEN
      term2 = ( LOG (rh/100.0)  / 2840.0 )
    ELSE
      term2 = ( LOG (0.001/1.0) / 2840.0 )
    END IF
    denom = term1 - term2
    tlcl = ( 1.0 / denom ) + REAL ( 55 ) 

  END FUNCTION TLCL



  
  
  
  
  
  
  
  
  
  
  FUNCTION Pwat  ( nz, qv, qc, dz8w, rho )

    IMPLICIT NONE

     
     
     INTEGER, INTENT ( IN ) :: nz          
     REAL, INTENT ( IN )    :: qv   ( nz ) 
     REAL, INTENT ( IN )    :: qc   ( nz ) 
     REAL, INTENT ( IN )    :: dz8w ( nz ) 
     REAL, INTENT ( IN )    :: rho  ( nz ) 
     REAL                   :: Pwat        
     INTEGER                :: k           

     
     
     Pwat=0
     DO k = 1, nz
       
        
        
        Pwat = Pwat + qv(k) * dz8w(k) * rho(k)
     ENDDO
             
  END FUNCTION Pwat
 


  
  
  
  
  
  
  
  
  
  
  
  

  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
    FUNCTION Buoyancy ( nz, tk, rh, p, hgt, sfc, cape, cin, zlcl, zlfc, zel, zmu,  &
                        parcel ) result (ostat) 
  
      IMPLICIT NONE
  
      INTEGER, INTENT ( IN )  :: nz          
      INTEGER, INTENT ( IN )  :: sfc         
      REAL,    INTENT ( IN )  :: tk   ( nz ) 
      REAL,    INTENT ( IN )  :: rh   ( nz ) 
      REAL,    INTENT ( IN )  :: p    ( nz ) 
      REAL,    INTENT ( IN )  :: hgt  ( nz ) 
      REAL,    INTENT ( OUT ) :: cape        
      REAL,    INTENT ( OUT ) :: cin         
      REAL,    INTENT ( OUT ) :: zlcl        
      REAL,    INTENT ( OUT ) :: zlfc        
      REAL,    INTENT ( OUT ) :: zel         
      REAL,    INTENT ( OUT ) :: zmu         
      INTEGER                 :: ostat       
                                             

      INTEGER, INTENT ( IN  ) :: parcel      
                                             
                                             
  
      
      
      REAL                    :: ws   ( nz ) 
      REAL                    :: w    ( nz ) 
      REAL                    :: etheta( nz )
      REAL                    :: dTvK ( nz ) 
      REAL                    :: buoy ( nz ) 
      REAL                    :: tlclK       
      REAL                    :: plcl        
      REAL                    :: pel         
      REAL                    :: nbuoy       
      REAL                    :: pbuoy       
  
      
      
      REAL                    :: srctK       
      REAL                    :: srcrh       
      REAL                    :: srcws       
      REAL                    :: srcw        
      REAL                    :: srcp        
      REAL                    :: srctheta    
      REAL                    :: srcthetaeK  
      INTEGER                 :: srclev      
      REAL                    :: spdiff      
      REAL                    :: srce        
   
      
      
      REAL                    :: ptK        
      REAL                    :: ptvK       
      REAL                    :: tvK        
      REAL                    :: pw         
  
      
      
      INTEGER                 :: i, j, k    
      INTEGER                 :: lfclevhigh 
      INTEGER                 :: lfclevlow  
      INTEGER                 :: ellevhigh  
      INTEGER                 :: ellevlow   
      INTEGER                 :: prcl       
      INTEGER                 :: mlev       
      INTEGER                 :: lyrcnt     
      LOGICAL                 :: flag       
      LOGICAL                 :: wflag      
      REAL                    :: freeze     
      REAL                    :: pdiff      
      REAL                    :: hdiff      
      REAL                    :: pm, pu, pd 
      REAL                    :: lidxu      
      REAL                    :: lidxd      
  
      
      
      REAL                    :: Rd         
         PARAMETER ( Rd = 287.058 )         
      REAL                    :: Cp         
         PARAMETER ( Cp = 1004.67 )         
      REAL                    :: g          
         PARAMETER ( g  = 9.80665 )         
      REAL                    :: RUNDEF
         PARAMETER ( RUNDEF = -9.999E30 )
  
      
      
      ostat  = 0
      CAPE   = RUNDEF
      CIN    = RUNDEF
      ZLFC   = RUNDEF
      ZLCL   = RUNDEF
      ZEL    = RUNDEF
      ZMU    = RUNDEF
  
      
      
      
      
      
      IF ( parcel > 3 .or. parcel < 1 ) THEN
         prcl = 1
      ELSE
         prcl =  parcel
      END IF
  
      

      
      
      
      
      
      
      
  
      
      
      DO k = sfc, nz
        ws  ( k )   = SaturationMixingRatio ( tK(k), p(k) )
        w   ( k )   = ( rh(k)/100.0 ) * ws ( k )
        
        etheta(k)   = Thetae( tK(k), p(k)/100.0, rh(k), w(k) ) 
      END DO
  
      srclev      = sfc
      srctK       = tK    ( sfc )
      srcrh       = rh    ( sfc )
      srcp        = p     ( sfc )
      srcws       = ws    ( sfc )
      srcw        = w     ( sfc )
      srctheta    = Theta ( tK(sfc), p(sfc)/100.0 )
      srce        = etheta (sfc) 
   
      
      
      mlev = sfc + 1
      DO k = sfc, nz 
   
         
         
         pdiff = ( p (sfc) - p (k) ) / REAL ( 100 )
         IF ( pdiff <= REAL (100) ) mlev = k


         
         IF ( p(k) <= REAL (50000) ) EXIT

         IF ( prcl == 1 ) THEN
            
            IF (etheta(k) > srce)  THEN
               srctheta = Theta ( tK(k), p(k)/100.0 )
               srcw = w ( k )
               srclev  = k
               srctK   = tK ( k )
               srcrh   = rh ( k )
               srcp    = p  ( k )
               srce = etheta(k)
            END IF
         END IF
   
      END DO
   
      
      
      
      lyrcnt =  mlev - sfc + 1
      IF ( prcl == 2 ) THEN
   
         srclev   = sfc
         srctK    = SUM ( tK (sfc:mlev) ) / REAL ( lyrcnt )
         srcw     = SUM ( w  (sfc:mlev) ) / REAL ( lyrcnt )
         srcrh    = SUM ( rh (sfc:mlev) ) / REAL ( lyrcnt )
         srcp     = SUM ( p  (sfc:mlev) ) / REAL ( lyrcnt )
         srctheta = Theta ( srctK, srcp/100. )
   
      END IF
   
      srcthetaeK = Thetae ( srctK, srcp/100.0, srcrh, srcw )
      
      
      tlclk = TLCL ( tK(srclev), rh(srclev))
      ZLCL = hgt(srclev) + Cp/g * (tk(srclev)-tlclk)
      
      
      
   
      buoy  = REAL ( 0 )
      pw    = srcw
      wflag = .false.
      DO k  = srclev, nz
         IF ( hgt (k) >= ZLCL ) THEN
   


            
            

            
            
   
            IF ( wflag ) THEN
               flag  = .false.
   
               
               
               
               
               
               
               
               ptK   = The2T ( srcthetaeK, p(k), flag )
   
               
               
               
               
               
               pw = SaturationMixingRatio ( ptK, p(k) )
   
               
               
               
               ptvK  = VirtualTemperature ( ptK, pw )
               tvK   = VirtualTemperature ( tK (k), w (k) )
               
               
               
               
               
               
               
   
               
               
               
               
   
               
               dTvK ( k ) = ptvK - tvK
               buoy ( k ) = g * ( dTvK ( k ) / tvK )
   
            ELSE
   
               
               
               
               
               ptK   = srctheta / ( 100000.0/p(k) )**(Rd/Cp)
   
               
               
               
               ptvK  = VirtualTemperature ( ptK, srcw )
   
               
               
               tvK   = VirtualTemperature ( tK (k), w (k) )
   
               
               
               dTvK ( k ) = ptvK - tvK
               buoy ( k ) = g * ( dtvK ( k ) / tvK )
   
               wflag = .true.
   
            END IF
   
         ELSE
   
            
            
            
            
            ptK   = srctheta / ( 100000.0/p(k) )**(Rd/Cp)
   
            
            
            
            ptvK  = VirtualTemperature ( ptK, srcw )
   
            
            
            tvK   = VirtualTemperature ( tK (k), w (k) )
   
            
            
            dTvK ( k ) = ptvK - tvK
            buoy ( k ) = g * ( dtvK ( k ) / tvK )
   
         END IF

         
         

   
      END DO
      
      
      
      flag = .false. 
      lfclevhigh = -1
      lfclevlow = -1
      ellevhigh = -1
      ellevlow = -1
      DO k = sfc, nz 
         
         
         IF (buoy (k) > REAL (0) .and. hgt(k) >= ZLCL .and. (.not.flag)) THEN 
            lfclevhigh = k
            lfclevlow = k - 1
            flag = .true. 
         END IF
         
         
         IF (k>1) THEN
            IF (flag .and. buoy (k) <= REAL (0) .and. buoy (k-1) > REAL (0)) THEN
                ellevhigh = k
                ellevlow = k - 1
            END IF
         END IF
         IF (buoy (k) < REAL (0) .and. flag) THEN 
            flag = .false. 
         END IF
      END DO
      
      
      IF ( lfclevhigh > 0 ) THEN
         
         IF (lfclevlow == 0) THEN
             ZLFC = hgt(lfclevhigh)
         ELSE IF (hgt(lfclevlow)<= ZLCL .and. buoy(lfclevlow) >= REAL(0)) THEN
             ZLFC = ZLCL 
         ELSE IF (buoy(lfclevlow) == REAL(0) .and. buoy(lfclevhigh) > REAL(0)) THEN
             ZLFC = hgt(lfclevlow) 
         ELSE IF (lfclevhigh > lfclevlow) THEN
             ZLFC = hgt(lfclevlow) + (hgt(lfclevhigh)-hgt(lfclevlow))/(buoy(lfclevhigh)-buoy(lfclevlow)) * (REAL(0)-buoy(lfclevlow))
         END IF
      END IF
      
      
      IF ( ellevhigh > 0 ) THEN
         
         IF (buoy(ellevhigh) == REAL(0)) THEN
             ZEL = hgt(ellevhigh) 
         ELSE IF (ellevhigh > ellevlow) THEN
             ZEL = hgt(ellevlow) + (hgt(ellevhigh)-hgt(ellevlow))/(buoy(ellevhigh)-buoy(ellevlow)) * (REAL(0)-buoy(ellevlow))
         END IF
      END IF
      
      
      IF (ellevhigh > 0) THEN 
         CAPE=REAL ( 0 ) 
         CIN=REAL ( 0 )
         DO k = sfc+1, nz
            
            IF ( hgt (k) > ZLFC .and. hgt (k) < ZEL) THEN
               CAPE = CAPE + MAX ( buoy (k), 0.0 ) * ( hgt (k) - hgt (k-1) )
            END IF
            
            IF ( hgt (k) > ZLCL .and. hgt (k) < ZLFC) THEN
               CIN  = CIN  + MIN ( buoy (k), 0.0 ) * ( hgt (k) - hgt (k-1) )
            END IF
         END DO
      END IF

      ZMU = hgt ( srclev )

      IF ( CAPE /= CAPE ) cape = RUNDEF
      IF ( CIN  /= CIN  ) cin  = RUNDEF
      
      






   
  END FUNCTION Buoyancy 












































































  FUNCTION MSLP ( zsfc, psfc, zlev1, qlev1, tlev1 )

      implicit none
     
     


      REAL,    INTENT ( IN )  :: zsfc         
      REAL,    INTENT ( IN )  :: psfc         
      REAL,    INTENT ( IN )  :: zlev1        
      REAL,    INTENT ( IN )  :: qlev1        
      REAL,    INTENT ( IN )  :: tlev1        
      real,PARAMETER :: G=9.81
      real,PARAMETER :: GI=1./G
      real,PARAMETER :: RD=287.0
      real,PARAMETER :: ZSL=0.0
      real,PARAMETER :: TAUCR=RD*GI*290.66,CONST=0.005*G/RD
      real,PARAMETER :: GORD=G/RD,DP=60.E2
      real,PARAMETER :: GAMMA=6.5E-3

      real MSLP,TVRT,TVRSFC,TAUSFC,TVRSL,TAUSL,TAUAVG




         MSLP = PSFC


         TVRT = TLEV1*(1.0+0.608*QLEV1)
         



         TVRSFC = TVRT + (ZLEV1 - ZSFC)*GAMMA
         TAUSFC = TVRSFC*RD*GI
         TVRSL  = TVRT + (ZLEV1 - ZSL)*GAMMA
         TAUSL  = TVRSL*RD*GI


         IF ((TAUSL.GT.TAUCR).AND.(TAUSFC.LE.TAUCR)) THEN
            TAUSL=TAUCR
         ELSEIF ((TAUSL.GT.TAUCR).AND.(TAUSFC.GT.TAUCR)) THEN
            TAUSL = TAUCR-CONST*(TAUSFC-TAUCR)**2
         ENDIF


         TAUAVG = 0.5*(TAUSL+TAUSFC)


         MSLP = PSFC*EXP(ZSFC/TAUAVG)

  END FUNCTION MSLP



  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  FUNCTION calc_fits ( p, tK, rh ) RESULT ( fits )
 
    implicit none

    
    
    real, intent ( in ) :: p               
    real, intent ( in ) :: tK              
    real, intent ( in ) :: rh              
    real                :: fits            
 
    
    
    real                :: twb             
    real                :: wbt
 
    
    

    
    
    fits = REAL ( 0 )

    
    
    twb =  WetBulbTemp ( p, tK, rh ) - 273.15 

    
    
    fits = 0.8281*twb + 0.3549*( tK - 273.15 ) + 5.08
 
    
    
    fits = fits + 273.15

  END FUNCTION calc_fits



  
  
  
  
  
  
  
  
  
  
  FUNCTION calc_wc ( tK, wspd ) RESULT ( wc )

    implicit none

    
    
    real, intent ( in  ) :: tK
    real, intent ( in  ) :: wspd

    real                 :: tF, wc, wspd_mph

    wspd_mph = wspd * 2.23693629 
    tF  = (( tK - 273.15 ) * ( REAL (9) / REAL (5) ) ) + REAL ( 32 )

    wc =    35.74                           &
       +  (  0.6215 * tF                  ) &
       -  ( 35.75   *      ( wspd_mph**0.16 ) ) &
       +  (  0.4275 * tF * ( wspd_mph**0.16 ) )

    wc = (( wc - REAL (32) ) * ( REAL (5) / REAL (9) ) ) + 273.15

  END FUNCTION calc_wc



  
  
  
  
  
  
  
  
  
  
  FUNCTION calc_hi ( Tk, RH ) result ( HI )

    implicit none

    
    
    real, intent ( in  ) :: Tk
    real, intent ( in  ) :: RH

    real :: tF, tF2, rh2, HI

    
    
    
    IF ( Tk > 294.26111 ) THEN

      tF   = ( (Tk - 273.15) * (REAL (9)/REAL (5))  ) + REAL ( 32 )
      tF2  = tF ** 2
      rh2  = RH ** 2

      HI =  -42.379 &
         +  (  2.04901523   * tF              ) &
         +  ( 10.14333127   * RH              ) &
         -  (  0.22475541   * tF  * RH        ) &
         -  (  6.83783E-03  * tF2             ) &
         -  (  5.481717E-02 * rh2             ) &
         +  (  1.22874E-03  * tF2 * RH        ) &
         +  (  8.5282E-04   * tF  * rh2       ) &
         -  (  1.99E-06     * tF2 * rh2       )

      HI = ((HI - REAL (32)) * (REAL (5)/REAL (9))) + 273.15
    ELSE
      HI = Tk
    END IF

  END FUNCTION calc_hi

  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  FUNCTION WetBulbTemp ( p, tK, rh) result( wbt )

    implicit none

    
    
    real, intent ( in ) :: p        
    real, intent ( in ) :: tK       
    real, intent ( in ) :: rh       
    real                :: wbt      
 
    
    
    real                :: tdK      
    real                :: tC       
    real                :: tdC      
    real                :: svapr    
    real                :: vapr     
    real                :: gamma    
    real                :: delta    

    
    
    
    
    wbt = REAL ( 0 )
    tC  = tK - 273.15

    
    
    svapr = calc_es ( tK ) * REAL ( 100 )

    
    
    vapr  = svapr * ( rh / REAL (100) )

    
    
    tdC = calc_Dewpoint ( tC, rh )
    tdK = tdC + 273.15

    
    
    gamma = 0.00066 * ( p / REAL (1000) )
    delta = REAL ( 4098 ) * ( vapr / REAL(1000) )  / ( (tC+237.3)**2 )

    
    
    wbt = ( ((gamma * tC) + (delta * tdC)) / (gamma + delta) ) + 273.15

  END FUNCTION WetBulbTemp


  
  
  
  
  
  
  
  
  
  FUNCTION calc_Dewpoint ( tC, rh) result( Dewpoint )

    implicit none

    
    
    real, intent ( in ) :: tC
    real, intent ( in ) :: rh
    real                :: Dewpoint
 
    real :: term, es, e1, e, logs, expon

    expon    = ( 7.5*tC ) / ( 237.7+tC )
    es       = 6.112 * ( 10**expon )     
    e        = es * ( rh/100.0 )         
    logs     = LOG10 ( e/6.112 )
    Dewpoint = ( 237.7*logs ) / ( 7.5-logs )

  END FUNCTION calc_Dewpoint


  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  FUNCTION calc_es ( tK ) result ( es )

    implicit none

    
    
    real, intent ( in ) :: tK
    real                :: es
 
    real                :: p1, p2, c1

    p1 = 11.344    - 0.0303998 * tK
    p2 = 3.49149   - 1302.8844 / tK
    c1 = 23.832241 - 5.02808   * ALOG10 ( tK )
    es = 10.**(c1-1.3816E-7*10.**p1+8.1328E-3*10.**p2-2949.076/tK)

  END FUNCTION calc_es



  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  FUNCTION CATTurbulence ( ugrdbot, ugrdtop, vgrdbot, vgrdtop &
                           ,defor11bot, defor11top, defor12bot, defor12top &
                           ,defor22bot, defor22top, zbot, ztop ) result ( ti )

    IMPLICIT NONE

    
    
    REAL,    INTENT ( IN )  :: ugrdbot       
    REAL,    INTENT ( IN )  :: ugrdtop       
    REAL,    INTENT ( IN )  :: vgrdbot       
    REAL,    INTENT ( IN )  :: vgrdtop       
    REAL,    INTENT ( IN )  :: defor11bot    
    REAL,    INTENT ( IN )  :: defor11top    
    REAL,    INTENT ( IN )  :: defor12bot    
    REAL,    INTENT ( IN )  :: defor12top    
    REAL,    INTENT ( IN )  :: defor22bot    
    REAL,    INTENT ( IN )  :: defor22top    
    REAL,    INTENT ( IN )  :: zbot          
    REAL,    INTENT ( IN )  :: ztop          
    REAL                    :: ti            

    
    
    REAL    :: dudx, dudx1, dudx2 
    REAL    :: dvdy, dvdy1, dvdy2
    REAL    :: dudz, dvdz

    REAL    :: depth, vws, conv    
    REAL    :: def, shear, stretch 

    
    
    ti = REAL ( 0 )

    
    
    depth = ABS ( ztop - zbot )
    dudz  = ( ugrdbot - ugrdtop ) / depth
    dvdz  = ( vgrdbot - vgrdtop ) / depth
    vws   = SQRT ( dudz**2 + dvdz**2  )

    dudx1 = defor11top / 2.
    dudx2 = defor11bot / 2.
    dudx  = ( dudx1 + dudx2 ) / REAL ( 2 )

    dvdy1 = defor22top / 2.
    dvdy2 = defor22bot / 2.
    dvdy  = ( dvdy1 + dvdy2 ) / REAL ( 2 )

    
    
    stretch = dudx - dvdy
    shear   = ( defor12top + defor12bot ) / REAL ( 2 )
    def     = SQRT ( stretch**2 + shear**2 )

    
    
    conv    = - ( dudx + dvdy )

    
    
    ti = vws * ( def + conv ) * 1.0E+07

    IF ( ti /= ti ) ti = REAL ( 0 )
    IF ( ti < 0   ) ti = REAL ( 0 )

  END FUNCTION CATTurbulence



  FUNCTION lin_interp ( x, f, y ) result ( g )

    
    
    
    
    
    

    
    

    implicit none

    real, intent(in), dimension(:) :: x  
    real, intent(in), dimension(:) :: f  
    real, intent(in) :: y                
    real :: g                            

    integer :: k  
    integer :: n  
    real    :: a

    n = size(x)

    
    

    if (y <= x(1)) then
      k = 1
    else if (y >= x(n)) then
      k = n - 1
    else
      k = 1
      do while (y > x(k+1) .and. k < n)
        k = k + 1
      end do
    end if
    
    a = (  f(k+1) - f(k) ) / ( x(k+1) - x(k) )
    g = f(k) + a * (y - x(k))

  END FUNCTION lin_interp



  
  
  
  
  
  
  
  
  
  
  FUNCTION LLT_Windspeed ( nlayer, u, v ) RESULT ( dynamic )
    IMPLICIT NONE
 
    
    
    INTEGER, INTENT ( IN )         :: nlayer
    REAL, INTENT ( IN )            :: u     ( nlayer )
    REAL, INTENT ( IN )            :: v     ( nlayer )
    REAL                           :: dynamic
 
    
    
    INTEGER           :: i
    REAL              :: this_windspeed ( nlayer )
    REAL              :: PI
       PARAMETER ( PI = 3.14159265359 )

    
    
    
    
    dynamic = REAL ( 0 )

    
    
    DO i = 1, nlayer
       this_windspeed ( i ) = SQRT ( u(i)**2 + v(i)**2 )
    END DO

    
    
    dynamic = ( this_windspeed(1)+this_windspeed(nlayer) ) / REAL (20)
    IF ( dynamic > REAL (2) ) dynamic = REAL ( 2 )
    dynamic = ( dynamic + REAL (3) ) / REAL ( 2 )
    dynamic = SIN ( dynamic*PI )
    dynamic = ( dynamic + REAL (1) ) / REAL ( 2 )


  END FUNCTION LLT_Windspeed


  
  
  
  
  
  
  
  
  
  
  FUNCTION LLT_Thermodynamic ( nlayer, tK, hgt ) RESULT ( thermodynamic )
  IMPLICIT NONE
 
    
    
    INTEGER, INTENT ( IN )         :: nlayer
    REAL, INTENT ( IN )            :: tK     ( nlayer ) 
    REAL, INTENT ( IN )            :: hgt    ( nlayer ) 
    REAL                           :: thermodynamic

    
    
    INTEGER :: i
    REAL    :: lapse
    REAL    :: PI
       PARAMETER ( PI = 3.14159265359 )

    
    

    
    
    thermodynamic = REAL ( 0 )

    
    
    
    lapse = ( tk(1) - tk(nlayer) ) * REAL ( 1000 )
    lapse = lapse / ( hgt(nlayer) - hgt(1) )

    
    
    thermodynamic = lapse / REAL ( 10 )
    thermodynamic = ( thermodynamic + REAL (3) ) / REAL ( 2 )
    thermodynamic = SIN ( thermodynamic * PI )
    thermodynamic = ( thermodynamic + REAL (1) ) / REAL ( 2 )

  END FUNCTION LLT_Thermodynamic


  
  
  
  
  
  
  
  
  
  
  FUNCTION LLT_MountainWave ( nlayer, tdx, tdy, u, v, tK, hgt) &
                                                         RESULT ( MountainWave )
    IMPLICIT NONE

    
    
    INTEGER, INTENT ( IN )           :: nlayer
    REAL, INTENT ( IN )              :: tdx               
    REAL, INTENT ( IN )              :: tdy               
    REAL, INTENT ( IN )              :: u   ( nlayer )    
    REAL, INTENT ( IN )              :: v   ( nlayer )    
    REAL, INTENT ( IN )              :: tK  ( nlayer )    
    REAL, INTENT ( IN )              :: hgt ( nlayer )    
    REAL                             :: MountainWave      
 
    
    
    REAL    :: u_term
    REAL    :: v_term
    REAL    :: uv_term
    REAL    :: lapse
    REAL    :: total_mw, this_total_mw
    REAL    :: this_uv_term
    REAL    :: min_uv_term, cross_terrain, max_total_mw
    INTEGER :: i, j, k

    REAL    :: PI
       PARAMETER ( PI = 3.14159265359 )

    
    

    
    
    MountainWave = REAL ( 0 )

    
    
    DO i = 2, nlayer

      
      
      u_term       = ( (u(i-1) + u(i) ) / REAL(2) ) * tdx
      v_term       = ( (v(i-1) + v(i) ) / REAL(2) ) * tdy
      this_uv_term = ( u_term + v_term ) * REAL ( -1 )
      
      IF ( min_uv_term < REAL (0) ) min_uv_term = REAL ( 0 )
      IF ( i == 2 ) THEN
        
        min_uv_term = this_uv_term
      ELSE
        
        IF ( this_uv_term < min_uv_term ) min_uv_term = this_uv_term
      END IF

      
      
      lapse = ( tK (i-1) - tK (i) ) * REAL ( 1000 )
      lapse  = lapse / ABS ( hgt(i)-hgt(i-1) )
      IF ( lapse > REAL (0) ) lapse = REAL ( 0 )
      lapse = lapse * REAL ( -1 )

      this_total_mw = this_uv_term * lapse * REAL ( 40000 )
      IF ( i == 2 ) THEN
        total_mw = this_total_mw
      ELSE
        IF ( this_total_mw > total_mw ) total_mw = this_total_mw
      END IF
 
    END DO

    
    cross_terrain = min_uv_term * REAL ( 500 )

    IF ( min_uv_term < 0.03 ) THEN
      cross_terrain = REAL ( 0 )
    END IF

    IF ( cross_terrain > REAL (50) ) cross_terrain = REAL ( 50 )

    
    
    IF ( total_mw > REAL (50) ) total_mw = REAL ( 50 )

    
    
    MountainWave = ( total_mw*(cross_terrain/50.) ) + cross_terrain
    MountainWave = MountainWave / REAL ( 100 )

  END FUNCTION LLT_MountainWave



  
  
  
  
  
  
  
  
  
  
  FUNCTION LLT_TrappedWave ( nlayer, u, v, p ) RESULT ( trapped )
     IMPLICIT NONE

     
     
     INTEGER, INTENT ( IN )         :: nlayer
     REAL, INTENT ( IN )            :: u     ( nlayer )
     REAL, INTENT ( IN )            :: v     ( nlayer )
     REAL, INTENT ( IN )            :: p     ( nlayer )
     REAL                           :: trapped

     
     
     INTEGER           :: i
     REAL              :: du, dv
     REAL              :: scale_fact, this_p
     REAL              :: dudv, this_dudv
     REAL              :: PI
        PARAMETER ( PI = 3.14159265359 )

     
     
     REAL, PARAMETER :: scale_950 = 0.050000  
     REAL, PARAMETER :: scale_925 = 0.040000  
     REAL, PARAMETER :: scale_900 = 0.025000  
     REAL, PARAMETER :: scale_850 = 0.010000  
     REAL, PARAMETER :: scale_800 = 0.005000  
     REAL, PARAMETER :: scale_750 = 0.002941  
     REAL, PARAMETER :: scale_700 = 0.001923  
     REAL, PARAMETER :: scale_650 = 0.001351  
     REAL, PARAMETER :: scale_600 = 0.001000  
     REAL, PARAMETER :: scale_550 = 0.000800  

     
     

     
     
     trapped = REAL ( 0 )

     
     
     dudv = REAL ( 0 )
     DO i = 2, nlayer

       
       
       du         = u ( i-1 ) - u ( i )
       dv         = v ( i-1 ) - v ( i )

       
       
       this_p = p ( i ) / REAL ( 100 )
       IF ( this_p > REAL (950) ) THEN
         scale_fact = scale_950
       ELSE IF ( this_p <= REAL (950) .AND. this_p > REAL (925) ) THEN
         scale_fact = scale_925
       ELSE IF ( this_p <= REAL (925) .AND. this_p > REAL (900) ) THEN
         scale_fact = scale_900
       ELSE IF ( this_p <= REAL (900) .AND. this_p > REAL (850) ) THEN
         scale_fact = scale_850
       ELSE IF ( this_p <= REAL (850) .AND. this_p > REAL (800) ) THEN
         scale_fact = scale_800
       ELSE IF ( this_p <= REAL (800) .AND. this_p > REAL (750) ) THEN
         scale_fact = scale_750
       ELSE IF ( this_p <= REAL (750) .AND. this_p > REAL (700) ) THEN
         scale_fact = scale_700
       ELSE IF ( this_p <= REAL (700) .AND. this_p > REAL (650) ) THEN
         scale_fact = scale_650
       ELSE IF ( this_p <= REAL (650) .AND. this_p > REAL (600) ) THEN
         scale_fact = scale_600
       ELSE IF ( this_p <= REAL (600) ) THEN
         scale_fact = scale_550
       END IF

       this_dudv = ( (du**2)*(dv**2) ) * scale_fact
       IF ( this_dudv > dudv ) dudv = this_dudv

    END DO

    trapped = dudv
    IF ( trapped > REAL ( 1 ) ) trapped = REAL ( 1 )
    trapped = trapped / REAL ( 4 )

  END FUNCTION LLT_TrappedWave
  
    FUNCTION diag_map ( nz, nj, ni, tk, rh, p, hgt, sfc, cape, cin, zlcl, zlfc, zel, zmu,  &
                        parcel ) result (ostat) 
  
      IMPLICIT NONE
  
      INTEGER, INTENT ( IN )  :: nz, nj, ni          
      INTEGER, INTENT ( IN )  :: sfc         
      REAL,    INTENT ( IN )  :: tk   ( nz, nj, ni ) 
      REAL,    INTENT ( IN )  :: rh   ( nz, nj, ni ) 
      REAL,    INTENT ( IN )  :: p    ( nz, nj, ni )   
      REAL,    INTENT ( IN )  :: hgt  ( nz, nj, ni ) 
      REAL,    INTENT ( OUT ) :: cape ( nj, ni )      
      REAL,    INTENT ( OUT ) :: cin  ( nj, ni )       
      REAL,    INTENT ( OUT ) :: zlcl ( nj, ni )       
      REAL,    INTENT ( OUT ) :: zlfc ( nj, ni )        
      REAL,    INTENT ( OUT ) :: zel  ( nj, ni )         
      REAL,    INTENT ( OUT ) :: zmu  ( nj, ni )       
      INTEGER                 :: ostat                                              
      INTEGER, INTENT ( IN  ) :: parcel
      
      INTEGER                 :: i, j 
                                                  
      ostat = 0
      DO i = 1, ni
         DO j = 1, nj
             ostat = Buoyancy (                         nz &
                                     ,        tk(1:nz,j,i) &
                                     ,        rh(1:nz,j,i) &
                                     ,         p(1:nz,j,i) &
                                     ,       hgt(1:nz,j,i) &
                                     ,                 sfc & 
                                     ,           cape(j,i) &
                                     ,            cin(j,i) &
                                     ,           zlcl(j,i) &
                                     ,           zlfc(j,i) &
                                     ,            zel(j,i) &
                                     ,           zmu(j,i)  &
                                     ,              parcel ) 
         END DO
      END DO
      ostat = 1
      
  END FUNCTION diag_map
  
  FUNCTION diag_row ( nz, ni, tk, rh, p, hgt, sfc, cape, cin, zlcl, zlfc, zel, zmu,  &
                        parcel ) result (ostat) 
  
      IMPLICIT NONE
  
      INTEGER, INTENT ( IN )  :: nz, ni          
      INTEGER, INTENT ( IN )  :: sfc         
      REAL,    INTENT ( IN )  :: tk   ( nz, ni ) 
      REAL,    INTENT ( IN )  :: rh   ( nz, ni ) 
      REAL,    INTENT ( IN )  :: p    ( nz, ni )   
      REAL,    INTENT ( IN )  :: hgt  ( nz, ni ) 
      REAL,    INTENT ( OUT ) :: cape ( ni )      
      REAL,    INTENT ( OUT ) :: cin  ( ni )       
      REAL,    INTENT ( OUT ) :: zlcl ( ni )       
      REAL,    INTENT ( OUT ) :: zlfc ( ni )        
      REAL,    INTENT ( OUT ) :: zel  ( ni )         
      REAL,    INTENT ( OUT ) :: zmu  ( ni )       
      INTEGER                 :: ostat                                              
      INTEGER, INTENT ( IN  ) :: parcel
      
      INTEGER                 :: i, j
                                                  
      ostat = 0
      DO i = 1, ni
         ostat = Buoyancy (                        nz &
                                 ,        tk(1:nz, i) &
                                 ,        rh(1:nz, i) &
                                 ,         p(1:nz, i) &
                                 ,       hgt(1:nz, i) &
                                 ,                sfc & 
                                 ,            cape(i) &
                                 ,             cin(i) &
                                 ,           zlcl(i)  &
                                 ,            zlfc(i) &
                                 ,           zel(i)   &
                                 ,            zmu(i)  &
                                 ,             parcel ) 
      END DO
      ostat = 1
  END FUNCTION diag_row

END MODULE diag_functions
