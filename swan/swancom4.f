!
!     SWAN/COMPU   file 4 of 5
!
!
!     PROGRAM SWANCOM4.FOR
!
!
!     This file SWANCOM4 of the main program SWAN
!     include the next subroutines
!
!     *** nonlinear 4 wave-wave interactions ***
!
!     FAC4WW (compute the constants for the nonlinear wave
!             interactions)
!     RANGE4 (compute the counters for the different types of
!             computations for the nonlinear wave interactions)
!     SWSNL1 (nonlinear four wave interactions; semi-implicit and computed
!             for all bins that fall within a sweep with DIA technique.
!             Interaction are calculated per sweep)
!     SWSNL2 (nonlinear four wave interactions; fully explicit and computed
!             for all bins that fall within a sweep with DIA technique.
!             Interaction are calculated per sweep)
!     SWSNL3 (calculate nonlinear four wave interactions fully explicitly
!             for the full circle per iteration by means of DIA approach
!             and store results in auxiliary array MEMNL4)
!     SWSNL4 (calculate nonlinear four wave interactions fully explicitly
!             for the full circle per iteration by means of MDIA approach
!             and store results in auxiliary array MEMNL4)
!     SWSNL8 (calculate nonlinear four wave interactions fully explicitly
!             for the full circle per iteration by means of DIA approach
!             and store results in auxiliary array MEMNL4. Neighbouring
!             interactions are interpolated in piecewise constant manner)
!     FILNL3 (fill main diagonal and right-hand side of the system with
!             results of array MEMNL4)
!
!     RIAM_SLW (calculate nonlinear four wave interactions by means
!               of the exact FD-RIAM technique)
!
!     SWINTFXNL (interface with SWAN model to compute nonlinear transfer
!                with the XNL method for given action density spectrum)
!
!     *** nonlinear 3 wave-wave interactions ***
!
!     SWLTA  (triad-wave interactions calculated with the Lumped Triad
!             Approximation of Eldeberky, 1996)
!
!----------------------------------------------------------------------
!
!******************************************************************
!
      SUBROUTINE FAC4WW (XIS   ,SNLC1 ,                                   40.41 34.00
     &                  DAL1  ,DAL2  ,DAL3         ,SPCSIG,               34.00
     &                  WWINT ,WWAWG ,WWSWG                )              40.17 34.00
!
!******************************************************************
!
      USE SWCOMM3                                                         40.41
      USE SWCOMM4                                                         40.41
      USE OCPCOMM4                                                        40.41
      USE M_SNL4                                                          40.17
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering                              |
!     | Fluid Mechanics Section                                   |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmers: H.L. Tolman, R.C. Ris                        |
!   --|-----------------------------------------------------------|--
!
!
!     SWAN (Simulating WAves Nearshore); a third generation wave model
!     Copyright (C) 2008  Delft University of Technology
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!     GNU General Public License for more details.
!
!     A copy of the GNU General Public License is available at
!     http://www.gnu.org/copyleft/gpl.html#SEC3
!     or by writing to the Free Software Foundation, Inc.,
!     59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!
!  0. Authors
!
!     30.72: IJsbrand Haagsma
!     40.17: IJsbrand Haagsma
!     40.41: Marcel Zijlema
!
!  1. Updates
!
!     30.72, Feb. 98: Introduced generic names XCGRID, YCGRID and SPCSIG for SWAN
!     40.17, Dec. 01: Implementation of Multiple DIA
!     40.41, Sep. 04: compute indices for interactions which will be
!                     interpolated in piecewise constant manner
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. Purpose :
!
!     Calculate interpolation constants for Snl.
!
!  3. Method :
!
!
!  4. Argument variables
!
!     SPCSIG: Relative frequencies in computational domain in sigma-space 30.72
!
      REAL    SPCSIG(MSC)                                                 30.72
!
!     INTEGERS:
!     ---------
!     MSC2,MSC1         Auxiliary variables
!     MSC,MDC           Maximum counters in spectral space
!     IDP,IDP1          Positive range for ID
!     IDM,IDM1          Negative range for ID
!     ISP,ISP1          idem for IS
!     ISM,ISM1          idem for IS
!     ISCLW,ISCHG       Minimum and maximum counter for discrete
!                       computations in frequency space
!     ISLOW,ISHGH       Minimum and maximum range in frequency space
!     IDLOW,IDHGH       idem in directional space
!     IS                Frequency counter
!     MSC4MI,MSC4MA     Array dimensions in frequency space
!     MDC4MI,MDC4MA     Array dimensions in direction space
!
!     REALS:
!     ------
!     LAMBDA            Coefficient set 0.25
!     GRAV              Gravitational acceleration
!     SNLC1             Coefficient for the subroutines SWSNLn
!     LAMM2,LAMP2
!     DELTH3,DELTH4     Angles between the interacting wavenumbers
!     DAL1,DAL2,DAL3    Coefficients for the non linear interactions
!     CIDP,CIDM
!     WIDP,WIDP1,WIDM,WIDM1  Weight factors
!     WISP,WISP1,WISM,WISM1  idem
!     AWGn              Interpolation weight factors
!     SWGn              Quadratic interpolation weight factors
!     XIS,XISLN         Difference between succeeding frequencies
!     PI                3.14
!     FREQ              Auxiliary frequency to fill scaling array
!     DDIR,RADE         band width in directional space and factor        34.00
!
!     ARRAYS
!     ------
!     AF11    1D   Scaling frequency
!     WWINT   1D   counters for 4WAVE interactions
!     WWAWG   1D   values for the interpolation
!     WWSWG   1D   vaules for the interpolation
!
!     WWINT ( 1 = IDP    WWAWG ( = AGW1    WWSWG ( = SWG1
!             2 = IDP1           = AWG2            = SWG2
!             3 = IDM            = AWG3            = SWG3
!             4 = IDM1           = AWG4            = SWG4
!             5 = ISP            = AWG5            = SWG5
!             6 = ISP1           = AWG6            = SWG6
!             7 = ISM            = AWG7            = SWG7
!             8 = ISM1           = AWG8 )          = SWG8  )
!             9 = ISLOW
!             10= ISHGH
!             11= ISCLW
!             12= ISCHG
!             13= IDLOW
!             14= IDHGH
!             15= MSC4MI
!             16= MSC4MA
!             17= MDC4MI
!             18= MDC4MA
!             19= MSCMAX
!             20= MDCMAX
!             21= IDPP
!             22= IDMM
!             23= ISPP
!             24= ISMM )
!
!  7. Common blocks used
!
!
!  9. Source code :
!
!     -----------------------------------------------------------------
!     Calculate :
!       1. counters for frequency and direction for NL-interaction
!       2. weight factors
!       3. the minimum and maximum counter in IS and ID space
!       4. the interpolation weights
!       5. the quadratic interpolation rates
!       6. fill the array for the frequency**11
!     ----------------------------------------------------------
!
!****************************************************************
!
      INTEGER     MSC2  ,MSC1  ,IS    ,IDP   ,IDP1  ,                     40.41 34.00
     &            IDM   ,IDM1  ,ISP   ,ISP1  ,ISM   ,ISM1  ,              34.00
     &            IDPP  ,IDMM  ,ISPP  ,ISMM  ,                            40.41
     &            ISLOW ,ISHGH ,ISCLW ,ISCHG ,IDLOW ,IDHGH ,
     &            MSCMAX,MDCMAX                                           34.00
!
      REAL        SNLC1 ,LAMM2 ,LAMP2 ,DELTH3,                            40.17 34.00
     &            AUX1  ,DELTH4,DAL1  ,DAL2  ,DAL3  ,CIDP  ,WIDP  ,
     &            WIDP1 ,CIDM  ,WIDM  ,WIDM1 ,XIS   ,XISLN ,WISP  ,
     &            WISP1 ,WISM  ,WISM1 ,AWG1  ,AWG2  ,AWG3  ,AWG4  ,
     &            AWG5  ,AWG6  ,AWG7  ,AWG8  ,SWG1  ,SWG2  ,SWG3  ,
     &            SWG4  ,SWG5  ,SWG6  ,SWG7  ,SWG8  ,FREQ  ,              34.00
     &            RADE                                                    34.00
!
      REAL       WWAWG(*)               ,                                 40.17
     &           WWSWG(*)
!
      INTEGER    WWINT(*)
!
      SAVE IENT
      DATA IENT/0/

      IF (LTRACE) CALL STRACE (IENT,'FAC4WW')

      IF (ALLOCATED(AF11)) DEALLOCATE(AF11)                               40.17

!     *** Compute frequency indices                               ***
!     *** XIS is the relative increment of the relative frequency ***
!
      MSC2   = INT ( FLOAT(MSC) / 2.0 )
      MSC1   = MSC2 - 1
      XIS    = SPCSIG(MSC2) / SPCSIG(MSC1)                                30.72
!
!     *** set values for the nonlinear four-wave interactions ***
!
      SNLC1  = 1. / GRAV**4                                               40.17 34.00
!
      LAMM2  = (1.-PQUAD(1))**2                                           40.17
      LAMP2  = (1.+PQUAD(1))**2                                           40.17
      DELTH3 = ACOS( (LAMM2**2+4.-LAMP2**2) / (4.*LAMM2) )
      AUX1   = SIN(DELTH3)
      DELTH4 = ASIN(-AUX1*LAMM2/LAMP2)
!
      DAL1   = 1. / (1.+PQUAD(1))**4                                      40.17
      DAL2   = 1. / (1.-PQUAD(1))**4                                      40.17
      DAL3   = 2. * DAL1 * DAL2
!
!     *** Compute directional indices in sigma and theta space ***
!
      CIDP   = ABS(DELTH4/DDIR)                                           40.00
      IDP   = INT(CIDP)
      IDP1  = IDP + 1
      WIDP   = CIDP - REAL(IDP)
      WIDP1  = 1.- WIDP
!
      CIDM   = ABS(DELTH3/DDIR)                                           40.00
      IDM   = INT(CIDM)
      IDM1  = IDM + 1
      WIDM   = CIDM - REAL(IDM)
      WIDM1  = 1.- WIDM
      XISLN  = LOG( XIS )
!
      ISP    = INT( LOG(1.+PQUAD(1)) / XISLN )                            40.17
      ISP1   = ISP + 1
      WISP   = (1.+PQUAD(1) - XIS**ISP) / (XIS**ISP1 - XIS**ISP)          40.17
      WISP1  = 1. - WISP
!
      ISM    = INT( LOG(1.-PQUAD(1)) / XISLN )                            40.17
      ISM1   = ISM - 1
      WISM   = (XIS**ISM -(1.-PQUAD(1))) / (XIS**ISM - XIS**ISM1)         40.17
      WISM1  = 1. - WISM
!
!     *** Range of calculations ***
!
      ISLOW =  1  + ISM1
      ISHGH = MSC + ISP1 - ISM1
      ISCLW =  1
      ISCHG = MSC - ISM1
      IDLOW = 1 - MDC - MAX(IDM1,IDP1)
      IDHGH = MDC + MDC + MAX(IDM1,IDP1)
!
      MSC4MI = ISLOW
      MSC4MA = ISHGH
      MDC4MI = IDLOW
      MDC4MA = IDHGH
      MSCMAX = MSC4MA - MSC4MI + 1
      MDCMAX = MDC4MA - MDC4MI + 1
!
!     *** Interpolation weights ***
!
      AWG1   = WIDP  * WISP
      AWG2   = WIDP1 * WISP
      AWG3   = WIDP  * WISP1
      AWG4   = WIDP1 * WISP1
!
      AWG5   = WIDM  * WISM
      AWG6   = WIDM1 * WISM
      AWG7   = WIDM  * WISM1
      AWG8   = WIDM1 * WISM1
!
!     *** quadratic interpolation ***
!
      SWG1   = AWG1**2
      SWG2   = AWG2**2
      SWG3   = AWG3**2
      SWG4   = AWG4**2
!
      SWG5   = AWG5**2
      SWG6   = AWG6**2
      SWG7   = AWG7**2
      SWG8   = AWG8**2
!
!     --- determine discrete counters for piecewise                       40.41
!         constant interpolation                                          40.41
!
      IF (AWG1.LT.AWG2) THEN
         IF (AWG2.LT.AWG3) THEN
            IF (AWG3.LT.AWG4) THEN
               ISPP=ISP
               IDPP=IDP
            ELSE
               ISPP=ISP
               IDPP=IDP1
            END IF
         ELSE IF (AWG2.LT.AWG4) THEN
            ISPP=ISP
            IDPP=IDP
         ELSE
            ISPP=ISP1
            IDPP=IDP
         END IF
      ELSE IF (AWG1.LT.AWG3) THEN
         IF (AWG3.LT.AWG4) THEN
            ISPP=ISP
            IDPP=IDP
         ELSE
            ISPP=ISP
            IDPP=IDP1
         END IF
      ELSE IF (AWG1.LT.AWG4) THEN
         ISPP=ISP
         IDPP=IDP
      ELSE
         ISPP=ISP1
         IDPP=IDP1
      END IF
      IF (AWG5.LT.AWG6) THEN
         IF (AWG6.LT.AWG7) THEN
            IF (AWG7.LT.AWG8) THEN
               ISMM=ISM
               IDMM=IDM
            ELSE
               ISMM=ISM
               IDMM=IDM1
            END IF
         ELSE IF (AWG6.LT.AWG8) THEN
            ISMM=ISM
            IDMM=IDM
         ELSE
            ISMM=ISM1
            IDMM=IDM
         END IF
      ELSE IF (AWG5.LT.AWG7) THEN
         IF (AWG7.LT.AWG8) THEN
            ISMM=ISM
            IDMM=IDM
         ELSE
            ISMM=ISM
            IDMM=IDM1
         END IF
      ELSE IF (AWG5.LT.AWG8) THEN
         ISMM=ISM
         IDMM=IDM
      ELSE
         ISMM=ISM1
         IDMM=IDM1
      END IF
!
!     *** fill the arrays *
!
      WWINT(1) = IDP
      WWINT(2) = IDP1
      WWINT(3) = IDM
      WWINT(4) = IDM1
      WWINT(5) = ISP
      WWINT(6) = ISP1
      WWINT(7) = ISM
      WWINT(8) = ISM1
      WWINT(9) = ISLOW
      WWINT(10)= ISHGH
      WWINT(11)= ISCLW
      WWINT(12)= ISCHG
      WWINT(13)= IDLOW
      WWINT(14)= IDHGH
      WWINT(15)= MSC4MI
      WWINT(16)= MSC4MA
      WWINT(17)= MDC4MI
      WWINT(18)= MDC4MA
      WWINT(19)= MSCMAX
      WWINT(20)= MDCMAX
      WWINT(21)= IDPP                                                     40.41
      WWINT(22)= IDMM                                                     40.41
      WWINT(23)= ISPP                                                     40.41
      WWINT(24)= ISMM                                                     40.41
!
      WWAWG(1) = AWG1
      WWAWG(2) = AWG2
      WWAWG(3) = AWG3
      WWAWG(4) = AWG4
      WWAWG(5) = AWG5
      WWAWG(6) = AWG6
      WWAWG(7) = AWG7
      WWAWG(8) = AWG8
!
      WWSWG(1) = SWG1
      WWSWG(2) = SWG2
      WWSWG(3) = SWG3
      WWSWG(4) = SWG4
      WWSWG(5) = SWG5
      WWSWG(6) = SWG6
      WWSWG(7) = SWG7
      WWSWG(8) = SWG8

      ALLOCATE (AF11(MSC4MI:MSC4MA))                                      40.17

!     *** Fill scaling array (f**11)                     ***
!     *** compute the radian frequency**11 for IS=1, MSC ***
!
      DO 100 IS=1, MSC
        AF11(IS) = ( SPCSIG(IS) / ( 2. * PI ) )**11                       30.72
 100  CONTINUE
!
!     *** compute the radian frequency for the IS = MSC+1, ISHGH ***
!
      FREQ   = SPCSIG(MSC) / ( 2. * PI )                                  30.72
      DO 110 IS = MSC+1, ISHGH
        FREQ   = FREQ * XIS
        AF11(IS) = FREQ**11
 110  CONTINUE
!
!     *** compute the radian frequency for IS = 0, ISLOW ***
!
      FREQ   = SPCSIG(1) / ( 2. * PI )                                    30.72
      DO 120 IS = 0, ISLOW, -1
        FREQ   = FREQ / XIS
        AF11(IS) = FREQ**11
 120  CONTINUE
!
!     *** test output ***
!
      IF (ISLOW .LT. MSC4MI .OR. ISHGH .GT. MSC4MA .OR.
     &    IDLOW .LT. MDC4MI .OR. IDHGH .GT. MDC4MA) THEN
        WRITE (PRINTF,900) IXCGRD(1), IYCGRD(1),
     &                     ISLOW, ISHGH, IDLOW, IDHGH,
     &                     MSC4MI,MSC4MA, MDC4MI, MDC4MA
 900    FORMAT ( ' ** Error : array bounds and maxima in subr FAC4WW, ',
     &           ' point ', 2I5,
     &         /,'            ISL,ISH : ',2I4, '   IDL,IDH : ',2I4,
     &         /,'            SMI,SMA : ',2I4, '   DMI,DMA : ',2I4)
      ENDIF
!
      IF (ITEST .GE. 40) THEN
        RADE = 360.0 / ( 2. * PI )
        WRITE(PRINTF,*)
        WRITE(PRINTF,*) ' FAC4WW subroutine '
        WRITE(PRINTF,9000) DELTH4*RADE, DELTH3*RADE, DDIR*RADE, XIS
 9000   FORMAT (' THET3 THET4 DDIR XIS  :',4E12.4)
        WRITE(PRINTF,9011) IDP, IDP1, IDM, IDM1
 9011   FORMAT (' IDP IDP1 IDM IDM1     :',4I5)
        WRITE(PRINTF,9012) WIDP, WIDP1, WIDM, WIDM1
 9012   FORMAT (' WIDP WIDP1 WIDM WIDM1 :',4E12.4)
        WRITE (PRINTF,9013) ISP, ISP1, ISM, ISM1
 9013   FORMAT (' ISP ISP1 ISM ISM1     :',4I5)
        WRITE (PRINTF,9014) WISP, WISP1, WISM, WISM1
 9014   FORMAT (' WISP WISP1 WISM WISM1 :',4E12.4)
        WRITE(PRINTF,9016) ISCLW, ISCHG
 9016   FORMAT (' ICLW ICHG             :',2I5)
        WRITE (PRINTF,9017) AWG1, AWG2, AWG3, AWG4
 9017   FORMAT (' AWG1 AWG2 AWG3 AWG4   :',4E12.4)
        WRITE (PRINTF,9018) AWG5, AWG6, AWG7, AWG8
 9018   FORMAT (' AWG5 AWG6 AWG7 AWG8   :',4E12.4)
        WRITE (PRINTF,9019) MSC4MI, MSC4MA, MDC4MI, MDC4MA
 9019   FORMAT (' S4MI S4MA D4MI D4MA   :',4I6)
        WRITE (PRINTF,9015) ISLOW, ISHGH, IDLOW,IDHGH
 9015   FORMAT (' ISLOW ISHG IDLOW IDHG :',4I5)
        WRITE(PRINTF,*)
      END IF
!
      RETURN
!     End of FAC4WW
      END
!
!******************************************************************
!
      SUBROUTINE RANGE4 (WWINT ,IDDLOW,IDDTOP)                            40.00
!
!******************************************************************
!
      USE SWCOMM3                                                         40.41
      USE SWCOMM4                                                         40.41
      USE OCPCOMM4                                                        40.41
!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering                              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmers: The SWAN team                                |
!   --|-----------------------------------------------------------|--
!
!
!     SWAN (Simulating WAves Nearshore); a third generation wave model
!     Copyright (C) 2008  Delft University of Technology
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!     GNU General Public License for more details.
!
!     A copy of the GNU General Public License is available at
!     http://www.gnu.org/copyleft/gpl.html#SEC3
!     or by writing to the Free Software Foundation, Inc.,
!     59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!
!  0. Authors
!
!     40.00: Nico Booij
!     40.10: IJsbrand Haagsma
!     40.41: Marcel Zijlema
!
!  1. Updates
!
!     40.10, Mar. 00: Made modification for exact quadruplets
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. Purpose :
!
!     calculate the minimum and maximum counters in frequency and
!     directional space which fall with the calculation for the
!     nonlinear wave-wave interactions.
!
!  3. Method :  review for the counters :
!
!                            Frequencies -->
!                 +---+---------------------+---------+- IDHGH
!              d  | 3 :          2          :    2    |
!              i  + - + - - - - - - - - - - + - - - - +- MDC
!              r  |   :                     :         |
!              e  | 3 :  original spectrum  :    1    |
!              c  |   :                     :         |
!              t. + - + - - - - - - - - - - + - - - - +- 1
!                 | 3 :          2          :    2    |
!                 +---+---------------------+---------+- IDLOW
!                 |   |                     |    ^    |
!             ISLOW   1                     MSC  |  ISHGH
!                     ^                          |
!                     |                          |
!                    ISCLW                     ISCHG
!              lowest discrete               highest discrete
!                central bin                   central bin
!
!
!       The directional counters depend on the numerical method that
!       is used.
!
!  4. Parameters :
!
!     INTEGER
!     -------
!     IQUAD         Counter for 4 wave interactions
!     ISLOW,ISHGH   Minimum and maximum counter in frequency space
!     ISCLW,ISCHG   idem for discrete computations
!     IDLOW,IDHGH   Minimum and maximum counters in directional space
!     MSC,MDC       Range of the original arrays
!     ISM1,ISP1,
!     IDM1,IDP1     see subroutine FAC4WW
!     IDDLOW        minimum counter of the bin that is propagated
!                   within a sweep
!     IDDTOP        minimum counter of the bin that is propagated
!                   within a sweep
!
!     array:
!     ------
!     WWINT         counters for the nonlinear interactions
!
!     WWINT ( 1  = IDP      2  = IDP1     3  = IDM     4  = IDM1
!             5  = ISP      6  = ISP1     7  = ISM     8  = ISM1
!             9  = ISLOW    10 = ISHGH    11 = ISCLW   12 = ISCHG
!             13 = IDLOW    14 = IDHGH    15 = MSC4MI  16 = MSC4MA
!             17 = MDC4MI   18 = MDC4MA
!             19 = MSCMAX   20 = MDCMAX )
!
!  5. Subroutines used :
!
!     ---
!
!  6. Called by :
!
!     SOURCE
!
!  7. Common blocks used
!
!
!  9. Source code :
!
!     -----------------------------------------------------------------
!     Calculate :
!       In absence of a current there are always four sectors
!         equal 90 degrees within a sweep that are propagated
!         Extend the boundaries to calculate the source term
!       In presence of a current and if IDTOT .eq. MDC then calculate
!         boundaries for calculation of interaction using the
!         unfolded area.
!     ----------------------------------------------------------
!
!****************************************************************
!
      INTEGER     IDDLOW,IDDTOP                                           40.00
!
      INTEGER     WWINT(*)
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'RANGE4')
!
!     *** Range in directional domain ***
!
      IF ( IQUAD .LT. 3 .AND. IQUAD .GT. 0 ) THEN                         40.10
!       *** counters based on bins which fall within a sweep ***
        WWINT(13) = IDDLOW - MAX( WWINT(4), WWINT(2) )
        WWINT(14) = IDDTOP + MAX( WWINT(4), WWINT(2) )
      ELSE
!       *** counters initially based on full circle ***
        WWINT(13) = 1   - MAX( WWINT(4), WWINT(2) )
        WWINT(14) = MDC + MAX( WWINT(4), WWINT(2) )
      END IF
!
!     *** error message ***
!
      IF (WWINT(9)  .LT. WWINT(15) .OR. WWINT(10) .GT. WWINT(16) .OR.
     &    WWINT(13) .LT. WWINT(17) .OR. WWINT(14) .GT. WWINT(18) ) THEN
        WRITE (PRINTF,900) IXCGRD(1), IYCGRD(1),
     &                     WWINT(9) ,WWINT(10) ,WWINT(13) ,WWINT(14),
     &                     WWINT(15),WWINT(16) ,WWINT(17) ,WWINT(18)
 900    FORMAT ( ' ** Error : array bounds and maxima in subr RANGE4, ',
     &           ' point ', 2I5,
     &         /,'            ISL,ISH : ',2I4, '   IDL,IDH : ',2I4,
     &         /,'            SMI,SMA : ',2I4, '   DMI,DMA : ',2I4)
        IF (ITEST.GE.50) WRITE (PRTEST, 901) MSC, MDC, IDDLOW, IDDTOP
 901    FORMAT (' MSC, MDC, IDDLOW, IDDTOP: ', 4I5)
      ENDIF
!
!     test output
!
      IF (TESTFL .AND. ITEST .GE. 60) THEN
        WRITE(PRTEST,911) WWINT(4), WWINT(2), WWINT(8), WWINT(6)
 911    FORMAT (' RANGE4: IDM1 IDP1 ISM1 ISP1    :',4I5)
        WRITE(PRTEST,916) WWINT(11), WWINT(12), IQUAD
 916    FORMAT (' RANGE4: ISCLW ISCHG IQUAD      :',3I5)
        WRITE (PRTEST,917) WWINT(9), WWINT(10), WWINT(13), WWINT(14)
 917    FORMAT (' RANGE4: ISLOW ISHGH IDLOW IDHGH:',4I5)
        WRITE (PRTEST,919) WWINT(15), WWINT(16), WWINT(17), WWINT(18)
 919    FORMAT (' RANGE4: MS4MI MS4MA MD4MI MD4MA:',4I5)
        WRITE(PRINTF,*)
      END IF
!
      RETURN
!     End of RANGE4
      END
!
!********************************************************************
!
      SUBROUTINE SWSNL1 (WWINT   ,WWAWG   ,WWSWG   ,                      34.00
     &                   IDCMIN  ,IDCMAX  ,UE      ,SA1     ,             40.17
     &                   SA2     ,DA1C    ,DA1P    ,DA1M    ,DA2C    ,
     &                   DA2P    ,DA2M    ,SPCSIG  ,SNLC1   ,KMESPC  ,    30.72
     &                   FACHFR  ,ISSTOP  ,DAL1    ,DAL2    ,DAL3    ,
     &                   SFNL    ,DSNL    ,DEP2    ,AC2     ,IMATDA  ,
     &                   IMATRA  ,PLNL4S  ,PLNL4D  ,                      34.00
     &                   IDDLOW  ,IDDTOP  )                               34.00
!
!********************************************************************
!
      USE SWCOMM3                                                         40.41
      USE SWCOMM4                                                         40.41
      USE OCPCOMM4                                                        40.41
      USE M_SNL4                                                          40.17
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering                              |
!     | Fluid Mechanics Section                                   |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmers: H.L. Tolman, R.C. Ris                        |
!   --|-----------------------------------------------------------|--
!
!
!     SWAN (Simulating WAves Nearshore); a third generation wave model
!     Copyright (C) 2008  Delft University of Technology
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!     GNU General Public License for more details.
!
!     A copy of the GNU General Public License is available at
!     http://www.gnu.org/copyleft/gpl.html#SEC3
!     or by writing to the Free Software Foundation, Inc.,
!     59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!
!  0. Authors
!
!     30.72: IJsbrand Haagsma
!     40.13: Nico Booij
!     40.17: IJsbrand Haagsma
!     40.23: Marcel Zijlema
!     40.41: Marcel Zijlema
!
!  1. Updates
!
!     30.72, Feb. 98: Introduced generic names XCGRID, YCGRID and SPCSIG for SWAN
!     40.17, Dec. 01: Implentation of Multiple DIA
!     40.23, Aug. 02: some corrections
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. Purpose
!
!     Calculate non-linear interaction using the discrete interaction
!     approximation (Hasselmann and Hasselmann 1985; WAMDI group 1988),
!     including the diagonal term for the implicit integration.
!
!     The interactions are calculated for all bin's that fall
!     within a sweep. No additional auxiliary array is required (see
!     SWSNL3)
!
!  3. Method
!
!     Discrete interaction approximation.
!
!     Since the domain in directional domain is by definition not
!     periodic, the spectral space can not beforehand
!     folded to the side angles. This can only be done if the
!     full circle has to be calculated
!
!
!                            Frequencies -->
!                 +---+---------------------+---------+- IDHGH
!              d  | 3 :          2          :    2    |
!              i  + - + - - - - - - - - - - + - - - - +- MDC
!              r  |   :                     :         |
!              e  | 3 :  original spectrum  :    1    |
!              c  |   :                     :         |
!              t. + - + - - - - - - - - - - + - - - - +- 1
!                 | 3 :          2          :    2    |
!                 +---+---------------------+---------+- IDLOW
!                 |   |                     |    ^    |
!             ISLOW   1                     MSC  |    ISHGH
!                     ^                          |
!                     |                          |
!                    ISCLW                     ISCHG
!              lowest discrete               highest discrete
!                central bin                   central bin
!
!                            1 : Extra tail added beyond MSC
!                            2 : Spectrum copied outside ID range
!                            3 : Empty bins at low frequencies
!
!     ISLOW =  1  + ISM1
!     ISHGH = MSC + ISP1 - ISM1
!     ISCLW =  1
!     ISCHG = MSC - ISM1
!     IDLOW =  IDDLOW - MAX(IDM1,IDP1)
!     IDHGH =  IDDTOP + MAX(IDM1,IDP1)
!
!     For the meaning of the counters on the right hand side of the
!     above equations see section 4.
!
!  4. Argument variables
!
!     SPCSIG: Relative frequencies in computational domain in sigma-space 30.72
!
      REAL    SPCSIG(MSC)                                                 30.72
!
!     Data in PARAMETER statements :
!     ----------------------------------------------------------------
!       DAL1    Real  LAMBDA dependend weight factors (see FAC4WW)
!       DAL2    Real
!       DAL3    Real
!       ITHP, ITHP1, ITHM, ITHM1, IFRP, IFRP1, IFRM, IFRM1
!               Int.  Counters of interpolation point relative to
!                     central bin, see figure below (set in FAC4WW).
!       NFRLOW, NFRHGH, NFRCHG, NTHLOW, NTHHGH
!               Int.  Range of calculations, see section 2.
!       AF11    R.A.  Scaling array (Freq**11).
!       AWGn    Real  Interpolation weights, see numbers in fig.
!       SWGn    Real  Id. squared.
!       UE      R.A.  "Unfolded" spectrum.
!       SA1     R.A.  Interaction constribution of first and second
!       SA2     R.A.    quadr. respectively (unfolded space).
!       DA1C, DA1P, DA1M, DA2C, DA2P, DA2M
!               R.A.  Idem for diagonal matrix.
!       PERCIR        full circle or sector
!     ----------------------------------------------------------------
!
!       Relative offsets of interpolation points around central bin
!       "#" and corresponding numbers of AWGn :
!
!               ISM1  ISM
!                5        7    T |
!          IDM1   +------+     H +
!                 |      |     E |      ISP      ISP1
!                 |   \  |     T |       3           1
!           IDM   +------+     A +        +---------+  IDP1
!                6       \8      |        |         |
!                                |        |  /      |
!                           \    +        +---------+  IDP
!                                |      /4           2
!                              \ |  /
!          -+-----+------+-------#--------+---------+----------+
!                                |           FREQ.
!
!  7. Common blocks used
!
!
!  8. Subroutines used
!
!     ---
!
!  9. Subroutines calling
!
!     SOURCE (in SWANCOM1)
!
! 12. Structure
!
!     -------------------------------------------
!       Initialisations.
!       Calculate proportionality constant.
!       Prepare auxiliary spectrum.
!       Calculate interactions :
!       -----------------------------------------
!         Energy at interacting bins
!         Contribution to interactions
!         Fold interactions to side angles
!       -----------------------------------------
!       Put source term together
!     -------------------------------------------
!
! 13. Source text
!
!*************************************************************
!
      INTEGER   IS     ,ID     ,I      ,J      ,                          34.00
     &          ISHGH  ,IDLOW  ,ISP    ,ISP1   ,IDP    ,IDP1   ,
     &          ISM    ,ISM1   ,IDHGH  ,IDM    ,IDM1   ,ISCLW  ,
     &          ISCHG  ,IDDLOW ,IDDTOP                                    34.00
!
      REAL      X      ,X2     ,CONS   ,FACTOR ,SNLCS1 ,SNLCS2 ,SNLCS3,
     &          E00    ,EP1    ,EM1    ,EP2    ,EM2    ,SA1A   ,SA1B  ,
     &          SA2A   ,SA2B   ,KMESPC ,FACHFR ,AWG1   ,AWG2   ,AWG3  ,
     &          AWG4   ,AWG5   ,AWG6   ,AWG7   ,AWG8   ,DAL1   ,DAL2  ,
     &          DAL3   ,SNLC1  ,SWG1   ,SWG2   ,SWG3   ,SWG4   ,SWG5  ,
     &          SWG6   ,SWG7   ,SWG8           ,JACOBI ,SIGPI             34.00
!
      REAL      AC2(MDC,MSC,MCGRD)                    ,
     &          DEP2(MCGRD)                           ,
     &          UE(MSC4MI:MSC4MA , MDC4MI:MDC4MA )    ,
     &          SA1(MSC4MI:MSC4MA , MDC4MI:MDC4MA )   ,
     &          SA2(MSC4MI:MSC4MA , MDC4MI:MDC4MA )   ,
     &          DA1C(MSC4MI:MSC4MA , MDC4MI:MDC4MA )  ,
     &          DA1P(MSC4MI:MSC4MA , MDC4MI:MDC4MA )  ,
     &          DA1M(MSC4MI:MSC4MA , MDC4MI:MDC4MA )  ,
     &          DA2C(MSC4MI:MSC4MA , MDC4MI:MDC4MA )  ,
     &          DA2P(MSC4MI:MSC4MA , MDC4MI:MDC4MA )  ,
     &          DA2M(MSC4MI:MSC4MA , MDC4MI:MDC4MA )  ,
     &          SFNL(MSC4MI:MSC4MA , MDC4MI:MDC4MA )  ,
     &          DSNL(MSC4MI:MSC4MA , MDC4MI:MDC4MA )  ,
     &          IMATDA(MDC,MSC)                       ,
     &          IMATRA(MDC,MSC)                       ,
     &          PLNL4S(MDC,MSC,NPTST)                 ,                   40.00
     &          PLNL4D(MDC,MSC,NPTST)                 ,
     &          WWAWG(*)                              ,
     &          WWSWG(*)
!
      INTEGER   IDCMIN(MSC)        ,
     &          IDCMAX(MSC)        ,
     &          WWINT(*)
!
      LOGICAL   PERCIR
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'SWSNL1')
!
      IDP    = WWINT(1)
      IDP1   = WWINT(2)
      IDM    = WWINT(3)
      IDM1   = WWINT(4)
      ISP    = WWINT(5)
      ISP1   = WWINT(6)
      ISM    = WWINT(7)
      ISM1   = WWINT(8)
      ISLOW  = WWINT(9)
      ISHGH  = WWINT(10)
      ISCLW  = WWINT(11)
      ISCHG  = WWINT(12)
      IDLOW  = WWINT(13)
      IDHGH  = WWINT(14)
!
      AWG1 = WWAWG(1)
      AWG2 = WWAWG(2)
      AWG3 = WWAWG(3)
      AWG4 = WWAWG(4)
      AWG5 = WWAWG(5)
      AWG6 = WWAWG(6)
      AWG7 = WWAWG(7)
      AWG8 = WWAWG(8)
!
      SWG1 = WWSWG(1)
      SWG2 = WWSWG(2)
      SWG3 = WWSWG(3)
      SWG4 = WWSWG(4)
      SWG5 = WWSWG(5)
      SWG6 = WWSWG(6)
      SWG7 = WWSWG(7)
      SWG8 = WWSWG(8)
!
!     *** Initialize auxiliary arrays per gridpoint ***
!
      DO ID = MDC4MI, MDC4MA
        DO IS = MSC4MI, MSC4MA
          UE(IS,ID)   = 0.
          SA1(IS,ID)  = 0.
          SA2(IS,ID)  = 0.
          SFNL(IS,ID) = 0.
          DA1C(IS,ID) = 0.
          DA1P(IS,ID) = 0.
          DA1M(IS,ID) = 0.
          DA2C(IS,ID) = 0.
          DA2P(IS,ID) = 0.
          DA2M(IS,ID) = 0.
          DSNL(IS,ID) = 0.
        ENDDO
      ENDDO
!
!     *** Calculate factor R(X) to calculate the NL wave-wave ***
!     *** interaction for shallow water                       ***
!     *** SNLC1 = 1/GRAV**4                                   ***         40.17
!
      SNLCS1 = PQUAD(3)                                                   34.00
      SNLCS2 = PQUAD(4)                                                   34.00
      SNLCS3 = PQUAD(5)                                                   34.00
      X      = MAX ( 0.75 * DEP2(KCGRD(1)) * KMESPC , 0.5 )
      X2     = MAX ( -1.E15, SNLCS3*X)
      CONS   = SNLC1 * ( 1. + SNLCS1/X * (1.-SNLCS2*X) * EXP(X2))
      JACOBI = 2. * PI
!
!     *** check whether the spectral domain is periodic in ***
!     *** directional space and if so, modify boundaries   ***
!
      PERCIR = .FALSE.
      IF ( IDDLOW .EQ. 1 .AND. IDDTOP .EQ. MDC ) THEN
!       *** periodic in theta -> spectrum can be folded    ***
!       *** (can only be present in presence of a current) ***
        IDCLOW = 1
        IDCHGH = MDC
        IIID   = 0
        PERCIR = .TRUE.
      ELSE
!       *** different sectors per sweep -> extend range with IIID ***
        IIID   = MAX ( IDM1 , IDP1 )
        IDCLOW = IDLOW
        IDCHGH = IDHGH
      ENDIF
!
!     *** Prepare auxiliary spectrum               ***
!     *** set action original spectrum in array UE ***
!
      DO IDDUM = IDLOW - IIID, IDHGH + IIID
        ID = MOD ( IDDUM - 1 + MDC , MDC ) + 1
        DO IS = 1, MSC
          UE(IS,IDDUM) = AC2(ID,IS,KCGRD(1)) * SPCSIG(IS) * JACOBI        30.72
        ENDDO
      ENDDO
!
!     *** set values in area 2 for IS > MSC+1  ***
!
      DO IS = MSC+1, ISHGH
        DO ID = IDLOW - IIID , IDHGH + IIID
          UE (IS,ID) = UE(IS-1,ID) * FACHFR
        ENDDO
      ENDDO
!
!     *** Calculate interactions      ***
!     *** Energy at interacting bins  ***
!
      DO IS = ISCLW, ISCHG
        DO ID = IDCLOW, IDCHGH
          E00    =        UE(IS      ,ID      )
          EP1    = AWG1 * UE(IS+ISP1,ID+IDP1) +
     &             AWG2 * UE(IS+ISP1,ID+IDP ) +
     &             AWG3 * UE(IS+ISP ,ID+IDP1) +
     &             AWG4 * UE(IS+ISP ,ID+IDP )
          EM1    = AWG5 * UE(IS+ISM1,ID-IDM1) +
     &             AWG6 * UE(IS+ISM1,ID-IDM ) +
     &             AWG7 * UE(IS+ISM ,ID-IDM1) +
     &             AWG8 * UE(IS+ISM ,ID-IDM )
!
          EP2    = AWG1 * UE(IS+ISP1,ID-IDP1) +
     &             AWG2 * UE(IS+ISP1,ID-IDP ) +
     &             AWG3 * UE(IS+ISP ,ID-IDP1) +
     &             AWG4 * UE(IS+ISP ,ID-IDP )
          EM2    = AWG5 * UE(IS+ISM1,ID+IDM1) +
     &             AWG6 * UE(IS+ISM1,ID+IDM ) +
     &             AWG7 * UE(IS+ISM ,ID+IDM1) +
     &             AWG8 * UE(IS+ISM ,ID+IDM )
!
!         *** Contribution to interactions                          ***
!         *** CONS is the shallow water factor for the NL interact. ***
!
          FACTOR = CONS * AF11(IS) * E00
!
          SA1A   = E00 * ( EP1*DAL1 + EM1*DAL2 ) * PQUAD(2)               40.17
          SA1B   = SA1A - EP1*EM1*DAL3 * PQUAD(2)                         40.17
          SA2A   = E00 * ( EP2*DAL1 + EM2*DAL2 ) * PQUAD(2)               40.17
          SA2B   = SA2A - EP2*EM2*DAL3 * PQUAD(2)                         40.17
!
          SA1 (IS,ID) = FACTOR * SA1B
          SA2 (IS,ID) = FACTOR * SA2B
!
          IF(ITEST.GE.100 .AND. TESTFL) THEN
            WRITE(PRINTF,9002) E00,EP1,EM1,EP2,EM2
 9002       FORMAT (' E00 EP1 EM1 EP2 EM2  :',5E11.4)
            WRITE(PRINTF,9003) SA1A,SA1B,SA2A,SA2B
 9003       FORMAT (' SA1A SA1B SA2A SA2B  :',4E11.4)
            WRITE(PRINTF,9004) IS,ID,SA1(IS,ID),SA2(IS,ID)
 9004       FORMAT (' IS ID SA1() SA2()    :',2I4,2E12.4)
            WRITE(PRINTF,9005) FACTOR
 9005       FORMAT (' FACTOR               : ',E12.4)
          END IF
!
          DA1C(IS,ID) = CONS * AF11(IS) * ( SA1A + SA1B )
          DA1P(IS,ID) = FACTOR * ( DAL1*E00 - DAL3*EM1 ) * PQUAD(2)       40.23
          DA1M(IS,ID) = FACTOR * ( DAL2*E00 - DAL3*EP1 ) * PQUAD(2)       40.23
!
          DA2C(IS,ID) = CONS * AF11(IS) * ( SA2A + SA2B )
          DA2P(IS,ID) = FACTOR * ( DAL1*E00 - DAL3*EM2 ) * PQUAD(2)       40.23
          DA2M(IS,ID) = FACTOR * ( DAL2*E00 - DAL3*EP2 ) * PQUAD(2)       40.23
        ENDDO
      ENDDO
!
!     *** Fold interactions to side angles if spectral domain ***
!     *** is periodic in directional space                    ***
!
      IF ( PERCIR ) THEN
        DO ID = 1, IDHGH - MDC
          ID0   = 1 - ID
          DO IS = ISCLW, ISCHG
            SA1 (IS,MDC+ID) = SA1 (IS,  ID   )
            SA2 (IS,MDC+ID) = SA2 (IS,  ID   )
            DA1C(IS,MDC+ID) = DA1C(IS,  ID   )
            DA1P(IS,MDC+ID) = DA1P(IS,  ID   )
            DA1M(IS,MDC+ID) = DA1M(IS,  ID   )
            DA2C(IS,MDC+ID) = DA2C(IS,  ID   )
            DA2P(IS,MDC+ID) = DA2P(IS,  ID   )
            DA2M(IS,MDC+ID) = DA2M(IS,  ID   )
!
            SA1 (IS,  ID0 ) = SA1 (IS, MDC+ID0)
            SA2 (IS,  ID0 ) = SA2 (IS, MDC+ID0)
            DA1C(IS,  ID0 ) = DA1C(IS, MDC+ID0)
            DA1P(IS,  ID0 ) = DA1P(IS, MDC+ID0)
            DA1M(IS,  ID0 ) = DA1M(IS, MDC+ID0)
            DA2C(IS,  ID0 ) = DA2C(IS, MDC+ID0)
            DA2P(IS,  ID0 ) = DA2P(IS, MDC+ID0)
            DA2M(IS,  ID0 ) = DA2M(IS, MDC+ID0)
          ENDDO
        ENDDO
      ENDIF
!
!     *** Put source term together (To save space I=IS and J=ID ***
!     *** is used)                                              ***
!
      PI3   = (2. * PI)**3
      DO I = 1, ISSTOP
        SIGPI = SPCSIG(I) * JACOBI                                        30.72
        DO J = IDCMIN(I), IDCMAX(I)
          ID = MOD ( J - 1 + MDC , MDC ) + 1
          SFNL(I,ID) =   - 2. * ( SA1(I,J) + SA2(I,J) )
     &        + AWG1 * ( SA1(I-ISP1,J-IDP1) + SA2(I-ISP1,J+IDP1) )
     &        + AWG2 * ( SA1(I-ISP1,J-IDP ) + SA2(I-ISP1,J+IDP ) )
     &        + AWG3 * ( SA1(I-ISP ,J-IDP1) + SA2(I-ISP ,J+IDP1) )
     &        + AWG4 * ( SA1(I-ISP ,J-IDP ) + SA2(I-ISP ,J+IDP ) )
     &        + AWG5 * ( SA1(I-ISM1,J+IDM1) + SA2(I-ISM1,J-IDM1) )
     &        + AWG6 * ( SA1(I-ISM1,J+IDM ) + SA2(I-ISM1,J-IDM ) )
     &        + AWG7 * ( SA1(I-ISM ,J+IDM1) + SA2(I-ISM ,J-IDM1) )
     &        + AWG8 * ( SA1(I-ISM ,J+IDM ) + SA2(I-ISM ,J-IDM ) )
!
          DSNL(I,ID) =   - 2. * ( DA1C(I,J) + DA2C(I,J) )
     &        + SWG1 * ( DA1P(I-ISP1,J-IDP1) + DA2P(I-ISP1,J+IDP1) )
     &        + SWG2 * ( DA1P(I-ISP1,J-IDP ) + DA2P(I-ISP1,J+IDP ) )
     &        + SWG3 * ( DA1P(I-ISP ,J-IDP1) + DA2P(I-ISP ,J+IDP1) )
     &        + SWG4 * ( DA1P(I-ISP ,J-IDP ) + DA2P(I-ISP ,J+IDP ) )
     &        + SWG5 * ( DA1M(I-ISM1,J+IDM1) + DA2M(I-ISM1,J-IDM1) )
     &        + SWG6 * ( DA1M(I-ISM1,J+IDM ) + DA2M(I-ISM1,J-IDM ) )
     &        + SWG7 * ( DA1M(I-ISM ,J+IDM1) + DA2M(I-ISM ,J-IDM1) )
     &        + SWG8 * ( DA1M(I-ISM ,J+IDM ) + DA2M(I-ISM ,J-IDM ) )
!
!         *** store results in IMATDA and IMATRA ***
!
          IF(TESTFL) THEN
            PLNL4S(ID,I,IPTST) = SFNL(I,ID) / SIGPI                       40.00
            PLNL4D(ID,I,IPTST) = -1. * DSNL(I,ID) / PI3                   40.00
          END IF
!
          IMATRA(ID,I) = IMATRA(ID,I) + SFNL(I,ID) / SIGPI
          IMATDA(ID,I) = IMATDA(ID,I) - DSNL(I,ID) / PI3
!
          IF(ITEST.GE.90 .AND. TESTFL) THEN
            WRITE(PRINTF,9006) I,J,SFNL(I,ID),DSNL(I,ID),
     &       SPCSIG(I)                                                    30.72
 9006       FORMAT (' IS ID SFNL DSNL SPCSIG:',2I4,3E12.4)                30.72
          END IF
!
        ENDDO
      ENDDO
!
!     *** test output ***
!
      IF (ITEST .GE. 50 .AND. TESTFL) THEN
        WRITE(PRINTF,*)
        WRITE(PRINTF,*) ' SWSNL1 subroutine '
        WRITE(PRINTF,9011) IDP, IDP1, IDM, IDM1
 9011   FORMAT (' IDP IDP1 IDM IDM1     :',4I5)
        WRITE (PRINTF,9013) ISP, ISP1, ISM, ISM1
 9013   FORMAT (' ISP ISP1 ISM ISM1     :',4I5)
        WRITE (PRINTF,9015) ISLOW, ISHGH, IDLOW,IDHGH
 9015   FORMAT (' ISLOW ISHGH IDLOW IDHG:',4I5)
        WRITE(PRINTF,9016) ISCLW, ISCHG, IDDLOW, IDDTOP
 9016   FORMAT (' ICLW ICHG IDDLOW IDDTO:',2I5)
        WRITE (PRINTF,9017) AWG1, AWG2, AWG3, AWG4
 9017   FORMAT (' AWG1 AWG2 AWG3 AWG4   :',4E12.4)
        WRITE (PRINTF,9018) AWG5, AWG6, AWG7, AWG8
 9018   FORMAT (' AWG5 AWG6 AWG7 AWG8   :',4E12.4)
        WRITE (PRINTF,9019) MSC4MI, MSC4MA, MDC4MI, MDC4MA
 9019   FORMAT (' S4MI S4MA D4MI D4MA   :',4I6)
        WRITE(PRINTF,9020) SNLC1,X,X2,CONS
 9020   FORMAT (' SNLC1  X  X2  CONS    :',4E12.4)
        WRITE(PRINTF,9021) DEP2(KCGRD(1)),KMESPC, FACHFR, PI
 9021   FORMAT (' DEPTH KMESPC FACHFR PI:',4E12.4)
        WRITE(PRINTF,9023) JACOBI
 9023   FORMAT (' JACOBI                :',E12.4)
        WRITE(PRINTF,*)
      END IF
!
      RETURN
!     End of the subroutine SWSNL1
      END
!
!*******************************************************************
!
      SUBROUTINE SWSNL2 (IDDLOW  ,IDDTOP  ,WWINT   ,                      34.00
     &                   WWAWG   ,UE      ,SA1     ,ISSTOP  ,             40.17
     &                   SA2     ,SPCSIG  ,SNLC1   ,DAL1    ,DAL2    ,    30.72
     &                   DAL3    ,SFNL    ,DEP2    ,AC2     ,KMESPC  ,
     &                                              IMATDA  ,IMATRA  ,    40.23 34.00
     &                   FACHFR  ,PLNL4S           ,IDCMIN  ,IDCMAX  )    34.00
!
!*******************************************************************
!
      USE SWCOMM3                                                         40.41
      USE SWCOMM4                                                         40.41
      USE OCPCOMM4                                                        40.41
      USE M_SNL4                                                          40.17
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering                              |
!     | Fluid Mechanics Section                                   |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmers: H.L. Tolman, R.C. Ris                        |
!   --|-----------------------------------------------------------|--
!
!
!     SWAN (Simulating WAves Nearshore); a third generation wave model
!     Copyright (C) 2008  Delft University of Technology
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!     GNU General Public License for more details.
!
!     A copy of the GNU General Public License is available at
!     http://www.gnu.org/copyleft/gpl.html#SEC3
!     or by writing to the Free Software Foundation, Inc.,
!     59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!
!  0. Authors
!
!     30.72: IJsbrand Haagsma
!     40.13: Nico Booij
!     40.17: IJsbrand Haagsma
!     40.23: Marcel Zijlema
!     40.41: Marcel Zijlema
!
!  1. Updates
!
!     30.72, Feb. 98: Introduced generic names XCGRID, YCGRID and SPCSIG for SWAN
!     40.17, Dec. 01: Implemented Multiple DIA
!     40.23, Aug. 02: rhs and main diagonal adjusted according to Patankar-rules
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. Purpose
!
!     Calculate non-linear interaction using the discrete interaction
!     approximation (Hasselmann and Hasselmann 1985; WAMDI group 1988)
!
!  3. Method
!
!     Discrete interaction approximation.
!
!                            Frequencies -->
!                 +---+---------------------+---------+- IDHGH
!              d  | 3 :          2          :    2    |
!              i  + - + - - - - - - - - - - + - - - - +- MDC
!              r  |   :                     :         |
!              e  | 3 :  original spectrum  :    1    |
!              c  |   :                     :         |
!              t. + - + - - - - - - - - - - + - - - - +- 1
!                 | 3 :          2          :    2    |
!                 +---+---------------------+---------+- IDLOW
!                 |   |                     |     ^   |
!              ISLOW  1                    MSC    |   ISHGH
!                     |                           |
!                   ISCLW                        ISCHG
!              lowest discrete               highest discrete
!                central bin                   central bin
!
!                            1 : Extra tail added beyond MSC
!                            2 : Spectrum copied outside ID range
!                            3 : Empty bins at low frequencies
!
!     ISLOW =  1  + ISM1
!     ISHGH = MSC + ISP1 - ISM1
!     ISCLW =  1
!     ISCHG = MSC - ISM1
!     IDLOW = IDDLOW - MAX(IDM1,IDP1)
!     IDHGH = IDDTOP + MAX(IDM1,IDP1)
!
!       Relative offsets of interpolation points around central bin
!       "#" and corresponding numbers of AWGn :
!
!               ISM1  ISM
!                5        7    T |
!          IDM1   +------+     H +
!                 |      |     E |      ISP      ISP1
!                 |   \  |     T |       3           1
!           IDM   +------+     A +        +---------+  IDP1
!                6       \8      |        |         |
!                                |        |  /      |
!                           \    +        +---------+  IDP
!                                |      /4           2
!                              \ |  /
!          -+-----+------+-------#--------+---------+----------+
!                                |           FREQ.
!
!
!  4. Argument variables
!
!     SPCSIG: Relative frequencies in computational domain in sigma-space 30.72
!
      REAL    SPCSIG(MSC)                                                 30.72
!
!  7. Common blocks used
!
!
!  8. Subroutines used
!
!     ---
!
!  9. Subroutines calling
!
!     SOURCE (in SWANCOM1)
!
! 12. Structure
!
!     -------------------------------------------
!       Initialisations.
!       Calculate proportionality constant.
!       Prepare auxiliary spectrum.
!       Calculate (unfolded) interactions :
!       -----------------------------------------
!         Energy at interacting bins
!         Contribution to interactions
!         Fold interactions to side angles
!       -----------------------------------------
!       Put source term together
!     -------------------------------------------
!
! 13. Source text
!
!*******************************************************************
!
      INTEGER   IS     ,ID     ,I      ,J                      ,ISHGH  ,  34.00
     &          ISSTOP ,ISP    ,ISP1   ,IDP    ,IDP1   ,ISM    ,ISM1   ,
     &          IDM    ,IDM1   ,ISCLW  ,ISCHG  ,                          34.00
     &                  IDLOW  ,IDHGH  ,IDDLOW ,IDDTOP ,IDCLOW ,IDCHGH    34.00
!
      REAL      X      ,X2     ,CONS   ,FACTOR ,SNLCS1 ,SNLCS2 ,SNLCS3 ,
     &          E00    ,EP1    ,EM1    ,EP2    ,EM2    ,SA1A   ,SA1B   ,
     &          SA2A   ,SA2B   ,KMESPC ,FACHFR ,AWG1   ,AWG2   ,AWG3   ,
     &          AWG4   ,AWG5   ,AWG6   ,AWG7   ,AWG8   ,DAL1   ,DAL2   ,
     &          DAL3           ,JACOBI ,SIGPI                             34.00
!
      REAL      AC2(MDC,MSC,MCGRD)                    ,                   30.21
     &          DEP2(MCGRD)                           ,                   30.21
     &          UE(MSC4MI:MSC4MA , MDC4MI:MDC4MA )    ,
     &          SA1(MSC4MI:MSC4MA , MDC4MI:MDC4MA )   ,
     &          SA2(MSC4MI:MSC4MA , MDC4MI:MDC4MA )   ,
     &          SFNL(MSC4MI:MSC4MA , MDC4MI:MDC4MA)   ,
     &          IMATRA(MDC,MSC)                       ,
     &          IMATDA(MDC,MSC)                       ,                   40.23
     &          PLNL4S(MDC,MSC,NPTST)                 ,                   40.00
     &          WWAWG(*)
!
      INTEGER   WWINT(*)         ,
     &          IDCMIN(MSC)      ,
     &          IDCMAX(MSC)
!
      LOGICAL   PERCIR
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'SWSNL2')
!
      IDP    = WWINT(1)
      IDP1   = WWINT(2)
      IDM    = WWINT(3)
      IDM1   = WWINT(4)
      ISP    = WWINT(5)
      ISP1   = WWINT(6)
      ISM    = WWINT(7)
      ISM1   = WWINT(8)
      ISLOW  = WWINT(9)
      ISHGH  = WWINT(10)
      ISCLW  = WWINT(11)
      ISCHG  = WWINT(12)
      IDLOW  = WWINT(13)
      IDHGH  = WWINT(14)
!
      AWG1 = WWAWG(1)
      AWG2 = WWAWG(2)
      AWG3 = WWAWG(3)
      AWG4 = WWAWG(4)
      AWG5 = WWAWG(5)
      AWG6 = WWAWG(6)
      AWG7 = WWAWG(7)
      AWG8 = WWAWG(8)
!
!     *** Initialize auxiliary arrays per gridpoint ***
!
      DO ID = MDC4MI, MDC4MA
        DO IS = MSC4MI, MSC4MA
          UE(IS,ID)   = 0.
          SA1(IS,ID)  = 0.
          SA2(IS,ID)  = 0.
          SFNL(IS,ID) = 0.
        ENDDO
      ENDDO
!
!     *** Calculate prop. constant.                           ***
!     *** Calculate factor R(X) to calculate the NL wave-wave ***
!     *** interaction for shallow water                       ***
!     *** SNLC1 = 1/GRAV**4                                   ***         40.17
!
      SNLCS1 = PQUAD(3)                                                   34.00
      SNLCS2 = PQUAD(4)                                                   34.00
      SNLCS3 = PQUAD(5)                                                   34.00
      X      = MAX ( 0.75 * DEP2(KCGRD(1)) * KMESPC , 0.5 )
      X2     = MAX ( -1.E15, SNLCS3*X)
      CONS   = SNLC1 * ( 1. + SNLCS1/X * (1.-SNLCS2*X) * EXP(X2))
      JACOBI = 2. * PI
!
!     *** check whether the spectral domain is periodic in ***
!     *** direction space and if so modify boundaries      ***
!
      PERCIR = .FALSE.
      IF ( IDDLOW .EQ. 1 .AND. IDDTOP .EQ. MDC ) THEN
!       *** periodic in theta -> spectrum can be folded  ***
!       *** (can only occur in presence of a current)    ***
        IDCLOW = 1
        IDCHGH = MDC
        IIID   = 0
        PERCIR = .TRUE.
      ELSE
!       *** different sectors per sweep -> extend range with IIID ***
        IIID   = MAX ( IDM1 , IDP1 )
        IDCLOW = IDLOW
        IDCHGH = IDHGH
      ENDIF
!
!     *** Prepare auxiliary spectrum               ***
!     *** set action original spectrum in array UE ***
!
      DO IDDUM = IDLOW - IIID , IDHGH + IIID
        ID = MOD ( IDDUM - 1 + MDC , MDC ) + 1
        DO IS = 1, MSC
          UE(IS,IDDUM) = AC2(ID,IS,KCGRD(1)) * SPCSIG(IS) * JACOBI        30.72
        ENDDO
      ENDDO
!
!     *** set values in the areas 2 for IS > MSC+1 ***
!
      DO IS = MSC+1, ISHGH
        DO ID = IDLOW - IIID , IDHGH + IIID
          UE (IS,ID) = UE(IS-1,ID) * FACHFR
        ENDDO
      ENDDO
!
!     *** Calculate interactions      ***
!     *** Energy at interacting bins  ***
!
      DO IS = ISCLW, ISCHG
        DO ID = IDCLOW , IDCHGH
          E00    =        UE(IS      ,ID      )
          EP1    = AWG1 * UE(IS+ISP1,ID+IDP1) +
     &             AWG2 * UE(IS+ISP1,ID+IDP ) +
     &             AWG3 * UE(IS+ISP ,ID+IDP1) +
     &             AWG4 * UE(IS+ISP ,ID+IDP )
          EM1    = AWG5 * UE(IS+ISM1,ID-IDM1) +
     &             AWG6 * UE(IS+ISM1,ID-IDM ) +
     &             AWG7 * UE(IS+ISM ,ID-IDM1) +
     &             AWG8 * UE(IS+ISM ,ID-IDM )
!
          EP2    = AWG1 * UE(IS+ISP1,ID-IDP1) +
     &             AWG2 * UE(IS+ISP1,ID-IDP ) +
     &             AWG3 * UE(IS+ISP ,ID-IDP1) +
     &             AWG4 * UE(IS+ISP ,ID-IDP )
          EM2    = AWG5 * UE(IS+ISM1,ID+IDM1) +
     &             AWG6 * UE(IS+ISM1,ID+IDM ) +
     &             AWG7 * UE(IS+ISM ,ID+IDM1) +
     &             AWG8 * UE(IS+ISM ,ID+IDM )
!
!         *** Contribution to interactions                          ***
!         *** CONS is the shallow water factor for the NL interact. ***
!
          FACTOR = CONS * AF11(IS) * E00
!
          SA1A   = E00 * ( EP1*DAL1 + EM1*DAL2 ) * PQUAD(2)               40.17
          SA1B   = SA1A - EP1*EM1*DAL3 * PQUAD(2)                         40.17
          SA2A   = E00 * ( EP2*DAL1 + EM2*DAL2 ) * PQUAD(2)               40.17
          SA2B   = SA2A - EP2*EM2*DAL3 * PQUAD(2)                         40.17
!
          SA1 (IS,ID) = FACTOR * SA1B
          SA2 (IS,ID) = FACTOR * SA2B
!
          IF(ITEST.GE.100 .AND. TESTFL) THEN
            WRITE(PRINTF,9002) E00,EP1,EM1,EP2,EM2
 9002       FORMAT (' E00 EP1 EM1 EP2 EM2  :',5E11.4)
            WRITE(PRINTF,9003) SA1A,SA1B,SA2A,SA2B
 9003       FORMAT (' SA1A SA1B SA2A SA2B  :',4E11.4)
            WRITE(PRINTF,9004) IS,ID,SA1(IS,ID),SA2(IS,ID)
 9004       FORMAT (' IS ID SA1() SA2()    :',2I4,2E12.4)
            WRITE(PRINTF,9005) FACTOR ,ISLOW
 9005       FORMAT (' FACTOR ISLOW         : ',E12.4,I4)
          END IF
!
        ENDDO
      ENDDO
!
!     *** Fold interactions to side angles if spectral domain ***
!     *** is periodic in directional space                    ***
!
      IF ( PERCIR ) THEN
        DO ID = 1, IDHGH - MDC
          ID0   = 1 - ID
          DO IS = ISCLW, ISCHG
            SA1 (IS,MDC+ID) = SA1 (IS ,  ID    )
            SA2 (IS,MDC+ID) = SA2 (IS ,  ID    )
            SA1 (IS,  ID0 ) = SA1 (IS , MDC+ID0)
            SA2 (IS,  ID0 ) = SA2 (IS , MDC+ID0)
          ENDDO
        ENDDO
      ENDIF
!
!     ***  Put source term together (To save space I=IS and J=ID ***
!     ***  is used)                                              ***
!
      DO I = 1, ISSTOP
        SIGPI = SPCSIG(I) * JACOBI                                        30.72
        DO J = IDCMIN(I), IDCMAX(I)
          ID = MOD ( J - 1 + MDC , MDC ) + 1
          SFNL(I,ID) =   - 2. * ( SA1(I,J) + SA2(I,J) )
     &        + AWG1 * ( SA1(I-ISP1,J-IDP1) + SA2(I-ISP1,J+IDP1) )
     &        + AWG2 * ( SA1(I-ISP1,J-IDP ) + SA2(I-ISP1,J+IDP ) )
     &        + AWG3 * ( SA1(I-ISP ,J-IDP1) + SA2(I-ISP ,J+IDP1) )
     &        + AWG4 * ( SA1(I-ISP ,J-IDP ) + SA2(I-ISP ,J+IDP ) )
     &        + AWG5 * ( SA1(I-ISM1,J+IDM1) + SA2(I-ISM1,J-IDM1) )
     &        + AWG6 * ( SA1(I-ISM1,J+IDM ) + SA2(I-ISM1,J-IDM ) )
     &        + AWG7 * ( SA1(I-ISM ,J+IDM1) + SA2(I-ISM ,J-IDM1) )
     &        + AWG8 * ( SA1(I-ISM ,J+IDM ) + SA2(I-ISM ,J-IDM ) )
!
!         *** store results in rhv ***
!         *** store results in rhs and main diagonal according ***        40.23
!         *** to Patankar-rules                                ***        40.23
!
          IF(TESTFL) PLNL4S(ID,I,IPTST) =  SFNL(I,ID) / SIGPI             40.00
          IF (SFNL(I,ID).GT.0.) THEN                                      40.23
             IMATRA(ID,I) = IMATRA(ID,I) + SFNL(I,ID) / SIGPI
          ELSE
             IMATDA(ID,I) = IMATDA(ID,I) - SFNL(I,ID) /
     &                      MAX(1.E-18,AC2(ID,I,KCGRD(1))*SIGPI)
          END IF
!
        ENDDO
      ENDDO
!
!     *** test output ***
!
      IF (ITEST .GE. 40 .AND. TESTFL) THEN
        WRITE(PRINTF,*) ' SWSNL2 subroutine '
        WRITE(PRINTF,9011) IDP, IDP1, IDM, IDM1
 9011   FORMAT (' IDP IDP1 IDM IDM1     :',4I5)
        WRITE (PRINTF,9013) ISP, ISP1, ISM, ISM1
 9013   FORMAT (' ISP ISP1 ISM ISM1     :',4I5)
        WRITE (PRINTF,9015) ISHGH, IDDLOW, IDDTOP
 9015   FORMAT (' ISHG IDDLOW IDDTOP    :',3I5)
        WRITE(PRINTF,9016) ISCLW, ISCHG, IDLOW, IDHGH
 9016   FORMAT (' ICLW ICHG IDLOW IDHGH :',4I5)
        WRITE (PRINTF,9017) AWG1, AWG2, AWG3, AWG4
 9017   FORMAT (' AWG1 AWG2 AWG3 AWG4   :',4E12.4)
        WRITE (PRINTF,9018) AWG5, AWG6, AWG7, AWG8
 9018   FORMAT (' AWG5 AWG6 AWG7 AWG8   :',4E12.4)
        WRITE (PRINTF,9019) MSC4MI, MSC4MA, MDC4MI, MDC4MA
 9019   FORMAT (' S4MI S4MA D4MI D4MA   :',4I6)
        WRITE(PRINTF,9020) SNLC1,X,X2,CONS
 9020   FORMAT (' SNLC1  X  X2  CONS    :',4E12.4)
        WRITE(PRINTF,9021) DEP2(KCGRD(1)),KMESPC, FACHFR,PI
 9021   FORMAT (' DEPTH KMESPC FACHFR PI:',4E12.4)
        WRITE(PRINTF,9023) JACOBI,ISLOW
 9023   FORMAT (' JACOBI  ISLOW         :',E12.4,I4)
        WRITE(PRINTF,*)
      END IF
!
      RETURN
!     End of SWSNL2
      END
!
!************************************************************
!

      SUBROUTINE SWSNL3 (                  WWINT   ,WWAWG   ,             40.17
     &                   UE      ,SA1     ,SA2     ,SPCSIG  ,SNLC1   ,    40.17
     &                   DAL1    ,DAL2    ,DAL3    ,SFNL    ,DEP2    ,    40.17
     &                   AC2     ,KMESPC  ,MEMNL4  ,FACHFR           )    40.17
!
!*******************************************************************
!
      USE SWCOMM3                                                         40.41
      USE SWCOMM4                                                         40.41
      USE OCPCOMM4                                                        40.41
      USE M_SNL4                                                          40.17
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering                              |
!     | Fluid Mechanics Section                                   |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmers: H.L. Tolman, R.C. Ris                        |
!   --|-----------------------------------------------------------|--
!
!
!     SWAN (Simulating WAves Nearshore); a third generation wave model
!     Copyright (C) 2008  Delft University of Technology
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!     GNU General Public License for more details.
!
!     A copy of the GNU General Public License is available at
!     http://www.gnu.org/copyleft/gpl.html#SEC3
!     or by writing to the Free Software Foundation, Inc.,
!     59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!
!  0. Authors
!
!     30.72: IJsbrand Haagsma
!     40.17: IJsbrand Haagsma
!     40.41: Marcel Zijlema
!
!  1. Updates
!
!     30.72, Feb. 98: Introduced generic names XCGRID, YCGRID and SPCSIG for SWAN
!     40.17, Dec. 01: Implemented Multiple DIA
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. Purpose
!
!     Calculate non-linear interaction using the discrete interaction
!     approximation (Hasselmann and Hasselmann 1985; WAMDI group 1988)
!     for the full circle (option if a current is present). Note: using
!     this subroutine requires an additional array with size
!     (MXC*MYC*MDC*MSC). This requires more internal memory but can
!     speed up the computations sigificantly if a current is present.
!
!  3. Method
!
!     Discrete interaction approximation. To make interpolation simple,
!     the interactions are calculated in a "folded" space.
!
!                            Frequencies -->
!                 +---+---------------------+---------+- IDHGH
!              d  | 3 :          2          :    2    |
!              i  + - + - - - - - - - - - - + - - - - +- MDC
!              r  |   :                     :         |
!              e  | 3 :  original spectrum  :    1    |
!              c  |   :                     :         |
!              t. + - + - - - - - - - - - - + - - - - +- 1
!                 | 3 :          2          :    2    |
!                 +---+---------------------+---------+- IDLOW
!                 |   |                     |     ^   |
!              ISLOW  1                    MSC    |   ISHGH
!                     |                           |
!                   ISCLW                        ISCHG
!              lowest discrete               highest discrete
!                central bin                   central bin
!
!                            1 : Extra tail added beyond MSC
!                            2 : Spectrum copied outside ID range
!                            3 : Empty bins at low frequencies
!
!     ISLOW =  1  + ISM1
!     ISHGH = MSC + ISP1 - ISM1
!     ISCLW =  1
!     ISCHG = MSC - ISM1
!     IDLOW =  1  - MAX(IDM1,IDP1)
!     IDHGH = MDC + MAX(IDM1,IDP1)
!
!       Relative offsets of interpolation points around central bin
!       "#" and corresponding numbers of AWGn :
!
!               ISM1  ISM
!                5        7    T |
!          IDM1   +------+     H +
!                 |      |     E |      ISP      ISP1
!                 |   \  |     T |       3           1
!           IDM   +------+     A +        +---------+  IDP1
!                6       \8      |        |         |
!                                |        |  /      |
!                           \    +        +---------+  IDP
!                                |      /4           2
!                              \ |  /
!          -+-----+------+-------#--------+---------+----------+
!                                |           FREQ.
!
!
!  4. Argument variables
!
!     MCGRD : number of wet grid points of the computational grid
!     MDC   : grid points in theta-direction of computational grid
!     MDC4MA: highest array counter in directional space (Snl4)
!     MDC4MI: lowest array counter in directional space (Snl4)
!     MSC   : grid points in sigma-direction of computational grid
!     MSC4MA: highest array counter in frequency space (Snl4)
!     MSC4MI: lowest array counter in frequency space (Snl4)
!     WWINT : counters for quadruplet interactions
!
      INTEGER WWINT(*)
!
!     AC2   : action density
!     AF11  : scaling frequency
!     DAL1  : coefficient for the quadruplet interactions
!     DAL2  : coefficient for the quadruplet interactions
!     DAL3  : coefficient for the quadruplet interactions
!     DEP2  : depth
!     FACHFR
!     KMESPC: mean average wavenumber over full spectrum
!     MEMNL4
!     PI    : circular constant
!     SA1   : interaction contribution of first quadruplet (unfolded space)
!     SA2   : interaction contribution of second quadruplet (unfolded space)
!     SFNL
!     SNLC1
!     SPCSIG: relative frequencies in computational domain in sigma-space
!     UE    : "unfolded" spectrum
!     WWAWG : weight coefficients for the quadruplet interactions
!
      REAL    DAL1, DAL2, DAL3, FACHFR, KMESPC, SNLC1                     40.17
      REAL    AC2(MDC,MSC,MCGRD)
      REAL    DEP2(MCGRD)
      REAL    MEMNL4(MDC,MSC,MCGRD)
      REAL    SA1(MSC4MI:MSC4MA,MDC4MI:MDC4MA)
      REAL    SA2(MSC4MI:MSC4MA,MDC4MI:MDC4MA)
      REAL    SFNL(MSC4MI:MSC4MA,MDC4MI:MDC4MA)
      REAL    SPCSIG(MSC)                                                 30.72
      REAL    UE(MSC4MI:MSC4MA,MDC4MI:MDC4MA)
      REAL    WWAWG(*)
!
!  7. Common blocks used
!
!
!  8. Subroutines used
!
!     ---
!
!  9. Subroutines calling
!
!     SOURCE (in SWANCOM1)
!
! 12. Structure
!
!     -------------------------------------------
!       Initialisations.
!       Calculate proportionality constant.
!       Prepare auxiliary spectrum.
!       Calculate (unfolded) interactions :
!       -----------------------------------------
!         Energy at interacting bins
!         Contribution to interactions
!         Fold interactions to side angles
!       -----------------------------------------
!       Put source term together
!     -------------------------------------------
!
! 13. Source text
!
!*******************************************************************
!
      INTEGER   IS      ,ID      ,ID0     ,I       ,J       ,
     &          ISHGH   ,IDLOW   ,IDHGH   ,ISP     ,ISP1    ,
     &          IDP     ,IDP1    ,ISM     ,ISM1    ,IDM     ,IDM1    ,
     &          ISCLW   ,ISCHG
!
      REAL      X       ,X2      ,CONS    ,FACTOR  ,SNLCS2  ,
     &          SNLCS3  ,E00     ,EP1     ,EM1     ,EP2     ,EM2     ,
     &          SA1A    ,SA1B    ,SA2A    ,SA2B    ,
     &          AWG1    ,AWG2    ,AWG3    ,AWG4    ,AWG5    ,AWG6    ,
     &          AWG7    ,AWG8    ,
     &          JACOBI  ,SIGPI
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'SWSNL3')
!
      IDP    = WWINT(1)
      IDP1   = WWINT(2)
      IDM    = WWINT(3)
      IDM1   = WWINT(4)
      ISP    = WWINT(5)
      ISP1   = WWINT(6)
      ISM    = WWINT(7)
      ISM1   = WWINT(8)
      ISLOW  = WWINT(9)
      ISHGH  = WWINT(10)
      ISCLW  = WWINT(11)
      ISCHG  = WWINT(12)
      IDLOW  = WWINT(13)
      IDHGH  = WWINT(14)
!
      AWG1 = WWAWG(1)
      AWG2 = WWAWG(2)
      AWG3 = WWAWG(3)
      AWG4 = WWAWG(4)
      AWG5 = WWAWG(5)
      AWG6 = WWAWG(6)
      AWG7 = WWAWG(7)
      AWG8 = WWAWG(8)
!
!     *** Initialize auxiliary arrays per gridpoint ***
!
      DO ID = MDC4MI, MDC4MA
        DO IS = MSC4MI, MSC4MA
          UE(IS,ID)   = 0.
          SA1(IS,ID)  = 0.
          SA2(IS,ID)  = 0.
          SFNL(IS,ID) = 0.
        ENDDO
      ENDDO
!
!     *** Calculate prop. constant.                           ***
!     *** Calculate factor R(X) to calculate the NL wave-wave ***
!     *** interaction for shallow water                       ***
!     *** SNLC1 = 1/GRAV**4                                   ***         40.17
!
      SNLCS1 = PQUAD(3)                                                   34.00
      SNLCS2 = PQUAD(4)                                                   34.00
      SNLCS3 = PQUAD(5)                                                   34.00
      X      = MAX ( 0.75 * DEP2(KCGRD(1)) * KMESPC , 0.5 )               30.21
      X2     = MAX ( -1.E15, SNLCS3*X)
      CONS   = SNLC1 * ( 1. + SNLCS1/X * (1.-SNLCS2*X) * EXP(X2))
      JACOBI = 2. * PI
!
!     *** extend the area with action density at periodic boundaries ***
!
      DO IDDUM = IDLOW, IDHGH
        ID = MOD ( IDDUM - 1 + MDC , MDC ) + 1
        DO IS=1, MSC
          UE (IS,IDDUM) = AC2(ID,IS,KCGRD(1)) * SPCSIG(IS) * JACOBI       30.72
        ENDDO
      ENDDO
!
      DO IS = MSC+1, ISHGH
        DO ID = IDLOW, IDHGH
          UE(IS,ID) = UE(IS-1,ID) * FACHFR
        ENDDO
      ENDDO
!
!     *** Calculate (unfolded) interactions ***
!     *** Energy at interacting bins        ***
!
      DO IS = ISCLW, ISCHG
        DO ID = 1, MDC
          E00    =        UE(IS      ,ID      )
          EP1    = AWG1 * UE(IS+ISP1,ID+IDP1) +
     &             AWG2 * UE(IS+ISP1,ID+IDP ) +
     &             AWG3 * UE(IS+ISP ,ID+IDP1) +
     &             AWG4 * UE(IS+ISP ,ID+IDP )
          EM1    = AWG5 * UE(IS+ISM1,ID-IDM1) +
     &             AWG6 * UE(IS+ISM1,ID-IDM ) +
     &             AWG7 * UE(IS+ISM ,ID-IDM1) +
     &             AWG8 * UE(IS+ISM ,ID-IDM )
          EP2    = AWG1 * UE(IS+ISP1,ID-IDP1) +
     &             AWG2 * UE(IS+ISP1,ID-IDP ) +
     &             AWG3 * UE(IS+ISP ,ID-IDP1) +
     &             AWG4 * UE(IS+ISP ,ID-IDP )
          EM2    = AWG5 * UE(IS+ISM1,ID+IDM1) +
     &             AWG6 * UE(IS+ISM1,ID+IDM ) +
     &             AWG7 * UE(IS+ISM ,ID+IDM1) +
     &             AWG8 * UE(IS+ISM ,ID+IDM )
!
!         Contribution to interactions
!
          FACTOR = CONS * AF11(IS) * E00
!
          SA1A   = E00 * ( EP1*DAL1 + EM1*DAL2 ) * PQUAD(2)               40.17
          SA1B   = SA1A - EP1*EM1*DAL3 * PQUAD(2)                         40.17
          SA2A   = E00 * ( EP2*DAL1 + EM2*DAL2 ) * PQUAD(2)               40.17
          SA2B   = SA2A - EP2*EM2*DAL3 * PQUAD(2)                         40.17
!
          SA1 (IS,ID) = FACTOR * SA1B
          SA2 (IS,ID) = FACTOR * SA2B
!
          IF(ITEST.GE.100 .AND. TESTFL) THEN
            WRITE(PRINTF,9002) E00,EP1,EM1,EP2,EM2
 9002       FORMAT (' E00 EP1 EM1 EP2 EM2  :',5E11.4)
            WRITE(PRINTF,9003) SA1A,SA1B,SA2A,SA2B
 9003       FORMAT (' SA1A SA1B SA2A SA2B  :',4E11.4)
            WRITE(PRINTF,9004) IS,ID,SA1(IS,ID),SA2(IS,ID)
 9004       FORMAT (' IS ID SA1() SA2()    :',2I4,2E12.4)
            WRITE(PRINTF,9005) FACTOR,JACOBI
 9005       FORMAT (' FACTOR JACOBI        : ',2E12.4)
          END IF
!
        ENDDO
      ENDDO
!
!     *** Fold interactions to side angles -> domain in theta is ***
!     *** periodic                                               ***
!
      DO ID = 1, IDHGH - MDC
        ID0   = 1 - ID
        DO IS = ISCLW, ISCHG
          SA1 (IS,MDC+ID) = SA1 (IS,  ID   )
          SA2 (IS,MDC+ID) = SA2 (IS,  ID   )
          SA1 (IS,  ID0 ) = SA1 (IS,MDC+ID0)
          SA2 (IS,  ID0 ) = SA2 (IS,MDC+ID0)
        ENDDO
      ENDDO
!
!     *** Put source term together (To save space I=IS and ***
!     *** J=MDC is used)  ----                             ***
!
      DO I = 1, MSC
        SIGPI = SPCSIG(I) * JACOBI                                        30.72
        DO J = 1, MDC
          SFNL(I,J) =   - 2. * ( SA1(I,J) + SA2(I,J) )
     &        + AWG1 * ( SA1(I-ISP1,J-IDP1) + SA2(I-ISP1,J+IDP1) )
     &        + AWG2 * ( SA1(I-ISP1,J-IDP ) + SA2(I-ISP1,J+IDP ) )
     &        + AWG3 * ( SA1(I-ISP ,J-IDP1) + SA2(I-ISP ,J+IDP1) )
     &        + AWG4 * ( SA1(I-ISP ,J-IDP ) + SA2(I-ISP ,J+IDP ) )
     &        + AWG5 * ( SA1(I-ISM1,J+IDM1) + SA2(I-ISM1,J-IDM1) )
     &        + AWG6 * ( SA1(I-ISM1,J+IDM ) + SA2(I-ISM1,J-IDM ) )
     &        + AWG7 * ( SA1(I-ISM ,J+IDM1) + SA2(I-ISM ,J-IDM1) )
     &        + AWG8 * ( SA1(I-ISM ,J+IDM ) + SA2(I-ISM ,J-IDM ) )
!
!         *** store value in auxiliary array and use values in ***
!         *** next four sweeps (see subroutine FILNL3)         ***
!
          MEMNL4(J,I,KCGRD(1)) = SFNL(I,J) / SIGPI                        30.21
        ENDDO
      ENDDO
!
!     *** test output ***
!
      IF (ITEST .GE. 50 .AND. TESTFL) THEN
        WRITE(PRINTF,*)
        WRITE(PRINTF,*) ' SWSNL3 subroutine '
        WRITE(PRINTF,9011) IDP, IDP1, IDM, IDM1
 9011   FORMAT (' IDP IDP1 IDM IDM1     :',4I5)
        WRITE (PRINTF,9013) ISP, ISP1, ISM, ISM1
 9013   FORMAT (' ISP ISP1 ISM ISM1     :',4I5)
        WRITE (PRINTF,9015) ISLOW, ISHGH, IDLOW, IDHGH
 9015   FORMAT (' ISLOW ISHG IDLOW IDHG :',4I5)
        WRITE(PRINTF,9016) ISCLW, ISCHG, JACOBI
 9016   FORMAT (' ICLW ICHG JACOBI      :',2I5,E12.4)
        WRITE (PRINTF,9017) AWG1, AWG2, AWG3, AWG4
 9017   FORMAT (' AWG1 AWG2 AWG3 AWG4   :',4E12.4)
        WRITE (PRINTF,9018) AWG5, AWG6, AWG7, AWG8
 9018   FORMAT (' AWG5 AWG6 AWG7 AWG8   :',4E12.4)
        WRITE (PRINTF,9019) MSC4MI, MSC4MA, MDC4MI, MDC4MA
 9019   FORMAT (' S4MI S4MA D4MI D4MA   :',4I6)
        WRITE(PRINTF,9020) SNLC1,X,X2,CONS
 9020   FORMAT (' SNLC1  X  X2  CONS    :',4E12.4)
        WRITE(PRINTF,9021) DEP2(KCGRD(1)),KMESPC,FACHFR,PI
 9021   FORMAT (' DEPTH KMESPC FACHFR PI:',4E12.4)
        WRITE(PRINTF,*)
!
!       *** value source term in every bin ***
!
        IF(ITEST.GE. 150 ) THEN
          DO I=1, MSC
            DO J=1, MDC
              WRITE(PRINTF,2006) I,J,MEMNL4(J,I,KCGRD(1)),SFNL(I,J),      30.21
     &                           SPCSIG(I)                                30.72
 2006         FORMAT (' I J MEMNL() SFNL() SPCSIG:',2I4,3E12.4)           30.72
            ENDDO
          ENDDO
        END IF
      END IF
!
      RETURN
!
      END SUBROUTINE SWSNL3
!
!*******************************************************************
!
      SUBROUTINE SWSNL4 (WWINT   ,WWAWG   ,
     &                   SPCSIG  ,SNLC1   ,
     &                   DAL1    ,DAL2    ,DAL3    ,DEP2    ,
     &                   AC2     ,KMESPC  ,MEMNL4  ,FACHFR  ,
     &                   IDIA    ,ITER    )
!
!*******************************************************************
!
      USE SWCOMM3                                                         40.41
      USE SWCOMM4                                                         40.41
      USE OCPCOMM4                                                        40.41
      USE M_SNL4
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering                              |
!     | Fluid Mechanics Section                                   |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmers: H.L. Tolman, R.C. Ris                        |
!   --|-----------------------------------------------------------|--
!
!
!     SWAN (Simulating WAves Nearshore); a third generation wave model
!     Copyright (C) 2008  Delft University of Technology
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!     GNU General Public License for more details.
!
!     A copy of the GNU General Public License is available at
!     http://www.gnu.org/copyleft/gpl.html#SEC3
!     or by writing to the Free Software Foundation, Inc.,
!     59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!
!  0. Authors
!
!     40.17: IJsbrand Haagsma
!     40.41: Marcel Zijlema
!
!  1. Updates
!
!     40.17, Dec. 01: New Subroutine based on SWSNL3
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. Purpose
!
!     Calculate non-linear interaction using the discrete interaction
!     approximation (Hasselmann and Hasselmann 1985; WAMDI group 1988)
!     for the full circle (option if a current is present). Note: using
!     this subroutine requires an additional array with size
!     (MXC*MYC*MDC*MSC). This requires more internal memory but can
!     speed up the computations sigificantly if a current is present.
!
!  3. Method
!
!     Discrete interaction approximation. To make interpolation simple,
!     the interactions are calculated in a "folded" space.
!
!                            Frequencies -->
!                 +---+---------------------+---------+- IDHGH
!              d  | 3 :          2          :    2    |
!              i  + - + - - - - - - - - - - + - - - - +- MDC
!              r  |   :                     :         |
!              e  | 3 :  original spectrum  :    1    |
!              c  |   :                     :         |
!              t. + - + - - - - - - - - - - + - - - - +- 1
!                 | 3 :          2          :    2    |
!                 +---+---------------------+---------+- IDLOW
!                 |   |                     |     ^   |
!              ISLOW  1                    MSC    |   ISHGH
!                     |                           |
!                   ISCLW                        ISCHG
!              lowest discrete               highest discrete
!                central bin                   central bin
!
!                            1 : Extra tail added beyond MSC
!                            2 : Spectrum copied outside ID range
!                            3 : Empty bins at low frequencies
!
!     ISLOW =  1  + ISM1
!     ISHGH = MSC + ISP1 - ISM1
!     ISCLW =  1
!     ISCHG = MSC - ISM1
!     IDLOW =  1  - MAX(IDM1,IDP1)
!     IDHGH = MDC + MAX(IDM1,IDP1)
!
!       Relative offsets of interpolation points around central bin
!       "#" and corresponding numbers of AWGn :
!
!               ISM1  ISM
!                5        7    T |
!          IDM1   +------+     H +
!                 |      |     E |      ISP      ISP1
!                 |   \  |     T |       3           1
!           IDM   +------+     A +        +---------+  IDP1
!                6       \8      |        |         |
!                                |        |  /      |
!                           \    +        +---------+  IDP
!                                |      /4           2
!                              \ |  /
!          -+-----+------+-------#--------+---------+----------+
!                                |           FREQ.
!
!
!  4. Argument variables
!
!     MCGRD : number of wet grid points of the computational grid
!     MDC   : grid points in theta-direction of computational grid
!     MDC4MA: highest array counter in directional space (Snl4)
!     MDC4MI: lowest array counter in directional space (Snl4)
!     MSC   : grid points in sigma-direction of computational grid
!     MSC4MA: highest array counter in frequency space (Snl4)
!     MSC4MI: lowest array counter in frequency space (Snl4)
!     WWINT : counters for quadruplet interactions
!
      INTEGER WWINT(*)
      INTEGER IDIA
!
!     AC2   : action density
!     AF11  : scaling frequency
!     DAL1  : coefficient for the quadruplet interactions
!     DAL2  : coefficient for the quadruplet interactions
!     DAL3  : coefficient for the quadruplet interactions
!     DEP2  : depth
!     FACHFR
!     KMESPC: mean average wavenumber over full spectrum
!     MEMNL4
!     PI    : circular constant
!     SA1   : interaction contribution of first quadruplet (unfolded space)
!     SA2   : interaction contribution of second quadruplet (unfolded space)
!     SFNL
!     SNLC1
!     SPCSIG: relative frequencies in computational domain in sigma-space
!     UE    : "unfolded" spectrum
!     WWAWG : weight coefficients for the quadruplet interactions
!
      REAL    DAL1, DAL2, DAL3, FACHFR, KMESPC, SNLC1
      REAL    AC2(MDC,MSC,MCGRD)
      REAL    DEP2(MCGRD)
      REAL    MEMNL4(MDC,MSC,MCGRD)
      REAL    SPCSIG(MSC)
      REAL    WWAWG(*)
!
!  6. Local variables
!
      REAL, ALLOCATABLE :: SA1(:,:), SA2(:,:), SFNL(:,:), UE(:,:)
!
!  7. Common blocks used
!
!
!  8. Subroutines used
!
!     ---
!
!  9. Subroutines calling
!
!     SOURCE (in SWANCOM1)
!
! 12. Structure
!
!     -------------------------------------------
!       Initialisations.
!       Calculate proportionality constant.
!       Prepare auxiliary spectrum.
!       Calculate (unfolded) interactions :
!       -----------------------------------------
!         Energy at interacting bins
!         Contribution to interactions
!         Fold interactions to side angles
!       -----------------------------------------
!       Put source term together
!     -------------------------------------------
!
! 13. Source text
!
!*******************************************************************
!
      INTEGER   IS      ,ID      ,ID0     ,I       ,J       ,
     &          ISHGH   ,IDLOW   ,IDHGH   ,ISP     ,ISP1    ,
     &          IDP     ,IDP1    ,ISM     ,ISM1    ,IDM     ,IDM1    ,
     &          ISCLW   ,ISCHG
!
      REAL      X       ,X2      ,CONS    ,FACTOR  ,SNLCS2  ,
     &          SNLCS3  ,E00     ,EP1     ,EM1     ,EP2     ,EM2     ,
     &          SA1A    ,SA1B    ,SA2A    ,SA2B    ,
     &          AWG1    ,AWG2    ,AWG3    ,AWG4    ,AWG5    ,AWG6    ,
     &          AWG7    ,AWG8    ,
     &          JACOBI  ,SIGPI
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'SWSNL4')

      ALLOCATE(SA1(MSC4MI:MSC4MA,MDC4MI:MDC4MA))
      ALLOCATE(SA2(MSC4MI:MSC4MA,MDC4MI:MDC4MA))
      ALLOCATE(SFNL(MSC4MI:MSC4MA,MDC4MI:MDC4MA))
      ALLOCATE(UE(MSC4MI:MSC4MA,MDC4MI:MDC4MA))
!
      IDP    = WWINT(1)
      IDP1   = WWINT(2)
      IDM    = WWINT(3)
      IDM1   = WWINT(4)
      ISP    = WWINT(5)
      ISP1   = WWINT(6)
      ISM    = WWINT(7)
      ISM1   = WWINT(8)
      ISLOW  = WWINT(9)
      ISHGH  = WWINT(10)
      ISCLW  = WWINT(11)
      ISCHG  = WWINT(12)
      IDLOW  = WWINT(13)
      IDHGH  = WWINT(14)
!
      AWG1 = WWAWG(1)
      AWG2 = WWAWG(2)
      AWG3 = WWAWG(3)
      AWG4 = WWAWG(4)
      AWG5 = WWAWG(5)
      AWG6 = WWAWG(6)
      AWG7 = WWAWG(7)
      AWG8 = WWAWG(8)
!
!     *** Initialize auxiliary arrays per gridpoint ***
!
      DO ID = MDC4MI, MDC4MA
        DO IS = MSC4MI, MSC4MA
          UE(IS,ID)   = 0.
          SA1(IS,ID)  = 0.
          SA2(IS,ID)  = 0.
          SFNL(IS,ID) = 0.
        ENDDO
      ENDDO
!
!     *** Calculate prop. constant.                           ***
!     *** Calculate factor R(X) to calculate the NL wave-wave ***
!     *** interaction for shallow water                       ***
!     *** SNLC1 = 1/GRAV**4                                   ***
!
      SNLCS1 = PQUAD(3)                                                   34.00
      SNLCS2 = PQUAD(4)                                                   34.00
      SNLCS3 = PQUAD(5)                                                   34.00
      X      = MAX ( 0.75 * DEP2(KCGRD(1)) * KMESPC , 0.5 )
      X2     = MAX ( -1.E15, SNLCS3*X)
      CONS   = SNLC1 * ( 1. + SNLCS1/X * (1.-SNLCS2*X) * EXP(X2))
      JACOBI = 2. * PI
!
!     *** extend the area with action density at periodic boundaries ***
!
      DO IDDUM = IDLOW, IDHGH
        ID = MOD ( IDDUM - 1 + MDC , MDC ) + 1
        DO IS=1, MSC
          UE (IS,IDDUM) = AC2(ID,IS,KCGRD(1)) * SPCSIG(IS) * JACOBI
        ENDDO
      ENDDO
!
      DO IS = MSC+1, ISHGH
        DO ID = IDLOW, IDHGH
          UE(IS,ID) = UE(IS-1,ID) * FACHFR
        ENDDO
      ENDDO
!
!     *** Calculate (unfolded) interactions ***
!     *** Energy at interacting bins        ***
!
      DO IS = ISCLW, ISCHG
        DO ID = 1, MDC
          E00    =        UE(IS      ,ID      )
          EP1    = AWG1 * UE(IS+ISP1,ID+IDP1) +
     &             AWG2 * UE(IS+ISP1,ID+IDP ) +
     &             AWG3 * UE(IS+ISP ,ID+IDP1) +
     &             AWG4 * UE(IS+ISP ,ID+IDP )
          EM1    = AWG5 * UE(IS+ISM1,ID-IDM1) +
     &             AWG6 * UE(IS+ISM1,ID-IDM ) +
     &             AWG7 * UE(IS+ISM ,ID-IDM1) +
     &             AWG8 * UE(IS+ISM ,ID-IDM )
          EP2    = AWG1 * UE(IS+ISP1,ID-IDP1) +
     &             AWG2 * UE(IS+ISP1,ID-IDP ) +
     &             AWG3 * UE(IS+ISP ,ID-IDP1) +
     &             AWG4 * UE(IS+ISP ,ID-IDP )
          EM2    = AWG5 * UE(IS+ISM1,ID+IDM1) +
     &             AWG6 * UE(IS+ISM1,ID+IDM ) +
     &             AWG7 * UE(IS+ISM ,ID+IDM1) +
     &             AWG8 * UE(IS+ISM ,ID+IDM )
!
!         Contribution to interactions
!
          FACTOR = CONS * AF11(IS) * E00
!
          SA1A   = E00 * ( EP1*DAL1 + EM1*DAL2 ) * CNL4_1(IDIA)
          SA1B   = SA1A - EP1*EM1*DAL3 * CNL4_2(IDIA)
          SA2A   = E00 * ( EP2*DAL1 + EM2*DAL2 ) * CNL4_1(IDIA)
          SA2B   = SA2A - EP2*EM2*DAL3 * CNL4_2(IDIA)
!

          SA1 (IS,ID) = FACTOR * SA1B
          SA2 (IS,ID) = FACTOR * SA2B
!
          IF(ITEST.GE.100 .AND. TESTFL) THEN
            WRITE(PRINTF,9002) E00,EP1,EM1,EP2,EM2
 9002       FORMAT (' E00 EP1 EM1 EP2 EM2  :',5E11.4)
            WRITE(PRINTF,9003) SA1A,SA1B,SA2A,SA2B
 9003       FORMAT (' SA1A SA1B SA2A SA2B  :',4E11.4)
            WRITE(PRINTF,9004) IS,ID,SA1(IS,ID),SA2(IS,ID)
 9004       FORMAT (' IS ID SA1() SA2()    :',2I4,2E12.4)
            WRITE(PRINTF,9005) FACTOR,JACOBI
 9005       FORMAT (' FACTOR JACOBI        : ',2E12.4)
          END IF
!
        ENDDO
      ENDDO
!
!     *** Fold interactions to side angles -> domain in theta is ***
!     *** periodic                                               ***
!
      DO ID = 1, IDHGH - MDC
        ID0   = 1 - ID
        DO IS = ISCLW, ISCHG
          SA1 (IS,MDC+ID) = SA1 (IS,  ID   )
          SA2 (IS,MDC+ID) = SA2 (IS,  ID   )
          SA1 (IS,  ID0 ) = SA1 (IS,MDC+ID0)
          SA2 (IS,  ID0 ) = SA2 (IS,MDC+ID0)
        ENDDO
      ENDDO
!
!     *** Put source term together (To save space I=IS and ***
!     *** J=MDC is used)                                   ***
!
      FAC = 1.                                                            40.17

      DO I = 1, MSC
        SIGPI = SPCSIG(I) * JACOBI                                        30.72
        DO J = 1, MDC
          SFNL(I,J) =   - 2. * ( SA1(I,J) + SA2(I,J) )
     &        + AWG1 * ( SA1(I-ISP1,J-IDP1) + SA2(I-ISP1,J+IDP1) )
     &        + AWG2 * ( SA1(I-ISP1,J-IDP ) + SA2(I-ISP1,J+IDP ) )
     &        + AWG3 * ( SA1(I-ISP ,J-IDP1) + SA2(I-ISP ,J+IDP1) )
     &        + AWG4 * ( SA1(I-ISP ,J-IDP ) + SA2(I-ISP ,J+IDP ) )
     &        + AWG5 * ( SA1(I-ISM1,J+IDM1) + SA2(I-ISM1,J-IDM1) )
     &        + AWG6 * ( SA1(I-ISM1,J+IDM ) + SA2(I-ISM1,J-IDM ) )
     &        + AWG7 * ( SA1(I-ISM ,J+IDM1) + SA2(I-ISM ,J-IDM1) )
     &        + AWG8 * ( SA1(I-ISM ,J+IDM ) + SA2(I-ISM ,J-IDM ) )
!
!         *** store value in auxiliary array and use values in ***
!         *** next four sweeps (see subroutine FILNL3)         ***
!
          IF (IDIA.EQ.1) THEN
            MEMNL4(J,I,KCGRD(1)) = FAC * SFNL(I,J) / SIGPI
          ELSE
            MEMNL4(J,I,KCGRD(1)) = MEMNL4(J,I,KCGRD(1)) +
     &                             FAC * SFNL(I,J) / SIGPI
          END IF
        ENDDO
      ENDDO
!
!     *** test output ***
!
      IF (ITEST .GE. 50 .AND. TESTFL) THEN
        WRITE(PRINTF,*)
        WRITE(PRINTF,*) ' SWSNL4 subroutine '
        WRITE(PRINTF,9011) IDP, IDP1, IDM, IDM1
 9011   FORMAT (' IDP IDP1 IDM IDM1     :',4I5)
        WRITE (PRINTF,9013) ISP, ISP1, ISM, ISM1
 9013   FORMAT (' ISP ISP1 ISM ISM1     :',4I5)
        WRITE (PRINTF,9015) ISLOW, ISHGH, IDLOW, IDHGH
 9015   FORMAT (' ISLOW ISHG IDLOW IDHG :',4I5)
        WRITE(PRINTF,9016) ISCLW, ISCHG, JACOBI
 9016   FORMAT (' ICLW ICHG JACOBI      :',2I5,E12.4)
        WRITE (PRINTF,9017) AWG1, AWG2, AWG3, AWG4
 9017   FORMAT (' AWG1 AWG2 AWG3 AWG4   :',4E12.4)
        WRITE (PRINTF,9018) AWG5, AWG6, AWG7, AWG8
 9018   FORMAT (' AWG5 AWG6 AWG7 AWG8   :',4E12.4)
        WRITE (PRINTF,9019) MSC4MI, MSC4MA, MDC4MI, MDC4MA
 9019   FORMAT (' S4MI S4MA D4MI D4MA   :',4I6)
        WRITE(PRINTF,9020) SNLC1,X,X2,CONS
 9020   FORMAT (' SNLC1  X  X2  CONS    :',4E12.4)
        WRITE(PRINTF,9021) DEP2(KCGRD(1)),KMESPC,FACHFR,PI
 9021   FORMAT (' DEPTH KMESPC FACHFR PI:',4E12.4)
        WRITE(PRINTF,*)
!
!       *** value source term in every bin ***
!
        IF(ITEST.GE. 150 ) THEN
          DO I=1, MSC
            DO J=1, MDC
              WRITE(PRINTF,2006) I,J,MEMNL4(J,I,KCGRD(1)),SFNL(I,J),
     &                           SPCSIG(I)
 2006         FORMAT (' I J MEMNL() SFNL() SPCSIG:',2I4,3E12.4)
            ENDDO
          ENDDO
        END IF
      END IF
!
      DEALLOCATE (SA1, SA2, SFNL, UE)

      RETURN
!
      END SUBROUTINE SWSNL4
!
!*********************************************************************
      SUBROUTINE SWSNL8 (WWINT   ,UE      ,SA1     ,SA2     ,SPCSIG  ,
     &                   SNLC1   ,DAL1    ,DAL2    ,DAL3    ,SFNL    ,
     &                   DEP2    ,AC2     ,KMESPC  ,MEMNL4  ,FACHFR  )
!*********************************************************************
!
      USE SWCOMM3                                                         40.41
      USE SWCOMM4                                                         40.41
      USE OCPCOMM4                                                        40.41
      USE M_SNL4
!
      IMPLICIT NONE
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering                              |
!     | Fluid Mechanics Section                                   |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmers: H.L. Tolman, R.C. Ris                        |
!   --|-----------------------------------------------------------|--
!
!
!     SWAN (Simulating WAves Nearshore); a third generation wave model
!     Copyright (C) 2008  Delft University of Technology
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!     GNU General Public License for more details.
!
!     A copy of the GNU General Public License is available at
!     http://www.gnu.org/copyleft/gpl.html#SEC3
!     or by writing to the Free Software Foundation, Inc.,
!     59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!
!  0. Authors
!
!     40.41: Marcel Zijlema
!
!  1. Updates
!
!     40.41, Sep. 04: piecewise constant interpolation instead
!                     of bi-linear one
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. Purpose
!
!     Calculate non-linear interaction using the discrete interaction
!     approximation (Hasselmann and Hasselmann 1985; WAMDI group 1988)
!     for the full circle.
!
!  3. Method
!
!     Discrete interaction approximation. To make interpolation simple,
!     the interactions are calculated in a "folded" space.
!
!                            Frequencies -->
!                 +---+---------------------+---------+- IDHGH
!              d  | 3 :          2          :    2    |
!              i  + - + - - - - - - - - - - + - - - - +- MDC
!              r  |   :                     :         |
!              e  | 3 :  original spectrum  :    1    |
!              c  |   :                     :         |
!              t. + - + - - - - - - - - - - + - - - - +- 1
!                 | 3 :          2          :    2    |
!                 +---+---------------------+---------+- IDLOW
!                 |   |                     |     ^   |
!              ISLOW  1                    MSC    |   ISHGH
!                     |                           |
!                   ISCLW                        ISCHG
!              lowest discrete               highest discrete
!                central bin                   central bin
!
!                            1 : Extra tail added beyond MSC
!                            2 : Spectrum copied outside ID range
!                            3 : Empty bins at low frequencies
!
!     ISLOW =  1  + ISM1
!     ISHGH = MSC + ISP1 - ISM1
!     ISCLW =  1
!     ISCHG = MSC - ISM1
!     IDLOW =  1  - MAX(IDM1,IDP1)
!     IDHGH = MDC + MAX(IDM1,IDP1)
!
!     Note: using this subroutine requires an additional array
!           with size MXC*MYC*MDC*MSC.
!
!  4. Argument variables
!
!     MCGRD : number of wet grid points of the computational grid
!     MDC   : grid points in theta-direction of computational grid
!     MDC4MA: highest array counter in directional space (Snl4)
!     MDC4MI: lowest array counter in directional space (Snl4)
!     MSC   : grid points in sigma-direction of computational grid
!     MSC4MA: highest array counter in frequency space (Snl4)
!     MSC4MI: lowest array counter in frequency space (Snl4)
!     WWINT : counters for quadruplet interactions
!
      INTEGER WWINT(*)
!
!     AC2   : action density
!     AF11  : scaling frequency
!     DAL1  : coefficient for the quadruplet interactions
!     DAL2  : coefficient for the quadruplet interactions
!     DAL3  : coefficient for the quadruplet interactions
!     DEP2  : depth
!     FACHFR
!     KMESPC: mean average wavenumber over full spectrum
!     MEMNL4
!     PI    : circular constant
!     SA1   : interaction contribution of first quadruplet (unfolded space)
!     SA2   : interaction contribution of second quadruplet (unfolded space)
!     SFNL
!     SNLC1
!     SPCSIG: relative frequencies in computational domain in sigma-space
!     UE    : "unfolded" spectrum
!
      REAL    DAL1, DAL2, DAL3, FACHFR, KMESPC, SNLC1
      REAL    AC2(MDC,MSC,MCGRD)
      REAL    DEP2(MCGRD)
      REAL    MEMNL4(MDC,MSC,MCGRD)
      REAL    SA1(MSC4MI:MSC4MA,MDC4MI:MDC4MA)
      REAL    SA2(MSC4MI:MSC4MA,MDC4MI:MDC4MA)
      REAL    SFNL(MSC4MI:MSC4MA,MDC4MI:MDC4MA)
      REAL    SPCSIG(MSC)
      REAL    UE(MSC4MI:MSC4MA,MDC4MI:MDC4MA)
!
!  7. Common blocks used
!
!
!  8. Subroutines used
!
!     ---
!
!  9. Subroutines calling
!
!     SOURCE (in SWANCOM1)
!
! 12. Structure
!
!     -------------------------------------------
!       Initialisations.
!       Calculate proportionality constant.
!       Prepare auxiliary spectrum.
!       Calculate (unfolded) interactions :
!       -----------------------------------------
!         Energy at interacting bins
!         Contribution to interactions
!         Fold interactions to side angles
!       -----------------------------------------
!       Put source term together
!     -------------------------------------------
!
! 13. Source text
!
!*******************************************************************
!
      INTEGER   IS      ,ID      ,ID0     ,I       ,J       ,
     &          ISLOW   ,ISHGH   ,IDLOW   ,IDHGH   ,
     &          IDPP    ,IDMM    ,ISPP    ,ISMM    ,
     &          ISCLW   ,ISCHG   ,IDDUM   ,IENT
!
      REAL      X       ,X2      ,CONS    ,FACTOR  ,SNLCS1  ,SNLCS2  ,
     &          SNLCS3  ,E00     ,EP1     ,EM1     ,EP2     ,EM2     ,
     &          SA1A    ,SA1B    ,SA2A    ,SA2B    ,
     &          JACOBI  ,SIGPI
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'SWSNL8')
!
      ISLOW  = WWINT(9)
      ISHGH  = WWINT(10)
      ISCLW  = WWINT(11)
      ISCHG  = WWINT(12)
      IDLOW  = WWINT(13)
      IDHGH  = WWINT(14)
      IDPP   = WWINT(21)
      IDMM   = WWINT(22)
      ISPP   = WWINT(23)
      ISMM   = WWINT(24)
!
!     *** Calculate prop. constant.                           ***
!     *** Calculate factor R(X) to calculate the NL wave-wave ***
!     *** interaction for shallow water                       ***
!     *** SNLC1 = 1/GRAV**4                                   ***
!
      SNLCS1 = PQUAD(3)
      SNLCS2 = PQUAD(4)
      SNLCS3 = PQUAD(5)
      X      = MAX ( 0.75 * DEP2(KCGRD(1)) * KMESPC , 0.5 )
      X2     = MAX ( -1.E15, SNLCS3*X)
      CONS   = SNLC1 * ( 1. + SNLCS1/X * (1.-SNLCS2*X) * EXP(X2))
      JACOBI = 2. * PI
!
!     *** extend the area with action density at periodic boundaries ***
!
      DO IDDUM = IDLOW, IDHGH
        ID = MOD ( IDDUM - 1 + MDC , MDC ) + 1
        DO IS=1, MSC
          UE (IS,IDDUM) = AC2(ID,IS,KCGRD(1)) * SPCSIG(IS) * JACOBI
        ENDDO
      ENDDO
!
      DO ID = IDLOW, IDHGH
        DO IS = MSC+1, ISHGH
          UE(IS,ID) = UE(IS-1,ID) * FACHFR
        ENDDO
      ENDDO
!
!     *** Calculate (unfolded) interactions ***
!     *** Energy at interacting bins        ***
!
      DO ID = 1, MDC
        DO IS = ISCLW, ISCHG
          E00 = UE(IS     ,ID     )
          EP1 = UE(IS+ISPP,ID+IDPP)
          EM1 = UE(IS+ISMM,ID-IDMM)
          EP2 = UE(IS+ISPP,ID-IDPP)
          EM2 = UE(IS+ISMM,ID+IDMM)
!
!         Contribution to interactions
!
          FACTOR = CONS * AF11(IS) * PQUAD(2) * E00
!
          SA1A   = E00 * ( EP1*DAL1 + EM1*DAL2 )
          SA1B   = SA1A - EP1*EM1*DAL3
          SA2A   = E00 * ( EP2*DAL1 + EM2*DAL2 )
          SA2B   = SA2A - EP2*EM2*DAL3
!
          SA1 (IS,ID) = FACTOR * SA1B
          SA2 (IS,ID) = FACTOR * SA2B
!
        ENDDO
      ENDDO
!
!     *** Fold interactions to side angles -> domain in theta is ***
!     *** periodic                                               ***
!
      DO ID = 1, IDHGH - MDC
        ID0   = 1 - ID
        DO IS = ISCLW, ISCHG
          SA1 (IS,MDC+ID) = SA1 (IS,  ID   )
          SA2 (IS,MDC+ID) = SA2 (IS,  ID   )
          SA1 (IS,  ID0 ) = SA1 (IS,MDC+ID0)
          SA2 (IS,  ID0 ) = SA2 (IS,MDC+ID0)
        ENDDO
      ENDDO
!
!     *** Put source term together (To save space I=IS and ***
!     *** J=MDC is used)  ----                             ***
!
      DO I = 1, MSC
        SIGPI = SPCSIG(I) * JACOBI
        DO J = 1, MDC
          SFNL(I,J) =   - 2. * ( SA1(I,J) + SA2(I,J) )
     &        + ( SA1(I-ISPP,J-IDPP) + SA2(I-ISPP,J+IDPP) )
     &        + ( SA1(I-ISMM,J+IDMM) + SA2(I-ISMM,J-IDMM) )
!
!         *** store value in auxiliary array and use values in ***
!         *** next four sweeps (see subroutine FILNL3)         ***
!
          MEMNL4(J,I,KCGRD(1)) = SFNL(I,J) / SIGPI
        ENDDO
      ENDDO
!
!     *** value source term in every bin ***
!
      IF ( ITEST.GE.150 .AND. TESTFL ) THEN
         DO I=1, MSC
            DO J=1, MDC
               WRITE(PRINTF,2006) I,J,MEMNL4(J,I,KCGRD(1)),SFNL(I,J),
     &                            SPCSIG(I)
 2006          FORMAT (' I J MEMNL() SFNL() SPCSIG:',2I4,3E12.4)
            ENDDO
         ENDDO
      END IF
!
      RETURN
!
      END SUBROUTINE SWSNL8
!
!*******************************************************************
!
      SUBROUTINE FILNL3 (IDCMIN  ,IDCMAX  ,IMATRA  ,IMATDA  ,AC2     ,    40.23
     &                   MEMNL4  ,PLNL4S  ,ISSTOP                    )    40.41
!
!*******************************************************************
!
      USE SWCOMM3                                                         40.41
      USE SWCOMM4                                                         40.41
      USE OCPCOMM4                                                        40.41
!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering                              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmers: The SWAN team                                |
!   --|-----------------------------------------------------------|--
!
!
!     SWAN (Simulating WAves Nearshore); a third generation wave model
!     Copyright (C) 2008  Delft University of Technology
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!     GNU General Public License for more details.
!
!     A copy of the GNU General Public License is available at
!     http://www.gnu.org/copyleft/gpl.html#SEC3
!     or by writing to the Free Software Foundation, Inc.,
!     59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!
!  0. Authors
!
!     40.23: Marcel Zijlema
!     40.41: Marcel Zijlema
!
!  1. Updates
!
!     40.23, Aug. 02: rhs and main diagonal adjusted according to Patankar-rules
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. Purpose
!
!     Fill the IMATRA/IMATDA arrays with the nonlinear wave-wave interaction
!     source term for a gridpoint ix,iy per sweep direction
!
!  3. Method
!
!
!  4. Argument variables
!
!  7. Common blocks used
!
!
!  8. Subroutines used
!
!     ---
!
!  9. Subroutines calling
!
!     SOURCE (in SWANCOM1)
!
! 12. Structure
!
!     -------------------------------------------
!     Do for every frequency and spectral direction within a sweep
!         fill IMATRA/IMATDA
!     -------------------------------------------
!     End of FILNL3
!     -------------------------------------------
!
! 13. Source text
!
!*******************************************************************
!
      INTEGER   IS      ,ID      ,ISSTOP
!
      REAL      IMATRA(MDC,MSC)           ,
     &          IMATDA(MDC,MSC)           ,                               40.23
     &          AC2(MDC,MSC,MCGRD)        ,                               40.23
     &          PLNL4S(MDC,MSC,NPTST)     ,                               40.00
     &          MEMNL4(MDC,MSC,MCGRD)                                     30.21
!
      INTEGER   IDCMIN(MSC)         ,
     &          IDCMAX(MSC)
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'FILNL3')
!
      DO 990 IS=1, ISSTOP
        DO 980 IDDUM = IDCMIN(IS), IDCMAX(IS)
          ID = MOD ( IDDUM - 1 + MDC , MDC ) + 1
          IF(TESTFL) PLNL4S(ID,IS,IPTST) = MEMNL4(ID,IS,KCGRD(1))         40.00
          IF (MEMNL4(ID,IS,KCGRD(1)).GT.0.) THEN                          40.23
             IMATRA(ID,IS) = IMATRA(ID,IS) + MEMNL4(ID,IS,KCGRD(1))
          ELSE
             IMATDA(ID,IS) = IMATDA(ID,IS) - MEMNL4(ID,IS,KCGRD(1)) /
     &                       MAX(1.E-18,AC2(ID,IS,KCGRD(1)))
          END IF
  980   CONTINUE
  990 CONTINUE
!
      IF ( TESTFL .AND. ITEST.GE.50 ) THEN
        WRITE(PRINTF,9000) IDCMIN(1),IDCMAX(1),MSC,ISSTOP
 9000   FORMAT(' FILNL3: ID_MIN ID_MAX MSC ISTOP :',4I6)
        IF ( ITEST .GE. 100 ) THEN
          DO IS=1, ISSTOP
            DO IDDUM = IDCMIN(IS), IDCMAX(IS)
              ID = MOD ( IDDUM - 1 + MDC , MDC ) + 1
              WRITE(PRINTF,6001) IS,ID,MEMNL4(ID,IS,KCGRD(1))             30.21
 6001         FORMAT(' FILNL3: IS ID MEMNL()          :',2I6,E12.4)
            ENDDO
          ENDDO
        ENDIF
      ENDIF
!
      RETURN
!     End of FILNL3
      END
!
!********************************************************************
!
!     <<  Numerical Computations of the Nonlinear Energy Transfer
!           of Gravity Wave Spectra in Finite Water Depth  >>
!
!                  << developed by Noriaki Hashimoto >>
!
!        References:  N. Hashimoto et al. (1998)
!                     Komatsu and Masuda (2000) (in Japanese)
!
!
!     SWAN (Simulating WAves Nearshore); a third generation wave model
!     Copyright (C) 2008  Delft University of Technology
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!     GNU General Public License for more details.
!
!     A copy of the GNU General Public License is available at
!     http://www.gnu.org/copyleft/gpl.html#SEC3
!     or by writing to the Free Software Foundation, Inc.,
!     59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!
      SUBROUTINE RIAM_SLW(LMAX,N,N2,G,H,DQ,DQ2,DT,DT2,W,P,ACT,SNL,MINT)   40.17

      USE M_SNL4                                                          40.41

!     LMAX
!     N     : number of directional bins                                  40.17
!     N2
!     G     : gravitational acceleration                                  40.17
!     H     : depth                                                       40.17
!     DQ
!     DQ2
!     DT    : size of the directional bins (Delta Theta)                  40.17
!     DT2
!     W     : discretised frequency array                                 40.17
!     P     : density of water                                            40.17
!     ACT   : action density                                              40.17
!     SNL   : Quadruplet source term                                      40.17
!     MINT

!     IW4   : counter for the 4th frequency                               40.17
!     W4    : frequency of the fourth quadruplet component                40.17
!     AK4   : wavenumber of the fourth quadruplet component               40.17
!     DNA   : coefficient in eq. 17 of Hashimoto (98)                     40.17
!             = 2 w4 k4 / Cg(k4)                                          40.17
!     IW3L
!     II
!     JJ
!     DI
!     DJ
!     CGK4  : group velocity for the fourth quadruplet component          40.17

      REAL :: W(LMAX), ACT(LMAX,N), SNL(LMAX,N)                           40.17

!     Initialisation of the quadruplet source term                        40.17

      SNL(:,:) = 0.                                                       40.17
!
!     =================
      DO 130 IW4=1,LMAX
!     =================
!
      W4=W(IW4)

!     WAVE converts nondimensional Sqr(w)d/g to nondimensional kd

      AK4=WAVE(W4**2*H/G)/H

!     CGCMP computes group velocity

      CALL CGCMP(G,AK4,H,CGK4)

!     Calculates the coefficient in equation (17) of Hashimoto (1998)

      DNA=2.*W4*AK4/CGK4

      IW3L=MAX0(1,NINT(IW4-ALOG(3.)/DQ))

!     ------------------------------------------------------------
      CALL PRESET(LMAX,N,IW3L,IW4,N2,G,H,W4,AK4,W,DQ,DQ2,DT,DT2,P)        40.41 40.17
!     ------------------------------------------------------------
!
!     ==============
      DO 130 IT4=1,N
!     ==============
!
!     ================
      CURRIAM => FRIAM                                                    40.41
      DO                !    (W1 - W3 - W4 - W2)                          40.41
!     ================
!
      K1=CURRIAM%II(1)+IW4                                                40.41
      K2=CURRIAM%II(2)+IW4                                                40.41
      K3=CURRIAM%II(3)+IW4                                                40.41
!
      IF(K1.LT.1   ) GO TO 120                                            40.17
      IF(K2.GT.LMAX) GO TO 120                                            40.17
!
      GG=CURRIAM%SSS*DNA                                                  40.41
!
      M1P= CURRIAM%JJ(1)+IT4                                              40.41
      M1N=-CURRIAM%JJ(1)+IT4                                              40.41
      M2P= CURRIAM%JJ(2)+IT4                                              40.41
      M2N=-CURRIAM%JJ(2)+IT4                                              40.41
      M3P= CURRIAM%JJ(3)+IT4                                              40.41
      M3N=-CURRIAM%JJ(3)+IT4                                              40.41
!
      M1P=MOD(M1P,N)
      M1N=MOD(M1N,N)
      M2P=MOD(M2P,N)
      M2N=MOD(M2N,N)
      M3P=MOD(M3P,N)
      M3N=MOD(M3N,N)
!
      IF(M1P.LT.1) M1P=M1P+N
      IF(M1N.LT.1) M1N=M1N+N
      IF(M2P.LT.1) M2P=M2P+N
      IF(M2N.LT.1) M2N=M2N+N
      IF(M3P.LT.1) M3P=M3P+N
      IF(M3N.LT.1) M3N=M3N+N
!
      A4=ACT(IW4,IT4)
!
      IF(MINT.EQ.1) THEN
!
        K01=0
        K02=0
        K03=0
        M01=0
        M02=0
        M03=0
        IF(CURRIAM%DI(1).LE.1.) K01= 1                                    40.41
        IF(CURRIAM%DI(2).LE.1.) K02= 1                                    40.41
        IF(CURRIAM%DI(3).LE.1.) K03= 1                                    40.41
        IF(CURRIAM%DJ(1).LE.1.) M01=-1                                    40.41
        IF(CURRIAM%DJ(2).LE.1.) M02=-1                                    40.41
        IF(CURRIAM%DJ(3).LE.1.) M03=-1                                    40.41
!
        K11=K1+K01
        K21=K2+K02
        K31=K3+K03
        IF(K11.GT.LMAX) K11=LMAX
        IF(K21.GT.LMAX) K21=LMAX
        IF(K31.GT.LMAX) K31=LMAX
!
        K111=K1+K01-1
        K211=K2+K02-1
        K311=K3+K03-1
        IF(K111.LT.1) K111=1
        IF(K211.LT.1) K211=1
        IF(K311.LT.1) K311=1
!
        M1P1=M1P+M01
        M2P1=M2P+M02
        M3P1=M3P+M03
        IF(M1P1.LT.1) M1P1=N
        IF(M2P1.LT.1) M2P1=N
        IF(M3P1.LT.1) M3P1=N
!
        M1N1=M1N-M01
        M2N1=M2N-M02
        M3N1=M3N-M03
        IF(M1N1.GT.N) M1N1=1
        IF(M2N1.GT.N) M2N1=1
        IF(M3N1.GT.N) M3N1=1
!
        M1P11=M1P+M01+1
        M2P11=M2P+M02+1
        M3P11=M3P+M03+1
        IF(M1P11.GT.N) M1P11=1
        IF(M2P11.GT.N) M2P11=1
        IF(M3P11.GT.N) M3P11=1
!
        M1N11=M1N-M01-1
        M2N11=M2N-M02-1
        M3N11=M3N-M03-1
        IF(M1N11.LT.1) M1N11=N
        IF(M2N11.LT.1) M2N11=N
        IF(M3N11.LT.1) M3N11=N
!
        DI1=CURRIAM%DI(1)                                                 40.41
        DI2=CURRIAM%DI(2)                                                 40.41
        DI3=CURRIAM%DI(3)                                                 40.41
        DJ1=CURRIAM%DJ(1)                                                 40.41
        DJ2=CURRIAM%DJ(2)                                                 40.41
        DJ3=CURRIAM%DJ(3)                                                 40.41
!
        A1P=(DI1*(DJ1*ACT(K11, M1P1)+ACT(K11, M1P11))                     40.41
     &          + DJ1*ACT(K111,M1P1)+ACT(K111,M1P11))                     40.41
     &          /((1.+DI1)*(1.+DJ1))                                      40.41
!
        A1N=(DI1*(DJ1*ACT(K11, M1N1)+ACT(K11, M1N11))                     40.41
     &          + DJ1*ACT(K111,M1N1)+ACT(K111,M1N11))                     40.41
     &          /((1.+DI1)*(1.+DJ1))                                      40.41
!
        A2P=(DI2*(DJ2*ACT(K21, M2P1)+ACT(K21, M2P11))                     40.41
     &          + DJ2*ACT(K211,M2P1)+ACT(K211,M2P11))                     40.41
     &          /((1.+DI2)*(1.+DJ2))                                      40.41
!
        A2N=(DI2*(DJ2*ACT(K21, M2N1)+ACT(K21, M2N11))                     40.41
     &          + DJ2*ACT(K211,M2N1)+ACT(K211,M2N11))                     40.41
     &          /((1.+DI2)*(1.+DJ2))                                      40.41
!
        A3P=(DI3*(DJ3*ACT(K31, M3P1)+ACT(K31, M3P11))                     40.41
     &          + DJ3*ACT(K311,M3P1)+ACT(K311,M3P11))                     40.41
     &          /((1.+DI3)*(1.+DJ3))                                      40.41
!
        A3N=(DI3*(DJ3*ACT(K31, M3N1)+ACT(K31, M3N11))                     40.41
     &          + DJ3*ACT(K311,M3N1)+ACT(K311,M3N11))                     40.41
     &          /((1.+DI3)*(1.+DJ3))                                      40.41
!
      ELSE
!
        IF(K1.LT.1) K1=1
        IF(K2.LT.1) K2=1
        IF(K3.LT.1) K3=1
!
        A1P=ACT(K1,M1P)
        A1N=ACT(K1,M1N)
        A2P=ACT(K2,M2P)
        A2N=ACT(K2,M2N)
        A3P=ACT(K3,M3P)
        A3N=ACT(K3,M3N)
!
      ENDIF
!
      W1P2P=A1P+A2P
      W1N2N=A1N+A2N
      S1P2P=A1P*A2P
      S1N2N=A1N*A2N
      W3P4 =A3P+A4
      W3N4 =A3N+A4
      S3P4 =A3P*A4
      S3N4 =A3N*A4
!
      XP=S1P2P*W3P4-S3P4*W1P2P
      XN=S1N2N*W3N4-S3N4*W1N2N
!
      SNL( K1,M1P)=SNL( K1,M1P)-XP*GG
      SNL( K1,M1N)=SNL( K1,M1N)-XN*GG
      SNL( K2,M2P)=SNL( K2,M2P)-XP*GG
      SNL( K2,M2N)=SNL( K2,M2N)-XN*GG
      SNL( K3,M3P)=SNL( K3,M3P)+XP*GG
      SNL( K3,M3N)=SNL( K3,M3N)+XN*GG
      SNL(IW4,IT4)=SNL(IW4,IT4)+(XP+XN)*GG
!
! ============
  120 CONTINUE
      IF (.NOT.ASSOCIATED(CURRIAM%NEXTRIAM)) EXIT                         40.41
      CURRIAM => CURRIAM%NEXTRIAM                                         40.41
      END DO                                                              40.41
! ============
!
! ============
  130 CONTINUE
! ============
!
      RETURN
      END

      SUBROUTINE PRESET(LMAX,N,IW3L,IW4,N2,G,H,W4,AK4,W,DQ,DQ2,DT,DT2,P)  40.41 40.17

!     T1A   : Angle between theta_1 and theta_a
!     T2A   : Angle between theta_2 and theta_a
!     T34   : Angle between theta_3 and theta_4
!     WA    : wa = w1 + w2 = w3 + w4 (resonance conditions)
!     AKA   : absolute value of ka = k1 + k2 = k3 + k4 (resonance conditions)
!     TA    : Angle theta_a representing vector ka

      REAL :: W(LMAX)                                                     40.17
!
!     ================
      DO 120 IT34=1,N2
!     ================
!
      T34=DT*(IT34-1)
      DT3=DT
      IF(IT34.EQ.1.OR.IT34.EQ.N2) DT3=DT2
!
!     ===================
      DO 110 IW3=IW3L,IW4
!     ===================
!
      W3=W(IW3)
      AK3=WAVE(W3**2*H/G)/H
!
      WA=W3+W4
      AKA=SQRT(AK3*AK3+AK4*AK4+2.*AK3*AK4*COS(T34))
      TA=ATAN2(AK3*SIN(T34),AK3*COS(T34)+AK4)
      R=SQRT(G*AKA*TANH(AKA*H/2.))/WA-1./SQRT(2.)
!
      TL=0.
      ACS=AKA/(2.*WAVE(WA**2*H/(4.*G))/H)
      ACS=SIGN(1.,ACS)*MIN(1.,ABS(ACS))                                   40.41
      IF(R.LT.0.) TL=ACOS(ACS)
      IT1S=NINT(TL/DT)+1
!
!     ===================
      DO 100 IT1A=IT1S,N2
!     ===================
!
      T1A=DT*(IT1A-1)
!
      DT1=DT
      IF(IT1A.EQ.1.OR.IT1A.EQ.N2) DT1=DT2
!
!     -----------------------------------------------------------------
      CALL FINDW1W2(LMAX,IW3,W,G,H,W1,W2,W3,WA,AK1,AK2,AKA,T1A,T2A,IND)   40.17
!     -----------------------------------------------------------------
      IF(IND.LT.0) GO TO 100
!
      IF(W1.LE.W3.AND.W3.LE.W4.AND.W4.LE.W2) THEN
        CALL KERNEL(N,G,H,W1,W2,W3,W4,AK1,AK2,AK3,AK4,AKA,
     &              T1A,T2A,T34,TA,DQ,DT,DT1,DT3,P)                       40.41 40.17
      ENDIF
!
! ============
  100 CONTINUE
  110 CONTINUE
  120 CONTINUE
! ============
!
      RETURN
      END

      SUBROUTINE FINDW1W2(LMAX,                                           40.17
     &                    IW3,W,G,H,W1,W2,W3,WA,AK1,AK2,AKA,T1A,T2A,IND)  40.17
      REAL  :: W(LMAX)                                                    40.17
!
      IND=0
      EPS=0.0005
!
      X1=0.005
      X2=W3
!
!     ---------------------------------------
  110 CALL FDW1(G,H,AKA,T1A,WA,X1,X2,X,EPS,M)
!     ---------------------------------------
!
      IF(M.EQ.0) GO TO 999
!
      W1=X
!
      AK1=WAVE(W1**2*H/G)/H
      AK2=SQRT(AKA*AKA+AK1*AK1-2.*AKA*AK1*COS(T1A))
      W2=SQRT(G*AK2*TANH(AK2*H))
      IF(W1.GT.W2) GO TO 999
      T2A=ATAN2(-AK1*SIN(T1A),AKA-AK1*COS(T1A))
      RETURN
!
  999 IND=-999
      RETURN
      END

      SUBROUTINE FDW1(G,H,AKA,T1A,WA,X1,X2,X,EPS,M)
!
      M=0
      F1=FUNCW1(G,H,X1,AKA,T1A,WA)
      F2=FUNCW1(G,H,X2,AKA,T1A,WA)
   20 IF(F1*F2) 18,19,21
   18 M=M+1
      X=X2-((X2-X1)/(F2-F1)*F2)
      IF((ABS(X-X2)/ABS(X))-EPS) 22,23,23
   23 IF((ABS(X-X1)/ABS(X))-EPS) 22,24,24
   22 RETURN
   24 F=FUNCW1(G,H,X,AKA,T1A,WA)
      FM=F*F1
      IF(FM) 25,22,26
   25 X2=X
      F2=F
      GO TO 20
   26 X1=X
      F1=F
      GO TO 20
   19 M=M+1
      IF(F1) 17,16,17
   16 X=X1
      GO TO 22
   17 X=X2
      GO TO 22
   21 M=0
      RETURN
      END

      FUNCTION FUNCW1(G,H,X,AKA,T1A,WA)
!
      AK1=WAVE(X**2*H/G)/H
      AK2=SQRT(AKA*AKA+AK1*AK1-2.*AKA*AK1*COS(T1A))
      FUNCW1=WA-SQRT(G*AK1*TANH(AK1*H))-SQRT(G*AK2*TANH(AK2*H))
!
      RETURN
      END

      SUBROUTINE KERNEL(N,G,H,W1,W2,W3,W4,AK1,AK2,AK3,AK4,AKA,
     &                  T1A,T2A,T34,TA,DQ,DT,DT1,DT3,P)                   40.41 40.17

      USE M_SNL4                                                          40.41
!
      LOGICAL, SAVE :: LRIAM = .FALSE.                                    40.41
      TYPE(RIAMDAT), POINTER :: RIAMTMP                                   40.41
!
      PI =4.*ATAN(1.)
      PI2=2.*PI
!
!     ------------------------
      CALL CGCMP(G,AK1,H,CGK1)
      CALL CGCMP(G,AK2,H,CGK2)
      CALL CGCMP(G,AK3,H,CGK3)
!     ------------------------

!     Calculate denominator S (eq. 14, Hashimoto (1998))

      SS0=ABS(1.+CGK2/CGK1*(AK1-AKA*COS(T1A))/AK2)
      IF(SS0.LE.1.E-7) RETURN

!                        k1 k3 w3     G
!     Kernel function  -------------  - dth3 dth1 dOm    (eq. 17)
!                      Cg(k1) Cg(k3)  S

!     where S defined in eq. 14 and G defined in eq. 5 (Hashimoto, 1998)

      CF=AK1*AK3*W3/(CGK1*CGK3)/SS0
     &   *9.*PI*G*G/(4.*P*P*W1*W2*W3*W4)*DT1*DT3*DQ

!     --------------------------------------
      CALL NONKDD(G,H,W3,W4,-W1,AK3,AK4,AK1,
     &                T34,0.,T1A+TA+PI,DDD1)
!     --------------------------------------
      GGG1=CF*DDD1*DDD1
!
!     --------------------------------------
      CALL NONKDD(G,H,W3,W4,-W1,AK3,AK4,AK1,
     &               T34,0.,-T1A+TA+PI,DDD2)
!     --------------------------------------
      GGG2=CF*DDD2*DDD2
!
      I1=1-INTW(W1/W4,DQ)
      I2=1-INTW(W2/W4,DQ)
      I3=1-INTW(W3/W4,DQ)
!
      DI1=WEIGW(W1/W4,DQ)
      DI2=WEIGW(W2/W4,DQ)
      DI3=0.
!
!     ----------------------------
!
      TT1=T1A+TA
      TT2=T2A+TA
      TT3=T34
!
      ALLOCATE(RIAMTMP)                                                   40.41
!
      RIAMTMP%SSS  =GGG1                                                  40.41
      RIAMTMP%II(1)=I1                                                    40.41
      RIAMTMP%II(2)=I2                                                    40.41
      RIAMTMP%II(3)=I3                                                    40.41
      RIAMTMP%DI(1)=DI1                                                   40.41
      RIAMTMP%DI(2)=DI2                                                   40.41
      RIAMTMP%DI(3)=DI3                                                   40.41
      RIAMTMP%JJ(1)=INTT(TT1,DT,PI2,N)-1                                  40.41
      RIAMTMP%JJ(2)=INTT(TT2,DT,PI2,N)-1                                  40.41
      RIAMTMP%JJ(3)=INTT(TT3,DT,PI2,N)-1                                  40.41
      RIAMTMP%DJ(1)=WEIGT(TT1,DT,PI2,N)                                   40.41
      RIAMTMP%DJ(2)=WEIGT(TT2,DT,PI2,N)                                   40.41
      RIAMTMP%DJ(3)=0.                                                    40.41
!
      NULLIFY(RIAMTMP%NEXTRIAM)                                           40.41
      IF ( .NOT.LRIAM ) THEN                                              40.41
         FRIAM = RIAMTMP                                                  40.41
         CURRIAM => FRIAM                                                 40.41
         LRIAM = .TRUE.                                                   40.41
      ELSE
         CURRIAM%NEXTRIAM => RIAMTMP                                      40.41
         CURRIAM => RIAMTMP                                               40.41
      END IF
!
!     ----------------------------
!
      TT1=-T1A+TA
      TT2=-T2A+TA
      TT3= T34
!
      ALLOCATE(RIAMTMP)                                                   40.41
!
      RIAMTMP%SSS  =GGG2                                                  40.41
      RIAMTMP%II(1)=I1                                                    40.41
      RIAMTMP%II(2)=I2                                                    40.41
      RIAMTMP%II(3)=I3                                                    40.41
      RIAMTMP%DI(1)=DI1                                                   40.41
      RIAMTMP%DI(2)=DI2                                                   40.41
      RIAMTMP%DI(3)=DI3                                                   40.41
      RIAMTMP%JJ(1)=INTT(TT1,DT,PI2,N)-1                                  40.41
      RIAMTMP%JJ(2)=INTT(TT2,DT,PI2,N)-1                                  40.41
      RIAMTMP%JJ(3)=INTT(TT3,DT,PI2,N)-1                                  40.41
      RIAMTMP%DJ(1)=WEIGT(TT1,DT,PI2,N)                                   40.41
      RIAMTMP%DJ(2)=WEIGT(TT2,DT,PI2,N)                                   40.41
      RIAMTMP%DJ(3)=0.                                                    40.41
!
      NULLIFY(RIAMTMP%NEXTRIAM)                                           40.41
      IF ( .NOT.LRIAM ) THEN                                              40.41
         FRIAM = RIAMTMP                                                  40.41
         CURRIAM => FRIAM                                                 40.41
         LRIAM = .TRUE.                                                   40.41
      ELSE
         CURRIAM%NEXTRIAM => RIAMTMP                                      40.41
         CURRIAM => RIAMTMP                                               40.41
      END IF
!
      RETURN
      END

      FUNCTION INTW(W,DQ)
      AA=ALOG(W)/DQ
      INTW=1+NINT(-AA)
      RETURN
      END

      FUNCTION WEIGW(W,DQ)
      AA=ALOG(W)/DQ
      L=1+INT(-AA)
      A=DQ*(1-L)-ALOG(W)
      IF(A.NE.0.) THEN
        B=ALOG(W)+DQ*L
        WEIGW=B/A
      ELSE
        WEIGW=1000.
      ENDIF
      RETURN
      END

      FUNCTION INTT(T,DT,PI2,MM)
      T=AMOD(T,PI2)
      IF(T.LT.0.) T=T+PI2
      INTT=NINT(T/DT)+1
      IF(INTT.GT.MM) INTT=INTT-MM
      IF(INTT.LT.1)  INTT=INTT+MM
      RETURN
      END

      FUNCTION WEIGT(T,DT,PI2,MM)
      T=AMOD(T,PI2)
      IF(T.LT.0.) T=T+PI2
      N=INT(T/DT)+1
      IF(N.GT.MM) N=N-MM
      IF(N.LT.1)  N=N+MM
      C=T-DT*(N-1)
      IF(C.NE.0.) THEN
        D=DT*N-T
        WEIGT=D/C
      ELSE
        WEIGT=1000.
      ENDIF
      RETURN
      END

      SUBROUTINE NONKDD(G,H,W1,W2,W3,AK1,AK2,AK3,T1,T2,T3,DDD)

!     See Herterich and Hasselmann (1980; p. 223, eq. B1)

      CALL NONKD(G,H,W1,W2,W3,AK1,AK2,AK3,T1,T2,T3,DD1)
      CALL NONKD(G,H,W2,W3,W1,AK2,AK3,AK1,T2,T3,T1,DD2)
      CALL NONKD(G,H,W3,W1,W2,AK3,AK1,AK2,T3,T1,T2,DD3)
!
      DDD=(DD1+DD2+DD3)/3.
!
      RETURN
      END

      SUBROUTINE NONKD(G,H,W1,W2,W3,AK1,AK2,AK3,T1,T2,T3,DDD)

!     GG    : square of acceleration of gravity (=g**2)
!     G2    : twice the square of acceleration of gravity (=2g**2)
!     WP23  : sum of w2 and w3 (=w2+w3)
!     W123  : sum of w1, w2 and w3 (=w1+w2+w3)
!     AKX1  : size of the x-component of wavenumber k1
!     AKY1  : size of the y-component of wavenumber k1
!     AKX2  : size of the x-component of wavenumber k2
!     AKY2  : size of the y-component of wavenumber k2
!     AKX3  : size of the x-component of wavenumber k3
!     AKY3  : size of the y-component of wavenumber k3
!     AK23  : dot product of wavenumbers k2 and k3
!     W23   : frequency corresponding to wavenumber k2+k3


      DDD=0.
      GG=G*G
      G2=2.*G*G
!
      AKX1=AK1*COS(T1)
      AKY1=AK1*SIN(T1)
!
      AKX2=AK2*COS(T2)
      AKY2=AK2*SIN(T2)
!
      AKX3=AK3*COS(T3)
      AKY3=AK3*SIN(T3)
!
      AK23=SQRT((AKX2+AKX3)**2+(AKY2+AKY3)**2)
      W23=SQRT(G*AK23*TANH(AK23*H))
!
      WP23=W2+W3
      CF0=W23**2-WP23**2
      IF(CF0.EQ.0.) RETURN
!
      W123=W1+W2+W3
      AKXY23=AKX2*AKX3+AKY2*AKY3
      AKXY123=AKX1*(AKX2+AKX3)+AKY1*(AKY2+AKY3)
      CF1=0.
      IF(AK1*H.LE.10.) CF1=AK1*AK1/COSH(AK1*H)**2
!
      CF2=0.
      CF3=0.
      IF(AK3*H.LE.10.) CF2=W2*AK3*AK3/COSH(AK3*H)**2
      IF(AK2*H.LE.10.) CF3=W3*AK2*AK2/COSH(AK2*H)**2

!     See Herterich and Hasselmann (1980), first eq. after eq. B2
!     Equation is divided by I, to avoid complex numbers

      DD=WP23*(AK2*AK3*TANH(AK2*H)*TANH(AK3*H)-AKXY23)-(CF2+CF3)/2.

!     See Herterich and Hasselmann (1980), second eq. after eq. B2

      EE=(AKXY23-W2*W3*(W2**2+W3**2+W2*W3)/GG)/(2.*G)

      CF4=0.
      IF(AK23*H.LE.10.) CF4=W1*AK23**2/COSH(AK23*H)**2

!     Equation B1 of Herterich and Hasselmann (1980), where DD is
!     multiplied by I, to compensate for the earlier division

      DDD=-DD/CF0*(2.*W123*(W1**2*W23**2/GG-AKXY123)-CF4-WP23*CF1)
     &    +DD*W1*(W1*W1+W23*W23)/GG
     &    +EE*(W1**3*WP23/G-G*AKXY123-G*CF1)
     &    +W1*AKXY23*(W123*(W2**2+W3**2)+W2*W3*WP23)/G2
     &    -W1*W2*W2*AK3*AK3*(W123+W3)/G2
     &    -W1*W3*W3*AK2*AK2*(W123+W2)/G2
!
      RETURN
      END

      SUBROUTINE CGCMP(G,AK,H,CG)

!     Calculates group velocity Cg based on depth and wavenumber
!     Includes a deep water limit for kd > 10.

!     G     : gravitational acceleration
!     AK    : wave number
!     H     : depth
!     CG    : group velocity
!     AKH   : depth x Wave number (kd)
!     RGK   : square root of gk
!     SECH2 : the square of the secant Hyperbolic of kd
!             = 1 - TANH(kd)**2

!     Calculation of group velocity Cg

      AKH=AK*H
      IF(AKH.LE.10.) THEN                                                 40.33

!     Shallow water:

!     Cg = (1/2) (1 + (2 kd) / Sinh (2 kd)) (w / k)
!        = (1/2) (1 + (kd SECH2) / Tanh (kd)) (w / k)
!        = (1/2) (1 + (kd SECH2) / Tanh (kd)) (Sqrt (gk Tanh(kd)) / k)
!        = (1/2) (Sqrt (gk) (Sqrt(Tanh (kd)) + kd SECH2 / Sqrt (Tanh (kd)) / k
!        = (1/2) (Sqrt (k) g (Tanh (kd) + kd SECH2) / (k Sqrt (g Tanh (kd)))
!        = (g Tanh(kd) + gkd SECH2) / (2 Sqrt(gk Tanh(kd))

        SECH2=1.-TANH(AKH)**2
        CG=(G*AKH*SECH2+G*TANH(AKH))/(2.*SQRT(G*AK*TANH(AKH)))

      ELSE                                                                40.17

!     Deep water:

!     Cg = w / (2 k)
!        = g / (2 Sqrt (gk))

        RGK=SQRT(G*AK)
        CG=G/(2.*RGK)

      ENDIF                                                               40.17
!
      RETURN

      END SUBROUTINE CGCMP

      SUBROUTINE DCGCMP(G,AK,H,DCG)
!
      AKH=AK*H
      IF(AKH.GT.10.) GO TO 100
      SECH2=1-TANH(AKH)**2
      SECH3=SECH2*SQRT(SECH2)
      DCG=G*H*SECH3*(COSH(AKH)-AKH*SINH(AKH))/SQRT(G*AK*TANH(AKH))
     &   -(G*AKH*SECH2+G*TANH(AKH))**2/(4.*SQRT(G*AK*TANH(AKH))**3)
      RETURN
!
  100 RGK=SQRT(G*AK)
      DCG=-G*G/(4.*RGK**3)
!
      RETURN
      END

      REAL FUNCTION WAVE(D)

!     Transforms nondimensional Sqr(w)d/g into nondimensional kd using the
!     dispersion relation and an iterative method

      IF(D-10.) 2,2,1
    1 XX=D
      GO TO 6
    2 IF(D-1.) 3,4,4
    3 X=SQRT(D)
      GO TO 5
    4 X=D
    5 COTHX=1.0/TANH(X)
      XX=X-(X-D*COTHX)/(1.+D*(COTHX**2-1.))
      E=1.-XX/X
      X=XX
      IF(ABS(E)-0.0005) 6,5,5
    6 WAVE=XX
      RETURN
      END
!
!----------------------------------------------------------------------
      SUBROUTINE SWINTFXNL ( ASWAN,SIGMA,DIR,NDIR,NSIG,NGRID,DEPTH,
     &                       IQTYPE,SNL,KCGRD,ICMAX,IERROR )
!----------------------------------------------------------------------
!
!   +-------+    ALKYON Hydraulic Consultancy & Research
!   |       |    Gerbrant Ph. van Vledder
!   |   +---+
!   |   | +---+  Last update:  9 September 2002
!   +---+ |   |  Release: 5.0
!         +---+
!
!
!     SWAN (Simulating WAves Nearshore); a third generation wave model
!     Copyright (C) 2008  Delft University of Technology
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!     GNU General Public License for more details.
!
!     A copy of the GNU General Public License is available at
!     http://www.gnu.org/copyleft/gpl.html#SEC3
!     or by writing to the Free Software Foundation, Inc.,
!     59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!
      USE M_PARALL
      USE serv_xnl4v5
      USE m_xnldata
!
      IMPLICIT NONE
!-----------------------------------------------------------------------------------
!
!  0. Update history
!
!     Date Modification
!
!     25/02/1999 Initial version
!     16/03/1999 Interface with SWAN updated
!     24/03/1999 Interface extended with KCGRD and ICMAX
!     22/06/1999 Parameter IERR added to interface
!     20/07/1999 Call to Q_QMAIN modified
!     25/07/1999 Output files updated
!     24/09/1999 Finite depth effects included
!     25/11/1999 Bug fixed in deleting files
!     27/12/1999 Interface extended with grav,rho and ftail
!     02/02/2001 Interface modified, was N(k), now A(sigma) like SWAN
!     06/11/2001 Bug fixed in initialisation of Snl
!     08/08/2002 Version 4
!     22/08/2002 Size of direction array modified to conform with SWAN 40.11
!     09/09/2002 Release 5
!     18/05/2004 Implemented in SWAN 40.41, with adapted values of IQTYPE
!
!
!  1. Purpose:
!
!     interface with SWAN model to compute nonlinear transfer with
!     the XNL method for given action density spectrum
!
!  2. Method
!
!     Resio/Tracy deep water geometric scaling
!     Rewritten by Gerbrant van Vledder
!
!  3. Parameter list:
!
! Type     I/O         Name      Description
!----------------------------------------------------------------------------------
      INTEGER, INTENT(IN) :: NDIR                    ! number of directions
      INTEGER, INTENT(IN) :: NSIG                    ! number of sigma values in array
      INTEGER, INTENT(IN) :: NGRID                   ! number of sea points in the comp. grid
      INTEGER, INTENT(IN) :: IQTYPE                  ! method of computing nonlinear quadruplet
!                                                      interactions
      REAL   , INTENT(IN) :: ASWAN(NDIR,NSIG,NGRID)  ! action density spectrum as a
!                                                    ! function of (sigma,dir)
      REAL   , INTENT(IN) :: SIGMA(NSIG)             ! Intrinsic frequencies
      REAL   , INTENT(IN) :: DIR(NDIR,6)             ! directions in radians (when second index =1)
      REAL   , INTENT(IN) :: DEPTH(NGRID)            ! depth array
      INTEGER, INTENT(IN) :: ICMAX                   ! number of points of the computational stencil
      INTEGER, INTENT(IN) :: KCGRD(ICMAX)            ! grid addresses for points of computational stencil
      REAL   , INTENT(OUT):: SNL(NDIR,NSIG,NGRID)    ! nonlinear quadruplet interaction computed with
!                                                    ! a certain exact method (sigma,dir)
      INTEGER, INTENT(OUT):: IERROR                  ! Error indicator. If no errors are detected IERR=0
!--------------------------------------------------------------------------------------------------
!
!  4. Error messages
!
!     An error message is produced within the QUAD system.
!     If no errors are detected IERROR=0
!     1, incorrect IQUAD
!     2, depth < 0
!
!  5. Subroutines calling
!
!     SOURCE
!
!  6. Subroutines used
!
!  7. Remarks
!
!     The SWAN spectrum is given as an action density spectrum
!     as a function of Sigma and Theta: ASWAN(itheta,isig)
!
!     SWINTFXNL is called for each active grid point in a stencil
!     and for each time the complete array with all grid points
!     is given. Related grid points are specified in the array
!     KCGRD and ICMAX.
!
!  8. Structure
!
!  9. Switches
!
! 10. Source code
!------------------------------------------------------------------------------
!     Local parameters
!
      INTEGER ISIG,IDIR            ! counters
      INTEGER IGRID                ! grid index
      INTEGER IQUAD                ! type of computational method for Xnl
!
!     --- assign arrays for intermediate storage of results
!
      REAL AQUAD(NSIG,NDIR)        ! action density spectrum A(sigma,dir)
      REAL XNL(NSIG,NDIR)          ! transfer rate dA/dt(sigma,dir)
      REAL DIAG(NSIG,NDIR)         ! diagonal term (dXnl/dA)
      REAL DIRR(NDIR)              ! single array with directions in radians
!------------------------------------------------------------------------------
!
!     --- initialisations
!
      IERROR  = 0
!
      IGRID   = KCGRD(1) ! set index of current grid index
      DIRR(:) = DIR(:,1) ! copy radian directions to single array
!
      SNL(:,:,IGRID)  = 0.
      DIAG            = 0.
!
!     --- check value of iquad
!
      IF ((IQTYPE.NE.51).AND.(IQTYPE.NE.52).AND.(IQTYPE.NE.53)) THEN
         IERROR = 1
         GOTO 9999
      END IF
!
!  N.B. Take care of different order of indices in QUAD compared to SWAN
!
!     --- compute nonlinear interactions per individual spectrum
!         switch order of indices
!
      DO ISIG = 1, NSIG
        DO IDIR = 1, NDIR
           AQUAD(ISIG,IDIR) = ASWAN(IDIR,ISIG,IGRID)
        END DO
      END DO
!
!     --- transform parameter iquad of SWAN to the parameter iq_quad as
!         needed by the QUAD suite
!
      IQUAD = IQTYPE - 50
!
!     IQTYPE/IQUAD = 51/1   ! deep water transfer
!     IQTYPE/IQUAD = 52/2   ! deep water transfer with WAM depth scaling
!     IQTYPE/IQUAD = 53/3   ! finite depth transfer
!
!     --- call of main subroutine to compute nonlinear quadruplet interactions
!         for a given action density spectrum on a given spectral grid
!
      XNL   = 0.
      DIAG  = 0.
      CALL XNL_MAIN( AQUAD,SIGMA,DIRR,NSIG,NDIR,DEPTH(IGRID),IQUAD,
     &               XNL,DIAG,INODE,IERROR )
!
      IF (IERROR.NE.0) GOTO 9999
!
!     --- convert nonlinear transfer to SWAN convention, only sequence of indices
!
      DO ISIG = 1, NSIG
        DO IDIR = 1, NDIR
           SNL(IDIR,ISIG,IGRID) = XNL(ISIG,IDIR)
        END DO
      END DO
!
 9999 CONTINUE
!
      RETURN
!
      END SUBROUTINE SWINTFXNL
!
!****************************************************************
!
      SUBROUTINE SWLTA ( AC2   , DEP2  , CGO   , SPCSIG,
     &                   KWAVE , IMATRA, IMATDA,
     &                   IDDLOW, IDDTOP, ISSTOP, IDCMIN, IDCMAX,
     &                   HS    , SMEBRK, PLTRI , URSELL )
!
!****************************************************************
!
      USE OCPCOMM4
      USE SWCOMM3
      USE SWCOMM4
!
      IMPLICIT NONE
!
!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering                              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmers: The SWAN team                                |
!   --|-----------------------------------------------------------|--
!
!
!     SWAN (Simulating WAves Nearshore); a third generation wave model
!     Copyright (C) 2008  Delft University of Technology
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!     GNU General Public License for more details.
!
!     A copy of the GNU General Public License is available at
!     http://www.gnu.org/copyleft/gpl.html#SEC3
!     or by writing to the Free Software Foundation, Inc.,
!     59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!
!  0. Authors
!
!     40.55: Marcel Zijlema
!
!  1. Updates
!
!     40.55, Feb. 06: New subroutine
!
!  2. Purpose
!
!     In this subroutine the triad-wave interactions are calculated
!     with the Lumped Triad Approximation of Eldeberky (1996). His
!     expression is based on a parametrization of the biphase (as
!     function of the Ursell number), is directionally uncoupled and
!     takes into account for the self-self interactions only.
!
!     For a full description of the equations reference is made
!     to PhD thesis of Eldeberky (1996). Here only the main expressions
!     are given.
!
!  3. Method
!
!     The parametrized biphase is given by (see eq. 3.19):
!
!                                  0.2
!     beta = - pi/2 + pi/2 tanh ( ----- )
!                                   Ur
!
!     The Ursell number is calculated in routine SINTGRL
!
!     The source term as function of frequency p is (see eq. 7.25):
!
!             +      -
!     S(p) = S(p) + S(p)
!
!     in which
!
!      +
!     S(p) = alpha Cp Cg,p (R(p/2,p/2))**2 sin (|beta|) ( E(p/2)**2 -2 E(p) E(p/2) )
!
!      -          +
!     S(p) = - 2 S(2p)
!
!     with alpha a tunable coefficient and R(p/2,p/2) is the interaction
!     coefficient of which the expression can be found in Eldeberky (1996);
!     see eq. 7.26.
!
!     Note that a slightly adapted formulation of the LTA is used in
!     in the SWAN model:
!
!     - Only positive contributions to higher harmonics are considered
!       here (no energy is transferred to lower harmonics).
!
!     - The mean frequency in the expression of the Ursell number
!       is calculated according to the first order moment over the
!       zeroth order moment (personal communication, Y.Eldeberky, 1997).
!
!     - The interactions are calculated up to 2.5 times the mean
!       frequency only.
!
!     - Since the spectral grid is logarithmically distributed in frequency
!       space, the interactions between central bin and interacting bin
!       are interpolated such that the distance between these bins is
!       factor 2 (nearly).
!
!     - The interactions are calculated in terms of energy density
!       instead of action density. So the action density spectrum
!       is firstly converted to the energy density grid, then the
!       interactions are calculated and then the spectrum is converted
!       to the action density spectrum back.
!
!     - To ensure numerical stability the Patankar rule is used.
!
!  4. Argument variables
!
!     AC2         action density
!     CGO         group velocity
!     DEP2        water depth
!     HS          significant wave height
!     IDCMIN      minimum counter in directional space
!     IDCMAX      maximum counter in directional space
!     IDDLOW      minimum direction that is propagated within a sweep
!     IDDTOP      maximum direction that is propagated within a sweep
!     IMATDA      main diagonal of the linear system
!     IMATRA      right-hand side of system of equations
!     ISSTOP      maximum frequency counter in a sweep
!     KWAVE       wave number
!     PLTRI       triad contribution in TEST points
!     SMEBRK      average (angular) frequency
!     SPCSIG      relative frequencies in computational domain in sigma-space
!     URSELL      Ursell number
!
      INTEGER IDDLOW, IDDTOP, ISSTOP
      INTEGER IDCMIN(MSC), IDCMAX(MSC)

      REAL :: HS, SMEBRK
      REAL :: AC2(MDC,MSC,MCGRD)

      REAL :: CGO(MSC,MICMAX)
      REAL :: DEP2(MCGRD)
      REAL :: IMATDA(MDC,MSC), IMATRA(MDC,MSC)
      REAL :: SPCSIG(MSC)
      REAL :: KWAVE(MSC,MICMAX)
      REAL :: PLTRI(MDC,MSC,NPTST)
      REAL :: URSELL(MCGRD)
!
!  6. Local variables
!
!     AUX1  :     auxiliary real
!     AUX2  :     auxiliary real
!     BIPH  :     parameterized biphase of the spectrum
!     C0    :     phase velocity at central bin
!     CM    :     phase velocity at interacting bin
!     DEP   :     water depth
!     DEP_2 :     water depth to power 2
!     DEP_3 :     water depth to power 3
!     E     :     energy density as function of frequency
!     E0    :     energy density at central bin
!     EM    :     energy density at interacting bin
!     FT    :     auxiliary real indicating multiplication factor
!                 for triad contribution
!     I1    :     auxiliary integer
!     I2    :     auxiliary integer
!     ID    :     counter
!     IDDUM :     loop counter in direction space
!     IENT  :     number of entries
!     II    :     loop counter
!     IS    :     loop counter in frequency space
!     ISM   :     negative range for IS
!     ISM1  :     negative range for IS
!     ISMAX :     maximum of the counter in frequency space for
!                 which the triad interactions are calculated (cut-off)
!     ISP   :     positive range for IS
!     ISP1  :     positive range for IS
!     RINT  :     interaction coefficient
!     SA    :     interaction contribution of triad
!     SIGPI :     frequency times 2pi
!     SINBPH:     absolute sine of biphase
!     STRI  :     total triad contribution
!     WISM  :     interpolation weight factor corresponding to lower harmonic
!     WISM1 :     interpolation weight factor corresponding to lower harmonic
!     WISP  :     interpolation weight factor corresponding to higher harmonic
!     WISP1 :     interpolation weight factor corresponding to higher harmonic
!     W0    :     radian frequency of central bin
!     WM    :     radian frequency of interacting bin
!     WN0   :     wave number at central bin
!     WNM   :     wave number at interacting bin
!     XIS   :     rate between two succeeding frequency counters
!     XISLN :     log of XIS
!
      INTEGER I1, I2, ID, IDDUM, IENT, II, IS, ISM, ISM1, ISMAX,
     &        ISP, ISP1
      REAL    AUX1, AUX2, BIPH, C0, CM, DEP, DEP_2, DEP_3, E0, EM,
     &        FT, RINT, SIGPI, SINBPH, STRI, WISM, WISM1, WISP, WISP1,
     &        W0, WM, WN0, WNM, XIS, XISLN
      REAL, ALLOCATABLE :: E(:), SA(:,:)
!
!  9. Subroutines calling
!
!     SOURCE
!
! 12. Structure
!
!     Determine resonance condition and the maximum discrete freq.
!     for which the interactions are calculated.
!
!     If Ursell number larger than prescribed value compute interactions
!        Calculate biphase
!        Do for each direction
!           Convert action density to energy density
!           Do for all frequencies
!             Calculate interaction coefficient and interaction factor
!             Compute interactions and store results in matrix
!
! 13. Source text
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'SWLTA')

      DEP   = DEP2(KCGRD(1))
      DEP_2 = DEP**2
      DEP_3 = DEP**3
!
!     --- compute some indices in sigma space
!
      I2     = INT (FLOAT(MSC) / 2.)
      I1     = I2 - 1
      XIS    = SPCSIG(I2) / SPCSIG(I1)
      XISLN  = LOG( XIS )

      ISP    = INT( LOG(2.) / XISLN )
      ISP1   = ISP + 1
      WISP   = (2. - XIS**ISP) / (XIS**ISP1 - XIS**ISP)
      WISP1  = 1. - WISP

      ISM    = INT( LOG(0.5) / XISLN )
      ISM1   = ISM - 1
      WISM   = (XIS**ISM -0.5) / (XIS**ISM - XIS**ISM1)
      WISM1  = 1. - WISM

      ALLOCATE (E (1:MSC))
      ALLOCATE (SA(1:MDC,1:MSC+ISP1))
      E  = 0.
      SA = 0.
!
!     --- compute maximum frequency for which interactions are calculated
!
      ISMAX = 1
      DO IS = 1, MSC
       IF ( SPCSIG(IS) .LT. ( PTRIAD(2) * SMEBRK) ) THEN
          ISMAX = IS
        ENDIF
      ENDDO
      ISMAX = MAX ( ISMAX , ISP1 )
!
!     --- compute 3 wave-wave interactions
!
      IF ( URSELL(KCGRD(1)).GE.PTRIAD(5) ) THEN
!
!       --- calculate biphase
!
        BIPH   = (0.5*PI)*(TANH(PTRIAD(4)/URSELL(KCGRD(1)))-1.)
        SINBPH = ABS( SIN(BIPH) )
!
        DO II = IDDLOW, IDDTOP
           ID = MOD ( II - 1 + MDC , MDC ) + 1
!
!          --- initialize array with E(f) for the direction considered
!
           DO IS = 1, MSC
              E(IS) = AC2(ID,IS,KCGRD(1)) * 2. * PI * SPCSIG(IS)
           END DO
!
           DO IS = 1, ISMAX

              E0  = E(IS)
              W0  = SPCSIG(IS)
              WN0 = KWAVE(IS,1)
              C0  = W0 / WN0

              IF ( IS.GT.-ISM1 ) THEN
                 EM  = WISM * E(IS+ISM1)       + WISM1 * E(IS+ISM)
                 WM  = WISM * SPCSIG(IS+ISM1)  + WISM1 * SPCSIG(IS+ISM)
                 WNM = WISM * KWAVE(IS+ISM1,1) + WISM1 * KWAVE(IS+ISM,1)
                 CM  = WM / WNM
              ELSE
                 EM  = 0.
                 WM  = 0.
                 WNM = 0.
                 CM  = 0.
              END IF

              AUX1 = WNM**2 * ( GRAV * DEP + 2.*CM**2 )
              AUX2 = WN0 * DEP * ( GRAV * DEP +
     &                             (2./15.) * GRAV * DEP_3 * WN0**2 -
     &                             (2./ 5.) * W0**2 * DEP_2 )
              RINT = AUX1 / AUX2
              FT = PTRIAD(1) * C0 * CGO(IS,1) * RINT**2 * SINBPH

              SA(ID,IS) = MAX(0., FT * ( EM * EM - 2. * EM * E0 ))

           END DO
        END DO
!
!        ---  put source term together
!
        DO IS = 1, ISSTOP
           SIGPI = SPCSIG(IS) * 2. * PI
           DO IDDUM = IDCMIN(IS), IDCMAX(IS)
              ID = MOD ( IDDUM - 1 + MDC , MDC ) + 1
!
              STRI = SA(ID,IS) - 2.*(WISP  * SA(ID,IS+ISP1) +
     &                               WISP1 * SA(ID,IS+ISP ))
!
!             --- store results in rhs and main diagonal according
!                 to Patankar-rules
!
              IF(TESTFL) PLTRI(ID,IS,IPTST) = STRI / SIGPI
              IF (STRI.GT.0.) THEN
                 IMATRA(ID,IS) = IMATRA(ID,IS) + STRI / SIGPI
              ELSE
                 IMATDA(ID,IS) = IMATDA(ID,IS) - STRI /
     &                           MAX(1.E-18,AC2(ID,IS,KCGRD(1))*SIGPI)
              END IF
           END DO
        END DO

      END IF
!
!     --- test output
!
      IF ( ITEST .GE. 5 .AND. TESTFL ) THEN
         WRITE(PRINTF,2000) KCGRD(1), ISMAX
 2000    FORMAT (' SWLTA: KCGRD ISMAX  :',2I4)
         WRITE(PRINTF,2001) GRAV, DEP, DEP_2, DEP_3
 2001    FORMAT (' SWLTA: G DEP DEP2 DEP3   :',4E12.4)
         WRITE(PRINTF,2002) PTRIAD(1), PTRIAD(2), URSELL(KCGRD(1))
 2002    FORMAT (' SWLTA: P(1) P(2) P4) URSELL  :',4E12.4)
         WRITE(PRINTF,2003) SMEBRK, HS, BIPH, ABS(SIN(BIPH))
 2003    FORMAT (' SWLTA: SMEBRK HS B |SIN(B)|:',4E12.4)
      END IF

      DEALLOCATE(E,SA)

      RETURN
      END SUBROUTINE SWLTA
!
!********************************************************************
!
      SUBROUTINE STRICL (ACLOC   ,DEPLOC  ,SPCSIG  ,KWAVE   ,
     &                   IDDLOW  ,IDDTOP  ,ANYBIN  ,IMATDA  ,IMATRA,
     &                   CGO     ,KMESPC  ,ETOT    ,SMEBRK          )   40.45
!
!********************************************************************
!
      USE OCPCOMM4
      USE SWCOMM3
      USE SWCOMM4

      IMPLICIT NONE
!
!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering                              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmers: The SWAN team                                |
!   --|-----------------------------------------------------------|--
!
!
!     SWAN (Simulating WAves Nearshore); a third generation wave model
!     Copyright (C) 2008  Delft University of Technology
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!     GNU General Public License for more details.
!
!     A copy of the GNU General Public License is available at
!     http://www.gnu.org/copyleft/gpl.html#SEC3
!     or by writing to the Free Software Foundation, Inc.,
!     59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!
!  0. Authors

!     40.45: Nico Booij
!     40.96: Matthijs B�nit

!  1. Updates

!     40.45, July 04: NEW SUBROUTINE

!  2. Purpose

!     In this subroutine the triad-wave interactions are calculated
!     with the empiric distributed colinear approximation.

!  3. Method

!     Transfer of action between two components under influence of a
!     third one is formulated as follows:

!                  P+1          P+1
!     M * N  (Sigma   N  - Sigma   N  )
!          3       1   1        2   2

!     in which:

!     M      is a dimensional coefficient
!     N      is action density
!     Sigma  spectral (angular) frequency

!     Dimensions:
!     ACL#  : m2s2
!     KW#   : 1/m
!     SIG#  : 1/s
!     SDIA# : 1/s

!  4. Argument variables

      REAL, INTENT(IN)      :: ACLOC(1:MDC,1:MSC)                       ! local action density spectrum
      REAL, INTENT(IN)      :: DEPLOC                                   ! Depth at gridpoint ix,iy (obtained from SWANCOM1)
      REAL, INTENT(IN)      :: SPCSIG(1:MSC)                            ! Relative frequencies in computational domain in sigma-space
      REAL                  :: KWAVE(1:MSC)                             ! Wave number in stencil points
      INTEGER, INTENT(IN)   :: IDDLOW                                   ! Minimum counter in directional space
      INTEGER, INTENT(IN)   :: IDDTOP                                   ! Maximum counter in directional space
      LOGICAL, INTENT(IN)   :: ANYBIN(1:MDC,1:MSC)                      ! if True this bin is going to be updated using the matrix
      REAL                  :: IMATDA(1:MDC,1:MSC)                      ! IMATDA: Diagonal of matrix
      REAL                  :: IMATRA(1:MDC,1:MSC)                      ! IMATRA: Right hand vector of matrix
      REAL, INTENT(IN)      :: CGO(1:MSC)                               ! Group velocities in stencil points
      REAL, INTENT(IN)      :: KMESPC                                   ! Mean wave number of the spectrum
      REAL, INTENT(IN)      :: ETOT                                     !
      REAL, INTENT(IN)      :: SMEBRK                                   !

!     Values from common

!     MDC       : Size of array in theta-direction
!     MSC       : Size of array in sigma-direction
!     ITRIAD    : indicates type of triad formulation
!     PI        : Circular constant Pi
!     MTRIAD    : Size of array containing triad-coefficients
!     PTRIAD    : Tunable coefficients for nonlinear triad sourceterms
!     PTRIAD(1) : Interaction coefficient (labda)
!     PTRIAD(2) : Power of the tail of the spectrum (p)
!     PTRIAD(4) : xxx (delta)

!  6. Local variables

      INTEGER, SAVE     :: IENT=0                                       ! Number of entries into this subroutine
      INTEGER           :: ID                                           ! Grid counter in spectral space (direction)
      INTEGER           :: II                                           ! Counter
      INTEGER           :: IS1, IS2, IS3                                ! Grid counter for spectral frequency
      REAL              :: SIG1, SIG2, SIG3                             ! frequencies of 3 components
      REAL              :: ACL1, ACL2, ACL3                             ! action densities of 2 components
      REAL              :: KW1, KW2, KW3                                ! wave numbers of 3 components
      REAL              :: CG1, CG2                                     ! wave group velocities at interacting frequencies
      REAL              :: RS3, SS                                      ! aux. var. for determining SIG3
      REAL              :: DSIG                                         ! frequency increment
      REAL              :: SIGMEAN                                      ! mean freq of the 3 components
      REAL              :: KMEAN                                        ! mean wave number
      REAL              :: SDIA1, SDIA2                                 ! source term to diagonal
      REAL              :: SRHS
      REAL              :: DISPC                                        ! dispersion coefficient
      REAL              :: BETA
      REAL              :: SINABS
      REAL              :: BIPH
      REAL              :: URSLOC                                       ! aux. coefficients


!  7. SUBROUTINES USED

!      INTEGER :: CNTSIG

!  8. SUBROUTINES CALLING

!     SOURCE

!  9. ERROR MESSAGES

!     ---

! 10. REMARKS

!     ---

! 11. STRUCTURE

!     -----------------------------------------------------------------
!     For all active directions do
!         For all first components do
!             determine frequency and wave number of component
!             For all second components do
!                 determine frequency and wave number of components
!                 determine frequency of third resonating component
!                 If this frequency is within spectral range
!                 Then determine source terms
!                      determine contributions to matrix
!     -----------------------------------------------------------------

! 13. Source text

      IF (LTRACE) CALL STRACE (IENT,'STRICL')

      URSLOC = (GRAV*2.*SQRT(ETOT))/(SQRT(2.)*SMEBRK**2*DEPLOC**2)
      BIPH   = 0.5 * PI * (TANH(PTRIAD(4)/URSLOC) - 1.)
      SINABS = ABS(SIN(BIPH))
      BETA   = PTRIAD(1) / DEPLOC / DEPLOC * SINABS *
     &         KMESPC**(1.-PTRIAD(2))

!     first loop over all directions
      DO II = IDDLOW, IDDTOP
        ID = MOD (II - 1 + MDC, MDC) + 1
!       first loop over all frequencies
        DO IS1 = 1, MSC
          SIG1 = SPCSIG(IS1)
!         second loop over all higher frequencies
          DO IS2 = IS1+1, MSC
            SIG2 = SPCSIG(IS2)
!           determine properties of 3rd component
            SIG3 = SIG2 - SIG1
            IF (SIG3.GT.SPCSIG(1)) THEN

              SS  = ALOG(SIG3/SPCSIG(1)) / FRINTF
              IS3 = INT(SS) + 1
              RS3 = SS - REAL(IS3 - 1)

              ACL1 = ACLOC(ID,IS1)
              ACL2 = ACLOC(ID,IS2)
              ACL3 = RS3 * ACLOC(ID,IS3+1) + (1. - RS3) * ACLOC(ID,IS3)
              KW1  = KWAVE(IS1)
              KW2  = KWAVE(IS2)
              KW3  = RS3 * KWAVE(IS3+1) + (1. - RS3) * KWAVE(IS3)
              CG1  = CGO(IS1)
              CG2  = CGO(IS2)
!             determine properties of triad
              KMEAN = (KW1 + KW2 + KW3)/3.
              DISPC = TANH(KMEAN * DEPLOC) / (KMEAN * DEPLOC)

              SDIA1 = BETA * ACL3 * SIG1 * CG1 * KW1**(PTRIAD(2))
     &                   * DISPC
              SDIA2 = BETA * ACL3 * SIG2 * CG2 * KW2**(PTRIAD(2))
     &                   * DISPC

              IF (ANYBIN(ID,IS1)) THEN
                DSIG = FRINTF * SIG2
                SRHS =  SDIA2 * ACL2
                IMATDA(ID,IS1) = IMATDA(ID,IS1) + SDIA1 * DSIG
                IMATRA(ID,IS1) = IMATRA(ID,IS1) +  SRHS * DSIG
              ENDIF
              IF (ANYBIN(ID,IS2)) THEN
                DSIG = FRINTF * SIG1
                SRHS =  SDIA1 * ACL1
                IMATDA(ID,IS2) = IMATDA(ID,IS2) + SDIA2 * DSIG
                IMATRA(ID,IS2) = IMATRA(ID,IS2) +  SRHS * DSIG
              ENDIF
            ENDIF
          ENDDO
        ENDDO
      ENDDO
      RETURN
      END subroutine STRICL
