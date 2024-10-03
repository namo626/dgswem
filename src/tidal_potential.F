C***********************************************************************
C
C     SUBROUTINE TIDAL_POTENTIAL()
C
C     This subroutine computes the tidal potential forcing terms at each
C     node.  Note that the tidal potential term is defined linearly over
C     an element regardless of the p used in the DG calculation.
C
C     (Taken from timestep.f)
C
C     Comments and general clean-up by Ethan Kubatko (06-02-2005)
C
C***********************************************************************

      module tidal_dg

      use sizes, only : sz
      use global, only : tip2, facet, sfea, slam, fft, saltpha, etrf, pert, amigt,
     &  saltamp, tpk
      use dg, only : dtvd

      contains

      SUBROUTINE TIDAL_POTENTIAL(irk, np, dtdp, timeh, ntif, rampdg)

      IMPLICIT NONE

      integer, value :: irk, np, ntif
      real(sz), value :: dtdp, timeh, rampdg

      INTEGER :: j,i, ncyc, argt, na
      real(sz) :: argsalt, saltmul, tpmul, argtp, ccsfea, s2sfea, timeh_dg
      
      TIMEH_DG = TIMEH - DTDP + DTVD(IRK)*DTDP
      
C.....Loop over the nodes

!$cuf kernel do <<< *, 128 >>>
      DO I = 1,NP
      
C.....Initialize tidal potential terms

        TIP2(I) = 0.D0

C.....Loop over the tidal potential constituents

        DO J = 1,NTIF

          IF (PERT(J).EQ.0.) THEN
            NCYC = 0
          ELSE
            NCYC = INT(TIMEH_DG/PERT(J))
          ENDIF
          ARGT    = AMIGT(J)*(TIMEH_DG - NCYC*PERT(J)) + FACET(J)
          TPMUL   = RAMPDG*ETRF(J)*TPK(J)*FFT(J)
          SALTMUL = RAMPDG*FFT(J)
          NA      = NINT(0.00014/AMIGT(J))
        
C.....Semi-diurnal species

          IF (NA.EQ.1) THEN

            ARGTP    = ARGT + 2.D0*SLAM(I)
            ARGSALT  = ARGT - SALTPHA(J,I)
            CCSFEA   = COS(SFEA(I))
            CCSFEA   = CCSFEA*CCSFEA
            TIP2(I) = TIP2(I) + TPMUL*CCSFEA*COS(ARGTP)
     &                        + SALTMUL*SALTAMP(J,I)*COS(ARGSALT)

          ENDIF
        
C.....Diurnal species
        
          IF (NA.EQ.2) THEN

            ARGTP    = ARGT + SLAM(I)
            ARGSALT  = ARGT - SALTPHA(J,I)
            S2SFEA   = SIN(2.D0*SFEA(I))
            TIP2(I) = TIP2(I) + TPMUL*S2SFEA*COS(ARGTP)
     &                        + SALTMUL*SALTAMP(J,I)*COS(ARGSALT)
     
          ENDIF

        ENDDO

      ENDDO
      
      RETURN
      END SUBROUTINE


      end module
