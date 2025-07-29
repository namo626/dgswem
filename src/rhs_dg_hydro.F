C***********************************************************************
C
C     SUBROUTINE RHS_DG_HYDRO()
C
C     This subroutine computes the area integrals for the DG hydro and
C     adds them into the RHS.
C
C     Written by Ethan Kubatko (06-11-2004)
C
C-----------------------------------------------------------------------
C
C     Feb 23, 2007, sb, Modified for better performance
C     Jan 02, 2007, sb, Modified for LDG
C     Aug xx, 2005, sb, Modified for wetting/drying
C     01-10-2011 - cem - adapted for p_enrichment and multicomponent
C
C***********************************************************************
      SUBROUTINE RHS_DG_HYDRO(IT)

C.....Use appropriate modules

      USE GLOBAL
      USE DG
      USE NodalAttributes, ONLY : TAU, IFLINBF, IFHYBF, HBREAK, FTHETA,
     &     FGAMMA,LoadManningsN,ManningsN,CF

      USE sizes, ONLY: myproc,layers
      USE PRECIPITATION, only : computeElemTRR, PREC3, PREC2

      IMPLICIT NONE

C.....Declare local variables
      ! for rain test
      REAL(SZ) SR1, SR2, SR3, avg_rain
      real(sz) :: ze_in,ze_ex,qx_in,qx_ex,qy_in,qy_ex,
     &  hb_in,hb_ex,sfac_in,sfac_ex, nx, ny
      real(sz) :: f_hat, g_hat, h_hat
      INTEGER IT
      INTEGER L,k,i,ll
      REAL(SZ) DPSIDX1, dpsidx2, dpsidx3, DPSIDY1, dpsidy2, dpsidy3
      REAL(SZ) AREA, IMASS, TKX, TKY, Xpart, Ypart,H_0,C_1
      REAL(SZ) PHI_AREA_KI,MN_IN, MassAction1st,MassAction2nd,fx,fy
      REAL(SZ) LZ_XX, LZ_XY, LZ_YX, LZ_YY, rate, s_mass, s_sed,b_0
      REAL(SZ) DEPTH, F1_NL, FU_NL, FV_NL, FG_NL, FH_NL, FW_NL
      REAL(SZ) HUU, HVV, HUV, GH2,MZ_X(layers),MZ_Y(layers), fgauss, sig
      REAL(SZ) DEPTH_C, FH_NL_C, UX_C, UY_C, UMAG_C, DTDPH,SFACQUAD
      Real(SZ) discharge_modelX_IN,discharge_modelY_IN
      Real(SZ) DH_X,DH_Y,phi_tot,C_0,HZ_X,HZ_Y,TZ_X,TZ_Y

      real(sz), dimension(3,NE) :: ze_area, qx_area, qy_area
      real(sz), dimension(3,NE) :: fxs,fys,gxs,gys,hxs,hys,source_xs,source_ys

      real(sz) source_max, source_r

      source_max = 0.0

      pa = 1
      DTDPH = 1.0/DTDP
C.......Compute ZE, QX, QY, and HB at each area Gauss quadrature point

!$acc parallel loop gang vector collapse(2) async(1) default(present)
      DO L = 1, NE
        DO I = 1,NAGP(pa)
          ze_area(i,l) = ze(1,L,irk)
          qx_area(i,l) = qx(1,L,irk)
          qy_area(i,l) = qy(1,L,irk)
          DO K = 2,DOFS(L)
            ZE_area(I,L) = ZE_area(I,L) + ZE(K,L,IRK)*PHI_AREA(K,I,pa)
            QX_area(I,L) = QX_area(I,L) + QX(K,L,IRK)*PHI_AREA(K,I,pa)
            QY_area(I,L) = QY_area(I,L) + QY(K,L,IRK)*PHI_AREA(K,I,pa)

          ENDDO
        ENDDO
      ENDDO

!$acc parallel loop gang vector collapse(2)  async(1) default(present)
      DO L = 1,MNE
        DO I = 1,3


C.......Retrieve the global node numbers for the element

        N1 = NM(L,1)
        N2 = NM(L,2)
        N3 = NM(L,3)

C.......Compute avaraged values
C.......These will be used later when bottom friction is computed


C.......Compute derivatives of Lagrange basis functions at nodes

        IF ((NWS.NE.0).OR.(NTIP.NE.0)) THEN
          DPSIDX1 = DRPSI(1)*DRDX(L) + DSPSI(1)*DSDX(L)
          DPSIDX2 = DRPSI(2)*DRDX(L) + DSPSI(2)*DSDX(L)
          DPSIDX3 = DRPSI(3)*DRDX(L) + DSPSI(3)*DSDX(L)
          DPSIDY1 = DRPSI(1)*DRDY(L) + DSPSI(1)*DSDY(L)
          DPSIDY2 = DRPSI(2)*DRDY(L) + DSPSI(2)*DSDY(L)
          DPSIDY3 = DRPSI(3)*DRDY(L) + DSPSI(3)*DSDY(L)
        ENDIF
          ze_in = ze_area(I,L)
          qx_in = qx_area(I,L)
          qy_in = qy_area(I,L)

          HB_IN = BATH(I,L,pa)
          DHB_X = DBATHDX(I,L,pa)
          DHB_Y = DBATHDY(I,L,pa)

          DEPTH = ZE_in + HB_IN

          SFACQUAD = SFAC_ELEM(I,L,pa)
C.........Compute continuity fluxes

          F1_NL = NLEQ + LEQ*HB_IN

          FX_IN = QX_IN*F1_NL*SFACQUAD

          FY_IN = QY_IN*F1_NL
C.........Compute momentum flux terms

          FU_NL = NLEQ*QX_IN
          FV_NL = NLEQ*QY_IN
          FG_NL = NLEQG*ZE_IN*WDFLG(L)
          FH_NL = 1.0/(NLEQ*DEPTH + LEQ)
          U_IN  = QX_IN*FH_NL
          V_IN  = QY_IN*FH_NL

          HUU = FU_NL*U_IN
          HVV = FV_NL*V_IN
          HUV = FU_NL*V_IN
          GH2 = FG_NL*(0.5D0*ZE_IN + HB_IN) + FG_L*ZE_IN

C.........Compute x momentum fluxes

          GX_IN = (HUU + GH2 + LZ_XX)*SFACQUAD
          GY_IN = HUV 


C.........Compute y momentum fluxes

          HX_IN = (HUV )*SFACQUAD
          HY_IN = HVV + GH2 


C.........Compute the friction factor
          if (LoadManningsN) then
            fric_el(L)=G*
     $              ((ManningsN(n1)+ManningsN(n2)+ManningsN(n3))/3.)**2
     $              /(DEPTH**(1.0/3.0))
            if (fric_el(L).lt.CF) fric_el(L)=CF
          endif
          TAU = FRIC_EL(L)
C     Modified to compute TAU using elemental averages.
C     This seems necessary to avoid exessive bottom friction
C     at wetting-drying fronts where the total column height is very
C     small. S.B. 9-Feb-2008

            IF (IFLINBF.EQ.0) THEN
               UMAG = SQRT( U_IN*U_IN + V_IN*V_IN )
c     cnd modified 4/23/10 to test friction
               TAU  = TAU*UMAG*FH_NL
               IF (IFHYBF.EQ.1) TAU = TAU*
     &              (1.0  + (HBREAK*FH_NL_C)**FTHETA)**(FGAMMA/FTHETA)
C     It is numerically probable that the bottom friction is large enoght
C     to reverse the direction of currents backward due to a too small column
C     height even though it does not happen in reality. To avoid this, the MIN
C     function bellow is added. It is expected that this MIN function upper-limits
C     TAU so the bottom friction force does not reverse the currents within
C     half a time step.  S.B. 9-Feb-2008
               TAU = MIN(TAU, .9D0*DTDPH)
            ENDIF

C.........Compute the x momentum source/sink terms

            SOURCE_X =

C.........1.) Friction term

     &           - TAU*QX_IN

C.........2.) Bathymetric slope term

     &           + FG_NL*DHB_X*SFACQUAD

C.........3.) Coriolis force

     &           + CORI_EL(L)*QY_IN

C.....Compute the y momentum source/sink terms

            SOURCE_Y =

C.........1.) Friction term

     &           - TAU*QY_IN

C.........2) Bathymetric slope term


     &           + FG_NL*DHB_Y

C.........3.) Coriolis force

     &           - CORI_EL(L)*QX_IN

C.........4.) Wind and pressure forcing (in x and y)

#if 0
            IF (NWS.NE.0) THEN
               FW_NL = 1.0/F1_NL
               SOURCE_X = SOURCE_X + FW_NL*( WSX2(N1)*PSI1(I,pa)
     &              + WSX2(N2)*PSI2(I,pa)  + WSX2(N3)*PSI3(I,pa) )
     &              - G*SFACQUAD*DEPTH
     $              *( PR2(N1)*DPSIDX1
     &              + PR2(N2)*DPSIDX2 + PR2(N3)*DPSIDX3)
               SOURCE_Y = SOURCE_Y + FW_NL*( WSY2(N1)*PSI1(I,pa)
     &              + WSY2(N2)*PSI2(I,pa)  + WSY2(N3)*PSI3(I,pa) )
     &              - G*DEPTH*( PR2(N1)*DPSIDY1
     &              + PR2(N2)*DPSIDY2 + PR2(N3)*DPSIDY3)
            ENDIF


C.........5) Tidal potential forcing (in x and y)
            IF (NTIP.NE.0) THEN

               SOURCE_X = SOURCE_X +
     $              RAMPDG*G*DEPTH*SFACQUAD*( DPSIDX1*TIP2(N1)
     &              + DPSIDX2*TIP2(N2) + DPSIDX3*TIP2(N3) )

               SOURCE_Y = SOURCE_Y + RAMPDG*G*DEPTH*( DPSIDY1*TIP2(N1)
     &              + DPSIDY2*TIP2(N2) + DPSIDY3*TIP2(N3) )
            ENDIF
            fxs(i,l) = fx_in*wdflg(l)
            fys(i,l) = fy_in*wdflg(l)
            gxs(i,l) = gx_in*wdflg(l)
            gys(i,l) = gy_in*wdflg(l)
            hxs(i,l) = hx_in*wdflg(l)
            hys(i,l) = hy_in*wdflg(l)
            source_xs(i,l) = source_x*wdflg(l)
            source_ys(i,l) = source_y*wdflg(l)
#endif
          DO k = 1,3
!$acc atomic
            rhs_ze(L,K,IRK) = rhs_ze(L,K,IRK) +
     &       (XFAC(K,I,L,pa)*FX_in  + YFAC(K,I,L,pa)*FY_in)
!$acc atomic
            rhs_qx(L,K,IRK) = rhs_qx(L,K,IRK) +
     &       (XFAC(K,I,L,pa)*GX_in
     &              + YFAC(K,I,L,pa)*GY_in +
     &         SRFAC(K,I,L,pa)*SOURCE_X)
!$acc atomic
            rhs_qy(L,K,IRK) = rhs_qy(L,K,IRK) +
     &       (XFAC(K,I,L,pa)*HX_IN
     &              + YFAC(K,I,L,pa)*HY_IN +
     &           SRFAC(K,I,L,pa)*SOURCE_Y)

          ENDDO
        ENDDO
      ENDDO

      END SUBROUTINE
