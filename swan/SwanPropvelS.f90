subroutine SwanPropvelS ( cad   , cas   , ux2   , uy2   , &
                          dep1  , dep2  , cax   , cay   , &
                          kwave , cgo   , spcsig, idcmin, &
                          idcmax, ecos  , esin  , coscos, &
                          sincos, sinsin, rdx   , rdy   , &
                          jc    )
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering and Geosciences              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmer: Marcel Zijlema                                |
!   --|-----------------------------------------------------------|--
!
!
!     SWAN (Simulating WAves Nearshore); a third generation wave model
!     Copyright (C) 2009  Delft University of Technology
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
!   Authors
!
!   40.80: Marcel Zijlema
!   41.02: Marcel Zijlema
!   41.06: Gerbrant van Vledder
!
!   Updates
!
!   40.80,     July 2007: New subroutine
!   41.02, February 2009: adaption of velocities in case of diffraction
!   41.06,    March 2009: limitation of velocity in theta-direction
!
!   Purpose
!
!   computes wave transport velocities of energy in spectral space
!
!   Modules used
!
    use ocpcomm4
    use swcomm2
    use swcomm3
    use swcomm4
    use m_diffr
    use SwanGriddata
    use SwanGridobjects
    use SwanCompdata
!
    implicit none
!
!   Argument variables
!
    integer, intent(in)                        :: jc     ! counter corresponding to current cell
    !
    integer, dimension(MSC), intent(in)        :: idcmax ! maximum frequency-dependent counter in directional space
    integer, dimension(MSC), intent(in)        :: idcmin ! minimum frequency-dependent counter in directional space
    !
    real, dimension(MDC,MSC), intent(out)      :: cad    ! wave transport velocity in theta-direction
    real, dimension(MDC,MSC), intent(out)      :: cas    ! wave transport velocity in sigma-direction
    real, dimension(MDC,MSC,ICMAX), intent(in) :: cax    ! wave transport velocity in x-direction
    real, dimension(MDC,MSC,ICMAX), intent(in) :: cay    ! wave transport velocity in y-direction
    real, dimension(MSC,ICMAX), intent(in)     :: cgo    ! group velocity
    real, dimension(MDC), intent(in)           :: coscos ! help array containing cosine to power 2 of spectral directions
    real, dimension(nverts), intent(in)        :: dep1   ! water depth at previous time level
    real, dimension(nverts), intent(in)        :: dep2   ! water depth at current time level
    real, dimension(MDC), intent(in)           :: ecos   ! help array containing cosine of spectral directions
    real, dimension(MDC), intent(in)           :: esin   ! help array containing sine of spectral directions
    real, dimension(MSC,ICMAX), intent(in)     :: kwave  ! wave number
    real, dimension(MDC), intent(in)           :: sincos ! help array containing sine * cosine of spectral directions
    real, dimension(MDC), intent(in)           :: sinsin ! help array containing sine to power 2 of spectral directions
    real, dimension(nverts), intent(in)        :: ux2    ! ambient velocity in x-direction at current time level
    real, dimension(nverts), intent(in)        :: uy2    ! ambient velocity in y-direction at current time level
    real, dimension(2), intent(in)             :: rdx    ! first component of contravariant base vector rdx(b) = a^(b)_1
    real, dimension(2), intent(in)             :: rdy    ! second component of contravariant base vector rdy(b) = a^(b)_2
    real, dimension(MSC), intent(in)           :: spcsig ! relative frequency bins
!
!   Local variables
!
    integer                               :: icell    ! cell index
    integer                               :: id       ! loop counter over direction bins
    integer                               :: iddum    ! counter in directional space for considered sweep
    integer, save                         :: ient = 0 ! number of entries in this subroutine
    integer                               :: is       ! loop counter over frequency bins
    integer                               :: iv1      ! first index in computational stencil
    integer                               :: iv2      ! second index in computational stencil
    integer                               :: iv3      ! third index in computational stencil
    integer                               :: jm       ! counter corresponding to cell below considered sweep
    integer                               :: jp       ! counter corresponding to cell above considered sweep
    integer                               :: k        ! loop counter
    integer, dimension(3)                 :: v        ! vertices in considered cell
    !
    real, dimension(28)                   :: cd       ! coefficients for computing cad
    real, dimension(14)                   :: cs       ! coefficients for computing cas
    real, dimension(3)                    :: dloc     ! local depth at vertices
    real                                  :: dpdx     ! derivative of dep2 to x in considered sweep
    real                                  :: dpdy     ! derivative of dep2 to y in considered sweep
    real                                  :: dpdxm    ! derivative of dep2 to x in sweep below considered sweep
    real                                  :: dpdym    ! derivative of dep2 to y in sweep below considered sweep
    real                                  :: dpdxp    ! derivative of dep2 to x in sweep above considered sweep
    real                                  :: dpdyp    ! derivative of dep2 to y in sweep above considered sweep
    real                                  :: duxdx    ! derivative of ux2 to x in considered sweep
    real                                  :: duxdy    ! derivative of ux2 to y in considered sweep
    real                                  :: duydx    ! derivative of uy2 to x in considered sweep
    real                                  :: duydy    ! derivative of uy2 to y in considered sweep
    real                                  :: duxdxm   ! derivative of ux2 to x in sweep below considered sweep
    real                                  :: duxdym   ! derivative of ux2 to y in sweep below considered sweep
    real                                  :: duydxm   ! derivative of uy2 to x in sweep below considered sweep
    real                                  :: duydym   ! derivative of uy2 to y in sweep below considered sweep
    real                                  :: duxdxp   ! derivative of ux2 to x in sweep above considered sweep
    real                                  :: duxdyp   ! derivative of ux2 to y in sweep above considered sweep
    real                                  :: duydxp   ! derivative of uy2 to x in sweep above considered sweep
    real                                  :: duydyp   ! derivative of uy2 to y in sweep above considered sweep
    real                                  :: fac      ! a factor
    real                                  :: frlim    ! frequency range in which limit on velocity in theta-direction is applied
    real                                  :: kd       ! help variable, wave number times water depth
    real                                  :: pp       ! power of the frequency dependent limiter on refraction
    real, dimension(2)                    :: rdxl     ! first component of local contravariant base vector rdx(b) = a^(b)_1
    real, dimension(2)                    :: rdyl     ! second component of local contravariant base vector rdy(b) = a^(b)_2
    !
    type(celltype), dimension(:), pointer :: cell     ! datastructure for cells with their attributes
    type(verttype), dimension(:), pointer :: vert     ! datastructure for vertices with their attributes
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'SwanPropvelS')
    !
    ! point to vertex and cell objects
    !
    vert => gridobject%vert_grid
    cell => gridobject%cell_grid
    !
    ! set velocities to zero
    !
    cad = 0.
    cas = 0.
    !
    ! if no frequency shift and no refraction, return
    !
    if ( (ITFRE == 0 .or. (.not.DYNDEP .and. ICUR == 0)) .and. IREFR == 0 ) return
    !
    ! initialize coefficients
    !
    cd = 0.
    cs = 0.
    !
    iv1 = vs(1)
    iv2 = vs(2)
    iv3 = vs(3)
    !
    dloc(1) = dep2(iv1)
    dloc(2) = dep2(iv2)
    dloc(3) = dep2(iv3)
    !
    ! if at least one vertex is dry, return
    !
    if ( dloc(1) <= DEPMIN .or. dloc(2) <= DEPMIN .or. dloc(3) <= DEPMIN ) return
    !
    ! compute the derivatives of the depth for the considered sweep
    !
    if ( IREFR /= 0 .or. ICUR /= 0 ) then
       !
       if ( IREFR == -1 ) then
          !
          ! limitation of depths in upwave vertices
          !
          dloc(2) = min( dep2(iv2), PNUMS(17)*dloc(1) )
          dloc(3) = min( dep2(iv3), PNUMS(17)*dloc(1) )
          !
       endif
       !
       dpdx = rdx(1) * (dloc(1)-dloc(2)) + rdx(2) * (dloc(1)-dloc(3))
       dpdy = rdy(1) * (dloc(1)-dloc(2)) + rdy(2) * (dloc(1)-dloc(3))
       !
    endif
    !
    ! compute the derivatives of the ambient current for the considered sweep
    !
    if ( ICUR /= 0 ) then
       !
       duxdx = rdx(1) * (ux2(iv1) - ux2(iv2)) + rdx(2) * (ux2(iv1) - ux2(iv3))
       duxdy = rdy(1) * (ux2(iv1) - ux2(iv2)) + rdy(2) * (ux2(iv1) - ux2(iv3))
       duydx = rdx(1) * (uy2(iv1) - uy2(iv2)) + rdx(2) * (uy2(iv1) - uy2(iv3))
       duydy = rdy(1) * (uy2(iv1) - uy2(iv2)) + rdy(2) * (uy2(iv1) - uy2(iv3))
       !
    endif
    !
    ! compute extra derivatives of depth and current belonging to neighbouring sweeps meant for refraction
    !
    if ( IREFR /= 0 ) then
       !
       dpdxm = dpdx
       dpdym = dpdy
       dpdxp = dpdx
       dpdyp = dpdy
       !
       if ( ICUR /= 0 ) then
          !
          duxdxm = duxdx
          duxdym = duxdy
          duydxm = duydx
          duydym = duydy
          duxdxp = duxdx
          duxdyp = duxdy
          duydxp = duydx
          duydyp = duydy
          !
       endif
       !
       jm=vert(iv1)%noc-mod(vert(iv1)%noc+1-jc,vert(iv1)%noc)        ! cell below considered sweep
       !
       icell = vert(iv1)%cell(jm)%atti(CELLID)
       !
       v(1) = cell(icell)%atti(CELLV1)
       v(2) = cell(icell)%atti(CELLV2)
       v(3) = cell(icell)%atti(CELLV3)
       !
       ! compute local contravariant base vectors at present vertex in cell below considered sweep
       !
       do k = 1, 3
          if ( v(k) == iv1 ) then
             !
             rdxl(1) = cell(icell)%geom(k)%rdx1
             rdxl(2) = cell(icell)%geom(k)%rdx2
             rdyl(1) = cell(icell)%geom(k)%rdy1
             rdyl(2) = cell(icell)%geom(k)%rdy2
             !
             iv2 = v(mod(k  ,3)+1)
             iv3 = v(mod(k+1,3)+1)
             !
             exit
          endif
       enddo
       !
       if ( KSPHER > 0 ) then
          do k = 1, 2
             rdxl(k) = rdxl(k) / (COSLAT(1) * LENDEG)
             rdyl(k) = rdyl(k) / LENDEG
          enddo
       endif
       !
       dloc(2) = dep2(iv2)
       dloc(3) = dep2(iv3)
       !
       ! if at least one vertex is dry, skip and go to other cell
       !
       if ( dloc(2) <= DEPMIN .or. dloc(3) <= DEPMIN ) goto 10
       !
       ! compute the derivatives of the depth
       !
       if ( IREFR == -1 ) then
          !
          ! limitation of depths in upwave vertices
          !
          dloc(2) = min( dep2(iv2), PNUMS(17)*dloc(1) )
          dloc(3) = min( dep2(iv3), PNUMS(17)*dloc(1) )
          !
       endif
       !
       dpdxm = rdxl(1) * (dloc(1)-dloc(2)) + rdxl(2) * (dloc(1)-dloc(3))
       dpdym = rdyl(1) * (dloc(1)-dloc(2)) + rdyl(2) * (dloc(1)-dloc(3))
       !
       ! compute the derivatives of the ambient current
       !
       if ( ICUR /= 0 ) then
          !
          duxdxm = rdxl(1) * (ux2(iv1) - ux2(iv2)) + rdxl(2) * (ux2(iv1) - ux2(iv3))
          duxdym = rdyl(1) * (ux2(iv1) - ux2(iv2)) + rdyl(2) * (ux2(iv1) - ux2(iv3))
          duydxm = rdxl(1) * (uy2(iv1) - uy2(iv2)) + rdxl(2) * (uy2(iv1) - uy2(iv3))
          duydym = rdyl(1) * (uy2(iv1) - uy2(iv2)) + rdyl(2) * (uy2(iv1) - uy2(iv3))
          !
       endif
       !
  10   continue
       !
       jp = mod(jc,vert(iv1)%noc)+1                                  ! cell above considered sweep
       !
       icell = vert(iv1)%cell(jp)%atti(CELLID)
       !
       v(1) = cell(icell)%atti(CELLV1)
       v(2) = cell(icell)%atti(CELLV2)
       v(3) = cell(icell)%atti(CELLV3)
       !
       ! compute local contravariant base vectors at present vertex in cell above considered sweep
       !
       do k = 1, 3
          if ( v(k) == iv1 ) then
             !
             rdxl(1) = cell(icell)%geom(k)%rdx1
             rdxl(2) = cell(icell)%geom(k)%rdx2
             rdyl(1) = cell(icell)%geom(k)%rdy1
             rdyl(2) = cell(icell)%geom(k)%rdy2
             !
             iv2 = v(mod(k  ,3)+1)
             iv3 = v(mod(k+1,3)+1)
             !
             exit
          endif
       enddo
       !
       if ( KSPHER > 0 ) then
          do k = 1, 2
             rdxl(k) = rdxl(k) / (COSLAT(1) * LENDEG)
             rdyl(k) = rdyl(k) / LENDEG
          enddo
       endif
       !
       dloc(2) = dep2(iv2)
       dloc(3) = dep2(iv3)
       !
       ! if at least one vertex is dry, skip
       !
       if ( dloc(2) <= DEPMIN .or. dloc(3) <= DEPMIN ) goto 20
       !
       ! compute the derivatives of the depth
       !
       if ( IREFR == -1 ) then
          !
          ! limitation of depths in upwave vertices
          !
          dloc(2) = min( dep2(iv2), PNUMS(17)*dloc(1) )
          dloc(3) = min( dep2(iv3), PNUMS(17)*dloc(1) )
          !
       endif
       !
       dpdxp = rdxl(1) * (dloc(1)-dloc(2)) + rdxl(2) * (dloc(1)-dloc(3))
       dpdyp = rdyl(1) * (dloc(1)-dloc(2)) + rdyl(2) * (dloc(1)-dloc(3))
       !
       ! compute the derivatives of the ambient current
       !
       if ( ICUR /= 0 ) then
          !
          duxdxm = rdxl(1) * (ux2(iv1) - ux2(iv2)) + rdxl(2) * (ux2(iv1) - ux2(iv3))
          duxdym = rdyl(1) * (ux2(iv1) - ux2(iv2)) + rdyl(2) * (ux2(iv1) - ux2(iv3))
          duydxm = rdxl(1) * (uy2(iv1) - uy2(iv2)) + rdxl(2) * (uy2(iv1) - uy2(iv3))
          duydym = rdyl(1) * (uy2(iv1) - uy2(iv2)) + rdyl(2) * (uy2(iv1) - uy2(iv3))
          !
       endif
       !
  20   continue
       !
    endif
    !
    ! compute wave transport velocity in sigma-direction
    !
    if ( ITFRE /= 0 .and. (DYNDEP .or. ICUR /= 0) ) then
       !
       ! compute time derivative of water depth
       !
       if ( DYNDEP ) cs(2) = ( dep2(iv1) - dep1(iv1) ) * RDTIM
       !
       ! compute coefficients depending on ambient currents
       !
       if ( ICUR /= 0 ) then
          !
          cs(3) = ux2(iv1) * dpdx
          cs(4) = uy2(iv1) * dpdy
          cs(6) = duxdx
          cs(7) = duxdy
          cs(8) = duydx
          cs(9) = duydy
          !
       endif
       !
       do is = 1, MSC
          !
          ! compute frequency-dependent coefficients
          !
          kd = min(30.,kwave(is,1) * dep2(iv1))
          !
          cs(1) =  kwave(is,1) * spcsig(is) / sinh (2.* kd)
          cs(5) = -cgo(is,1) * kwave(is,1)
          if ( IDIFFR /= 0 ) cs(5) = cs(5)*DIFPARAM(iv1)
          !
          cs(10) = cs(1) * cs(2)
          cs(11) = cs(1) * (cs(3)+cs(4))
          cs(12) = cs(5) * cs(6)
          cs(13) = cs(5) * (cs(7)+cs(8))
          cs(14) = cs(5) * cs(9)
          !
          do iddum = idcmin(is), idcmax(is)
             id = mod ( iddum - 1 + MDC , MDC ) + 1
             !
             cas(id,is) = cs(10)
             if ( ICUR /= 0 ) cas(id,is) = cas(id,is) + cs(11) + coscos(id)*cs(12) + sincos(id)*cs(13) + sinsin(id)*cs(14)
!            !
          enddo
       enddo
       !
    endif
    !
    ! compute wave transport velocity in theta-direction
    !
    if ( IREFR /= 0 ) then
       !
       ! compute coefficients depending on depth
       !
       cd( 2) = dpdx
       cd( 3) = dpdy
       !
       cd(11) = dpdxm
       cd(12) = dpdym
       !
       cd(20) = dpdxp
       cd(21) = dpdyp
       !
       ! compute coefficients depending on ambient currents
       !
       if ( ICUR /= 0 ) then
          !
          cd( 4) = duxdx
          cd( 5) = duydy
          cd( 6) = duydx
          cd( 7) = duxdy
          !
          cd(13) = duxdxm
          cd(14) = duydym
          cd(15) = duydxm
          cd(16) = duxdym
          !
          cd(22) = duxdxp
          cd(23) = duydyp
          cd(24) = duydxp
          cd(25) = duxdyp
          !
       endif
       !
       do is = 1, MSC
          !
          ! compute frequency-dependent coefficients
          !
          kd = min(30.,kwave(is,1) * dep2(iv1))
          !
          cd(1) = spcsig(is) / sinh (2.* kd)
          !
          cd( 8) = cd(1) * cd(2)
          cd( 9) = cd(1) * cd(3)
          cd(10) = cd(4) - cd(5)
          !
          do iddum = idcmin(is), idcmax(is)
             id = mod ( iddum - 1 + MDC , MDC ) + 1
             !
             cad(id,is) = esin(id)*cd(8) - ecos(id)*cd(9)
             if ( IDIFFR /= 0 ) cad(id,is) = cad(id,is)*DIFPARAM(iv1) - DIFPARDX(iv1)*cgo(is,1)*esin(id) + DIFPARDY(iv1)*cgo(is,1)*ecos(id)
             if ( ICUR   /= 0 ) cad(id,is) = cad(id,is) + sincos(id)*cd(10) + sinsin(id)*cd(6) - coscos(id)*cd(7)
             !
          enddo
          !
          cd(17) = cd( 1) * cd(11)
          cd(18) = cd( 1) * cd(12)
          cd(19) = cd(13) - cd(14)
          !
          id = mod ( idcmin(is) - 2 + MDC , MDC ) + 1  ! this direction belongs to sweep below considered sweep
          !
          cad(id,is) = esin(id)*cd(17) - ecos(id)*cd(18)
          if ( IDIFFR /= 0 ) cad(id,is) = cad(id,is)*DIFPARAM(iv1) - DIFPARDX(iv1)*cgo(is,1)*esin(id) + DIFPARDY(iv1)*cgo(is,1)*ecos(id)
          if ( ICUR   /= 0 ) cad(id,is) = cad(id,is) + sincos(id)*cd(19) + sinsin(id)*cd(15) - coscos(id)*cd(16)
          !
          cd(26) = cd( 1) * cd(20)
          cd(27) = cd( 1) * cd(21)
          cd(28) = cd(22) - cd(23)
          !
          id = mod ( idcmax(is) + MDC , MDC ) + 1      ! this direction belongs to sweep above considered sweep
          !
          cad(id,is) = esin(id)*cd(26) - ecos(id)*cd(27)
          if ( IDIFFR /= 0 ) cad(id,is) = cad(id,is)*DIFPARAM(iv1) - DIFPARDX(iv1)*cgo(is,1)*esin(id) + DIFPARDY(iv1)*cgo(is,1)*ecos(id)
          if ( ICUR   /= 0 ) cad(id,is) = cad(id,is) + sincos(id)*cd(28) + sinsin(id)*cd(24) - coscos(id)*cd(25)
          !
       enddo
       !
       ! adapt velocity in case of spherical coordinates
       !
       if ( KSPHER > 0 ) then
          !
          cd(1) = tan(DEGRAD*(vert(iv1)%attr(VERTY) + YOFFS)) / REARTH
          !
          do id = 1, MDC
             cd(2) = cd(1) * ecos(id)
             do is = 1, MSC
                cad(id,is) = cad(id,is) - cd(2)*(cax(id,is,1)*ecos(id) + cay(id,is,1)*esin(id))
             enddo
          enddo
          !
       endif
       !
       ! limit velocity in some frequency range if requested
       !
       if ( IREFR == -2 ) then
          !
          frlim = PI2*PNUMS(26)
          pp    =     PNUMS(27)
          !
          do is = 1, MSC
             !
             fac = min(1.,(spcsig(is)/frlim)**pp)
             !
             do id = 1, MDC
                cad(id,is) = fac*cad(id,is)
             enddo
          enddo
          !
       endif
       !
    endif
    !
end subroutine SwanPropvelS
