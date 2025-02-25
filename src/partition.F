      module partition
        use sizes

        implicit none

        type domain
          integer :: NP, NE, MNP, MNE
          integer, allocatable, dimension(:,:) :: NM
          integer, allocatable, dimension(:) :: WDFLG, DOFS

          real(sz), allocatable, dimension(:,:,:) :: ze,qx,qy
          real(sz), allocatable, dimension(:,:,:) :: rhs_ze,rhs_qx,rhs_qy
          real(sz), allocatable, dimension(:,:,:,:) :: xfac, yfac
          real(sz), allocatable, dimension(:,:,:) :: HB
          real(sz), allocatable, dimension(:) :: FRIC_EL, manningsN
        end type

        integer :: numDomains
        type(domain), allocatable :: domains(:)

        contains

          subroutine initDomains()
          end subroutine


      end module

