module FORT_DG

   ! Module containing subroutines to read the fort.dg file
   !
   !  - Supports both fixed and keyword formats
   !
   !  - Keyword format is intended to increase flexibility in adding/depreciating features
   !    while maintaining forward compatibility and some degree of backward compatibility
   !    for the fort.dg file.
   !
   !      * The keywords are configured in the FORT_DG_SETUP subroutine.
   !        (See subrotine header for details.)
   !
   !      * keyword fort.dg file format rules are:
   !
   !          1) options are assigned in keyword = value format (e.g. fluxtype = 1)
   !          2) one option per line
   !          3) options can be specified in any order
   !          4) line beginning with ! indicates a comment
   !          5) blank lines are skipped
   !          6) comments may follow an assignment (e.g. fluxtype = 1 !comment)
   !          7) unrecognized keyword assignments are skipped
   !          8) unassigned options that are not required will use
   !             default values specified in FORT_DG_SETUP
   !
   !  - Subroutines contained are:
   !
   !      1) READ_FIXED_FORT_DG
   !
   !         * reads old fixed format fort.dg used in dgswem v11.13/dg-adcirc v22
   !
   !      2) READ_KEYWORD_FORT_DG
   !
   !         * reads keyword format fort.dg described above
   !
   !      3) CHECK_ERRORS
   !
   !         * Handles missing options
   !         * Terminates if required options are missing
   !         * Warns that default values are used for missing optional options and continues
   !
   !      4) FORT_DG_SETUP
   !
   !         * Responsible for configuring fort.dg options
   !         * MODIFICATIONS FOR ADDITION/REMOVAL OF FORT.DG OPTIONS SHOULD BE DONE HERE

   use sizes, only: sz

   type :: key_val
      character(15) :: key ! keyword
      real(SZ), pointer :: rptr ! pointer to real target
      integer, pointer :: iptr ! pointer to integer target
      character(100), pointer :: cptr ! pointer to character target

      integer :: vartype ! target type indicator: 1=integer, 2=real, 3=character

      integer :: required ! required/optional flag

      integer :: flag ! successful read flag
   end type key_val

   integer, parameter :: maxopt = 100 ! maximum allowable fort.dg options
   type(key_val), dimension(maxopt) :: fortdg

   integer :: nopt ! number of valid options in fortdg structure
   integer, dimension(maxopt) :: fortdg_ind ! indicies of valid options in fortdg structure

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine READ_FIXED_FORT_DG()

      use global, only: dgswe, dg_to_cg, sedflag, reaction_rate, sed_equationX, sed_equationY, &
                        rhowat0, vertexslope
      use sizes, only: myproc, layers, dirname
      use dg, only: padapt, pflag, gflag, diorism, pl, ph, px, slimit, plimit, &
                    pflag2con1, pflag2con2, lebesgueP, fluxtype, rk_stage, rk_order, &
                    modal_ic, dghot, dghotspool, slopeflag, slope_weight, porosity, &
                    sevdm, mnes, artdif, kappa, s0, uniform_dif, tune_by_hand, &
                    sl2_m, sl2_nyu, sl3_md, rainfall

      implicit none

      integer :: i
      character(256) :: LINE

      call FORT_DG_SETUP()

      open (25, FILE=DIRNAME//'/'//'fort.dg', POSITION="rewind")

      if (myproc == 0) then
         print *, ""
         print("(A)"), "READING FIXED FORMAT FORT.DG..."
         print *, ""
      end if

      read (25, *) DGSWE
      read (25, *) padapt, pflag
      read (25, *) gflag, diorism
      read (25, *) pl, ph, px
      read (25, *) slimit
      read (25, *) plimit
      read (25, *) pflag2con1, pflag2con2, lebesgueP
      read (25, *) FLUXTYPE
      read (25, *) RK_STAGE, RK_ORDER
      read (25, *) DG_TO_CG
      read (25, *) MODAL_IC
      read (25, *) DGHOT, DGHOTSPOOL
      read (25, "(A256)") LINE
      read (LINE, *) SLOPEFLAG
      if (SLOPEFLAG == 2) then
         read (LINE, *) SLOPEFLAG, SL2_M, SL2_NYU
      end if
      if (SLOPEFLAG == 3) then
         read (LINE, *) SLOPEFLAG, SL2_M, SL2_NYU, SL3_MD
      end if
      if (SLOPEFLAG == 4) then
         read (LINE, *) SLOPEFLAG, slope_weight
         vertexslope = .true.
      end if
      if (SLOPEFLAG == 5) then
         read (LINE, *) SLOPEFLAG
         vertexslope = .true.
      end if
      if (SLOPEFLAG == 6) then
         read (LINE, *) SLOPEFLAG, slope_weight
         vertexslope = .true.
      end if
      if (SLOPEFLAG == 7) then
         read (LINE, *) SLOPEFLAG, slope_weight
         vertexslope = .true.
      end if
      if (SLOPEFLAG == 8) then
         read (LINE, *) SLOPEFLAG, slope_weight
         vertexslope = .true.
      end if
      if (SLOPEFLAG == 9) then
         read (LINE, *) SLOPEFLAG, slope_weight
         vertexslope = .true.
      end if
      if (SLOPEFLAG == 10) then
         read (LINE, *) SLOPEFLAG, slope_weight
         vertexslope = .true.
      end if
      read (25, *) SEDFLAG, porosity, SEVDM, layers
      read (25, *) reaction_rate
      read (25, *) MNES
      read (25, *) artdif, kappa, s0, uniform_dif, tune_by_hand
      read (25, '(a)') sed_equationX
      read (25, '(a)') sed_equationY

      if (FLUXTYPE /= 1 .and. FLUXTYPE /= 2 .and. FLUXTYPE /= 3 .and. FLUXTYPE /= 4) then
         if (myproc == 0) then
            print *, 'SPECIFIED FLUXTYPE (=', FLUXTYPE, ') IS NOT ALLOWED.'
            print *, 'EXECUTION WILL BE TERMINATED.'
         end if

         stop
      end if

      ! print inputs
      if (myproc == 0) then
         do i = 1, maxopt
            if (associated(fortdg(i)%iptr)) then
               print("(A,A,I8)"), fortdg(i)%key, " = ", fortdg(i)%iptr
            end if

            if (associated(fortdg(i)%rptr)) then
               print("(A,A,E21.8)"), fortdg(i)%key, " = ", fortdg(i)%rptr
            end if

            if (associated(fortdg(i)%cptr)) then
               print("(A,A,A)"), fortdg(i)%key, " = ", fortdg(i)%cptr
            end if
         end do
         print *, " "
      end if

      close (25)

      return
   end subroutine READ_FIXED_FORT_DG

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine READ_KEYWORD_FORT_DG()

      use sizes, only: myproc, dirname
      use global, only: nfover

      implicit none

      integer :: i, j, opt
      integer :: read_stat
      integer :: opt_read
      integer :: comment, blank
      integer :: eqind, exind
      integer :: found
      character(100) :: temp, line, empty
      character(15) :: test_opt
      character(100) :: test_val

      ! initialize the fortdg option structure
      call FORT_DG_SETUP()

      opt_read = 0
      comment = 0
      blank = 0

      open (25, FILE=DIRNAME//'/'//'fort.dg', POSITION="rewind")
      if (myproc == 0) then
         print *, ""
         print("(A)"), "READING KEYWORD FORMAT FORT.DG..."
         print *, ""
      end if

      do while (opt_read < nopt)

         read (25, "(A100)", IOSTAT=read_stat) temp
         if (read_stat /= 0) then ! check for end-of-file
            exit
         end if

         line = adjustl(temp)

         if (index(line, "!") == 1) then ! lines beginning with ! are skipped

            comment = comment + 1

         else if (len(trim(line)) == 0) then ! blank lines are skipped

            blank = blank + 1

         else

            ! determine keyword and assignment value
            eqind = index(line, "=")
            exind = index(line, "!")
            test_opt = line(1:eqind - 1)
            if (exind > 0) then ! handle trailing comment
               test_val = adjustl(line(eqind + 1:exind - 1)) ! (only necessary if there is no space between value and the !)
            else
               test_val = adjustl(line(eqind + 1:))
            end if

            ! Look for a match for the keyword
            found = 0
            test: do opt = 1, nopt

               i = fortdg_ind(opt)

               if (test_opt == fortdg(i)%key) then

                  ! Set variables equal to value from fort.dg through pointer using an internal read
                  select case (fortdg(i)%vartype)
                  case (1)
                     read (test_val, *) fortdg(i)%iptr
                     if (myproc == 0) print("(A,A,I8)"), test_opt, " = ", fortdg(i)%iptr
                  case (2)
                     read (test_val, *) fortdg(i)%rptr
                     if (myproc == 0) print("(A,A,E21.8)"), test_opt, " = ", fortdg(i)%rptr
                  case (3)
                     fortdg(i)%cptr = trim(test_val)
                     if (myproc == 0) print("(A,A,A)"), test_opt, " = ", fortdg(i)%cptr
                  end select

                  found = 1 ! flag match
                  opt_read = opt_read + 1
                  fortdg(i)%flag = 1 ! flag option as found

                  exit test

               end if
            end do test

            if (myproc == 0) then
               if (found == 0 .and. eqind > 0) then
                  ! unmatched lines with an equal sign are either incorrect or no longer supported
                  print("(3A)"), "*** WARNING: ", test_opt, " is an incorrect or depreciated value ***"
               else if (found == 0) then
                  ! unmatched lines without an equal sign are ignored
                  print("(A)"), "*** WARNING: non-comment line does not contain a keyword assignment***"
               end if
            end if

         end if
      end do

      if (myproc == 0) print *, ""

      call CHECK_ERRORS(opt_read)

      if (myproc == 0) print *, ""
      close (25)

   end subroutine READ_KEYWORD_FORT_DG

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine CHECK_ERRORS(opt_read)

      use sizes, only: myproc

      implicit none

      integer :: i, j, opt
      integer :: opt_read
      integer :: quit

      if (opt_read /= nopt) then

         ! check for required options that are unspecifed
         quit = 0
         do opt = 1, nopt
            i = fortdg_ind(opt)
            if ((fortdg(i)%flag == 0) .and. (fortdg(i)%required == 1)) then
               quit = 1 ! flag fatal error
            end if
         end do

         if (quit == 1) then

            if (myproc == 0) then
               print("(A)"), "*** ERROR: There are missing required options in the fort.dg file ***"
               print("(A)"), "           The following options must be specified: "
               j = 0
               do opt = 1, nopt
                  i = fortdg_ind(opt)
                  if ((fortdg(i)%flag == 0) .and. (fortdg(i)%required == 1)) then
                     j = j + 1
                     print "(A,I3,2A)", "              ", j, ") ", fortdg(i)%key
                  end if
               end do

               print("(A)"), "!!!!!! EXECUTION WILL NOW BE TERMINATED !!!!!!"
            end if

            stop

         else

            if (myproc == 0) then
               print("(A)"), "*** WARNING: There are missing optional options in the fort.dg file ***"
               print("(A)"), "             The following default values will be used: "
               j = 0
               do opt = 1, nopt
                  i = fortdg_ind(opt)
                  if ((fortdg(i)%flag == 0) .and. (fortdg(i)%required == 0)) then

                     j = j + 1
                     select case (fortdg(i)%vartype)
                     case (1)
                        print("(A,I3,A,A,A,I8)"), "              ", j, ") ", fortdg(i)%key, " = ", fortdg(i)%iptr
                     case (2)
                        print("(A,I3,A,A,A,E21.8)"), "              ", j, ") ", fortdg(i)%key, " = ", fortdg(i)%rptr
                     case (3)
                        print("(A,I3,A,A,A,A)"), "              ", j, ") ", fortdg(i)%key, " = ", fortdg(i)%cptr
                     end select

                  end if
               end do

               print("(A)"), '!!!!!! EXECUTION WILL CONTINUE !!!!!!!!'
            end if

         end if

      end if

      return
   end subroutine CHECK_ERRORS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine FORT_DG_SETUP()

      ! Subroutine that configures the fort.dg options
      !
      !   This subroutine is meant to add flexibility in adding/depreciating
      !   features while maintaining forward (and some degree of backward) compatibility
      !
      !   - Options can be added to the fort.dg file by:
      !       1) Specifying a keyword in a unused index (<= maxopt) of the fortdg structure
      !       2) Associating the appropriate pointer with the corresponding variable
      !          Note: pointer must agree with the associated variable type
      !                (iptr=integer, rptr=real, cptr=character)
      !          Note: the associated variable must be declared using the TARGET attribute
      !       3) Specifying whether the variable is required (1 = yes, 0 = no)
      !       4) Providing a default value
      !
      !   - Options can be removed from the fort.dg file by:
      !       1) Commenting out or deleting an existing entry in the fortdg structure
      !          Note: re-indexing subsequent entries is not necessary (see fortdg(17) below)
      !
      !       OR
      !
      !       2) Setting the fortdg(i)%required variable to 0
      !
      !   - New features should be added as fortdg(i)%required = 0 as much as possible
      !     to maintain backward compatibility, older fort.dg files not containing these
      !     options will cause provided default values to be used (these should be set so
      !     the feature is turned off)
      !
      !   - fort.dg files containing new feature options can still be used for previous
      !     versions of the code because the new options will be ignored

      use global, only: dgswe, dg_to_cg, sedflag, reaction_rate, sed_equationX, sed_equationY
      use sizes, only: myproc, layers
      use dg, only: padapt, pflag, gflag, diorism, pl, ph, px, slimit, plimit, &
                    pflag2con1, pflag2con2, lebesgueP, fluxtype, rk_stage, rk_order, &
                    modal_ic, dghot, dghotspool, slopeflag, slope_weight, porosity, &
                    sevdm, mnes, artdif, kappa, s0, uniform_dif, tune_by_hand, rainfall

      implicit none

      integer :: i
      integer :: ncheck
      character(15) :: empty
      character(28) :: sedXdef, sedYdef

      ! initialize fortdg structure
      do i = 1, maxopt
         nullify (fortdg(i)%iptr)
         nullify (fortdg(i)%rptr)
         nullify (fortdg(i)%cptr)

         fortdg(i)%key = empty
         fortdg(i)%flag = 0
      end do

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Configure fort.dg options here:
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      sedXdef = "(ZE_ROE+bed_ROE)**-1 *QX_ROE"
      sedYdef = "(ZE_ROE+bed_ROE)**-1 *QY_ROE"

      !    keywords                         target variables                      requirement                 default values
      fortdg(1)%key = "dgswe"; fortdg(1)%iptr => dgswe; fortdg(1)%required = 1; fortdg(1)%iptr = 1
      fortdg(2)%key = "padapt"; fortdg(2)%iptr => padapt; fortdg(2)%required = 1; fortdg(2)%iptr = 0
      fortdg(3)%key = "pflag"; fortdg(3)%iptr => pflag; fortdg(3)%required = 1; fortdg(3)%iptr = 2
      fortdg(4)%key = "gflag"; fortdg(4)%iptr => gflag; fortdg(4)%required = 1; fortdg(4)%iptr = 1
      fortdg(5)%key = "dis_tol"; fortdg(5)%rptr => diorism; fortdg(5)%required = 1; fortdg(5)%rptr = 8
      fortdg(6)%key = "pl"; fortdg(6)%iptr => pl; fortdg(6)%required = 1; fortdg(6)%iptr = 1
      fortdg(7)%key = "ph"; fortdg(7)%iptr => ph; fortdg(7)%required = 1; fortdg(7)%iptr = 1
      fortdg(8)%key = "px"; fortdg(8)%iptr => px; fortdg(8)%required = 1; fortdg(8)%iptr = 1
      fortdg(9)%key = "slimit"; fortdg(9)%rptr => slimit; fortdg(9)%required = 1; fortdg(9)%rptr = 0.00005
      fortdg(10)%key = "plimit"; fortdg(10)%rptr => plimit; fortdg(10)%required = 1; fortdg(10)%rptr = 10
      fortdg(11)%key = "k"; fortdg(11)%rptr => pflag2con1; fortdg(11)%required = 1; fortdg(11)%rptr = 1
      fortdg(12)%key = "ks"; fortdg(12)%rptr => pflag2con2; fortdg(12)%required = 1; fortdg(12)%rptr = 0.5
      fortdg(13)%key = "L"; fortdg(13)%iptr => lebesgueP; fortdg(13)%required = 1; fortdg(13)%iptr = 2
      fortdg(14)%key = "fluxtype"; fortdg(14)%iptr => fluxtype; fortdg(14)%required = 1; fortdg(14)%iptr = 1
      fortdg(15)%key = "rk_stage"; fortdg(15)%iptr => rk_stage; fortdg(15)%required = 1; fortdg(15)%iptr = 2
      fortdg(16)%key = "rk_order"; fortdg(16)%iptr => rk_order; fortdg(16)%required = 1; fortdg(16)%iptr = 2
!       fortdg(17)%key = "dg_to_cg";      fortdg(17)%iptr => dg_to_cg;     fortdg(17)%required = 1;  fortdg(17)%iptr = 1
      fortdg(18)%key = "modal_ic"; fortdg(18)%iptr => modal_ic; fortdg(18)%required = 1; fortdg(18)%iptr = 0
      fortdg(19)%key = "dghot"; fortdg(19)%iptr => dghot; fortdg(19)%required = 1; fortdg(19)%iptr = 0
      fortdg(20)%key = "dghotspool"; fortdg(20)%iptr => dghotspool; fortdg(20)%required = 1; fortdg(20)%iptr = 86400
      fortdg(21)%key = "slopeflag"; fortdg(21)%iptr => slopeflag; fortdg(21)%required = 1; fortdg(21)%iptr = 5
      fortdg(22)%key = "weight"; fortdg(22)%rptr => slope_weight; fortdg(22)%required = 1; fortdg(22)%rptr = 1
      fortdg(23)%key = "sedflag"; fortdg(23)%iptr => sedflag; fortdg(23)%required = 1; fortdg(23)%iptr = 0
      fortdg(24)%key = "porosity"; fortdg(24)%rptr => porosity; fortdg(24)%required = 1; fortdg(24)%rptr = 0.0001
      fortdg(25)%key = "sevdm"; fortdg(25)%rptr => sevdm; fortdg(25)%required = 1; fortdg(25)%rptr = 0.00001
      fortdg(26)%key = "layers"; fortdg(26)%iptr => layers; fortdg(26)%required = 0; fortdg(26)%iptr = 1
      fortdg(27)%key = "rxn_rate"; fortdg(27)%rptr => reaction_rate; fortdg(27)%required = 1; fortdg(27)%rptr = 1.0
      fortdg(28)%key = "nelem"; fortdg(28)%iptr => mnes; fortdg(28)%required = 1; fortdg(28)%iptr = 23556
      fortdg(29)%key = "artdif"; fortdg(29)%iptr => artdif; fortdg(29)%required = 1; fortdg(29)%iptr = 0
      fortdg(30)%key = "kappa"; fortdg(30)%rptr => kappa; fortdg(30)%required = 1; fortdg(30)%rptr = -1.0
      fortdg(31)%key = "s0"; fortdg(31)%rptr => s0; fortdg(31)%required = 1; fortdg(31)%rptr = 0.0
      fortdg(32)%key = "uniform_dif"; fortdg(32)%rptr => uniform_dif; fortdg(32)%required = 1; fortdg(32)%rptr = 2.5e-6
      fortdg(33)%key = "tune_by_hand"; fortdg(33)%iptr => tune_by_hand; fortdg(33)%required = 1; fortdg(33)%iptr = 0
      fortdg(34)%key = "sed_equationX"; fortdg(34)%cptr => sed_equationX; fortdg(34)%required = 0; fortdg(34)%cptr = sedXdef
      fortdg(35)%key = "sed_equationY"; fortdg(35)%cptr => sed_equationY; fortdg(35)%required = 0; fortdg(35)%cptr = sedYdef
      fortdg(36)%key = "rainfall"; fortdg(36)%iptr => rainfall; fortdg(36)%required = 0; fortdg(36)%iptr = 0

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! End configuration
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      nopt = 0
      ncheck = 0
      do i = 1, maxopt

         ! find and keep track of populated indicies
         if (fortdg(i)%key /= empty) then
            nopt = nopt + 1
            fortdg_ind(nopt) = i
         end if

         ! determine target variable type by checking association status
         fortdg(i)%vartype = 0

         if (associated(fortdg(i)%iptr)) then ! integer
            ncheck = ncheck + 1
            fortdg(i)%vartype = 1
         end if

         if (associated(fortdg(i)%rptr)) then ! real
            ncheck = ncheck + 1
            fortdg(i)%vartype = 2
         end if

         if (associated(fortdg(i)%cptr)) then ! character
            ncheck = ncheck + 1
            fortdg(i)%vartype = 3
         end if
      end do

!       PRINT*, "Number of options = ", nopt
!       PRINT*, "Number of pointer associations = ", ncheck

      ! ensure user has associated each keyword pointer
      if (nopt /= ncheck) then
         if (myproc == 0) then
            print("(A)"), "*** ERROR: fort.dg option pointer association error ***"
            print("(A)"), "           check keyword configuration in fort_dg_setup subroutine"
         end if

         stop
      end if

      return
   end subroutine FORT_DG_SETUP

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module FORT_DG
