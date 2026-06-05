!C----------------------------------------------------------------------------
!C
!C                           MODULE DECOMP
!C
!C----------------------------------------------------------------------------
!C
!C                  For use with ADCPREP Version 2.0a (  7/07/05 )
!C
!C                     current for hp_DG_ADCIRC v8.1-mp   7/07/2005
!C----------------------------------------------------------------------------
!C
subroutine DECOMP()

   use pre_global

   !IMPLICIT NONE
   !C
   !C--------------------------------------------------------------------------C
   !C                    (  Serial Version 1.1  5/04/99  )                     C
   !C                                                                          C
   !C  Decomposes the ADCIRC grid into NPROC subdomains.                       C
   !C  The Decomposition Variables are defined in the include file adcprep.inc C
   !C  This version is compatible with ADCIRC version 34.03                    C
   !C                                                                          C
   !C  12/14/98 vjp  Added interface to METIS 4.0                              C
   !C   3/10/99 vjp  Rewritten to allow Weir-node pairs to both be ghost nodes C
   !C   4/05/99 vjp  Fixed bugs in metis interface routine                     C
   !C  ??/??/05 sb   Accomodate hp_DG_ADCIRC                                   C
   !C  02/26/06 sb   Speed up decomp process                                   C
   !C                                                                          C
   !C--------------------------------------------------------------------------C
   !C
   integer :: N1, N2, N3, VTMAX
   integer :: I, J, K, L, ITEMP, ITEMP2, IPR
   !Csb-DG1
   integer :: M1, M2, A, B
   integer :: INEIELG, NELPI
   integer, allocatable :: COMM_PE_FLAG(:)
   !C--
   !Csb--02/26/06
   integer :: NM1, NM2, NM3, IK, JJ
   !C--
   integer :: ITOT, IEL, IELG
   integer :: IG1, PE1, PE2, PE3
   integer :: NCOUNT
   integer :: INDEX, INDEX2
   integer, allocatable :: ITVECT(:)
   character(8) :: PE
   character(6) :: GLOBAL

   !Cmm   added 5/19/09
   !integer :: nn1, nn2, nn3
   !integer :: npe1, npe2, npe3

   !C
   VTMAX = 24*MNP
   allocate (ITVECT(VTMAX))

   !C   STEP 1:
   !C-- Use Partition of nodes to compute the number of Resident Nodes
   !C   to be assigned to each processor.
   !C-- Then construct Local-to-Global and Global-to-Local Node maps
   !C   for resident nodes: IMAP_NOD_LG(I,PE),  IMAP_NOD_GL(1:2,I)

   do I = 1, NPROC ! Use METIS 4.0  Partition
      NOD_RES_TOT(I) = 0
   end do
   do J = 1, NNODG
      NCOUNT = NOD_RES_TOT(PROC(J)) + 1
      IMAP_NOD_GL(1, J) = PROC(J)
      IMAP_NOD_GL(2, J) = NCOUNT
      IMAP_NOD_LG(NCOUNT, PROC(J)) = J
      NOD_RES_TOT(PROC(J)) = NCOUNT
   end do
   !c      DO I = 1, NNODG
   !c         if (proc(i) == 4) then
   !c            write(7,*) I, proc(I)
   !c         endif
   !c      ENDDO
   !c      do i = 1, nnodg
   !c         if (proc(i) == 8) then
   !c            write(7,*) i, proc(i)
   !c         endif
   !c      enddo

   !c     cmm added test here
   !c      do i = 1,nproc
   !c         write(*,*) "NOD_RES_TOT(",i,") = ", nod_res_tot(i)
   !c      enddo

   !C STEP 2:
   !C  Construct Local-to-Global Element Map: IMAP_EL_LG(:,PE)
   !C  Add an element to the map if it has an resident node

   do I = 1, NPROC
      NELP(I) = 0
      do K = 1, NELG
         N1 = NNEG(1, K)
         N2 = NNEG(2, K)
         N3 = NNEG(3, K)
         PE1 = IMAP_NOD_GL(1, N1) ! Is any vertex a resident node?
         PE2 = IMAP_NOD_GL(1, N2)
         PE3 = IMAP_NOD_GL(1, N3)
         if ((PE1 == I) .or. (PE2 == I) .or. (PE3 == I)) then
            NELP(I) = NELP(I) + 1
            IMAP_EL_LG(NELP(I), I) = K
         end if
      end do
      if (NELP(I) > MNEP) stop 'NELP(I) > MNEP'
   end do
   !C cmm added test here
   !c      write(*,*) "CMM test:"
   !c      do i = 1,nproc
   !c         write(*,*) "NELP(",i,") = ", NELP(i)
   !c      enddo

   !C STEP 3:
   !C--Using Local-to-Global Element map
   !C  Construct Local-to-Global Node map:  IMAP_NOD_LG(I,PE)
   !C  and reconstruct Global-to-Local map for resident nodes
   !C
   do I = 1, NPROC
      ITOT = 0
      do J = 1, NELP(I)
         IEL = IMAP_EL_LG(J, I)
         do L = 1, 3
            ITOT = ITOT + 1
            ITVECT(ITOT) = NNEG(L, IEL)
         end do
      end do
      ITEMP = ITOT
      if (ITOT > VTMAX) stop 'step3 decomp'
      call SORT(ITEMP, ITVECT) ! Sort and remove multiple occurrences
      ITOT = 1
      IMAP_NOD_LG(1, I) = ITVECT(1)
      if (IMAP_NOD_GL(1, ITVECT(1)) == I) then
         IMAP_NOD_GL(2, ITVECT(1)) = 1
      end if
      do J = 2, ITEMP
         if (ITVECT(J) /= IMAP_NOD_LG(ITOT, I)) then
            ITOT = ITOT + 1
            IMAP_NOD_LG(ITOT, I) = ITVECT(J)
            if (IMAP_NOD_GL(1, ITVECT(J)) == I) then
               IMAP_NOD_GL(2, ITVECT(J)) = ITOT
            end if
         end if
      end do
      NNODP(I) = ITOT
      if (NNODP(I) > MNPP) stop 'NNODP > MNPP'
   end do
   !c     print *, "Number of Nodes Assigned to PEs"
   !c     DO I=1, NPROC
   !c        print *, I-1, NNODP(I)
   !c        DO J=1, NNODP(I)
   !c           print *, J,IMAP_NOD_LG(J,I)
   !c        ENDDO
   !c     ENDDO

   !C STEP 4:
   !C--If there are any global Weir-node pairs, construct
   !C  Local-to-Global Weir Node maps: WEIRP_LG(:,PE), WEIRDP_LG(:,PE)
   !C  Rule: if a global Weir node is assigned ( either as a resident or ghost node )
   !C        then make it and its dual a local Weir-node pair

   if (NWEIR > 0) then
      ITEMP2 = 0 ! Seizo Add
      do I = 1, NPROC
         ITOT = 0
         ! Seizo 2008.07.14
         do J = 1, NWEIR
            INDEX = WEIR(J)
            INDEX2 = WEIRD(J)
            do K = 1, NNODP(I)
               N1 = IMAP_NOD_LG(K, I)
               if ((N1 == INDEX) .or. (N1 == INDEX2)) then
                  ITOT = ITOT + 1
                  ITVECT(ITOT) = J
                  exit
               end if
            end do
         end do
         !
         ! Seizo     DO J = 1,NNODP(I)
         ! Seizo        CALL SEARCH(WEIR,NWEIR,IMAP_NOD_LG(J,I),INDEX)
         ! Seizo        IF (INDEX.NE.0) THEN
         ! Seizo          ITOT = ITOT+1
         ! Seizo          ITVECT(ITOT) = INDEX
         ! Seizo        ENDIF
         ! Seizo        CALL SEARCH(WEIRD,NWEIR,IMAP_NOD_LG(J,I),INDEX2)
         ! Seizo        IF (INDEX2.NE.0) THEN
         ! Seizo          ITOT = ITOT+1
         ! Seizo          ITVECT(ITOT) = INDEX2
         ! Seizo        ENDIF
         ! Seizo     ENDDO
         NWEIRP(I) = 0
         ITEMP = ITOT
         if (ITOT > VTMAX) stop 'step4 decomp'
         if (ITEMP > 1) then
            call SORT(ITEMP, ITVECT)
            ITOT = 1
            INDEX = ITVECT(1)
            WEIRP_LG(ITOT + ITEMP2) = WEIR(INDEX)
            WEIRDP_LG(ITOT + ITEMP2) = WEIRD(INDEX)
            ! Seizo      WEIRP_LG(ITOT,I)  = WEIR(INDEX)
            ! Seizo      WEIRDP_LG(ITOT,I) = WEIRD(INDEX)
            do J = 2, ITEMP
               if (ITVECT(J) /= INDEX) then
                  INDEX = ITVECT(J)
                  ITOT = ITOT + 1
                  WEIRP_LG(ITOT + ITEMP2) = WEIR(INDEX)
                  WEIRDP_LG(ITOT + ITEMP2) = WEIRD(INDEX)
                  ! Seizo           WEIRP_LG(ITOT,I)  = WEIR(INDEX)
                  ! Seizo           WEIRDP_LG(ITOT,I) = WEIRD(INDEX)
               end if
            end do
            NWEIRP(I) = ITOT
            ITEMP2 = ITEMP2 + ITOT ! Seizo Add
         end if
         !c       DO J = 1, NWEIRP(I)
         !c          print *, J, WEIRP_LG(J,I),WEIRDP_LG(J,I)
         !c       ENDDO
         !c       print *, "decomp: Number of WEIR node-pairs on PE",I-1,
         !c    &           " = ",NWEIRP(I)
      end do
   else
      do I = 1, NPROC
         NWEIRP(I) = 0
         !c          print *, "decomp: Number of WEIR node-pairs on PE",I-1,
         !c    &              " = ",NWEIRP(I)
      end do
   end if

   !C STEP 5:
   !C--If there are any global Weir-node pairs,
   !C  Re-construct Local-to-Global Element Map: IMAP_EL_LG(:,PE)
   !C  Rule:  Add an element if it has an resident node or
   !C         has the dual Weir node of a resident or ghost node

   if (NWEIR > 0) then
      ITEMP2 = 1 ! Seizo Add
      do I = 1, NPROC
         NELP(I) = 0
         !c          print *, "PE = ",I-1
         do K = 1, NELG
            N1 = NNEG(1, K)
            N2 = NNEG(2, K)
            N3 = NNEG(3, K)
            PE1 = IMAP_NOD_GL(1, N1) ! Is any vertex a resident node?
            PE2 = IMAP_NOD_GL(1, N2)
            PE3 = IMAP_NOD_GL(1, N3) ! belong to a Weir-node pair ?
            call SEARCH3(WEIRP_LG(ITEMP2), NWEIRP(I), N1, N2, N3, INDEX)
            call SEARCH3(WEIRDP_LG(ITEMP2), NWEIRP(I), N1, N2, N3, INDEX2)
            ! Seizo       CALL SEARCH3(WEIRP_LG(1,I),NWEIRP(I),N1,N2,N3,INDEX)
            ! Seizo       CALL SEARCH3(WEIRDP_LG(1,I),NWEIRP(I),N1,N2,N3,INDEX2)
            if ((PE1 == I) .or. (PE2 == I) .or. (PE3 == I) &
                .or. (INDEX /= 0) .or. (INDEX2 /= 0)) then
               !c                print *, K, PE1,PE2,PE3,INDEX,INDEX2
               NELP(I) = NELP(I) + 1
               IMAP_EL_LG(NELP(I), I) = K
            end if
         end do
         if (NELP(I) > MNEP) stop 'NELP(I) > MNEP'
         !c       print *, "Number of elements on PE",I-1," = ",NELP(I)
         !c       DO J = 1, NELP(I)
         !c          print *, J, IMAP_EL_LG(J,I)
         !c       ENDDO
         ITEMP2 = ITEMP2 + NWEIRP(I) ! Seizo Add
      end do
   end if

   !C STEP 5-DG-1:
   !C--Alloc memory 1
   !C
   call ALLOC_DG1()

   !C STEP 5-DG-2
   !C--Count maximum of the number of the elements associated with a node
   NNDEL(:) = 0
   do IK = 1, MNE
      NNDEL(NNEG(1, IK)) = NNDEL(NNEG(1, IK)) + 1
      NNDEL(NNEG(2, IK)) = NNDEL(NNEG(2, IK)) + 1
      NNDEL(NNEG(3, IK)) = NNDEL(NNEG(3, IK)) + 1
   end do
   MNNDEL = 0
   do IK = 1, MNP
      if (NNDEL(IK) > MNNDEL) MNNDEL = NNDEL(IK)
   end do

   !C STEP 5-DG-3:
   !C--Alloc memory 1B
   !C
   call ALLOC_DG1B()

   !C STEP 5-DG-4:
   !C--Make node-to-elements table
   NNDEL(:) = 0

   do IK = 1, MNE
      NM1 = NNEG(1, IK)
      NM2 = NNEG(2, IK)
      NM3 = NNEG(3, IK)

      NNDEL(NM1) = NNDEL(NM1) + 1
      NDEL(NM1, NNDEL(NM1)) = IK

      NNDEL(NM2) = NNDEL(NM2) + 1
      NDEL(NM2, NNDEL(NM2)) = IK

      NNDEL(NM3) = NNDEL(NM3) + 1
      NDEL(NM3, NNDEL(NM3)) = IK
   end do

   !C STEP 5-DG-5:
   !C--Make a index of the neighboring element of each edge
   !C
   IMAP_NEIGHEDG(:, :, :) = 0

   do I = 1, NELG
      do A = 1, 3 !cmm removed 2001 label
         if (IMAP_NEIGHEDG(1, A, I) /= 0) cycle !GOTO 2001

         N1 = NNEG(mod(A + 0, 3) + 1, I) ! Endnode 1 of edge A
         N2 = NNEG(mod(A + 1, 3) + 1, I) ! Endnode 2 of edge A

         do JJ = 1, NNDEL(N1)
            J = NDEL(N1, JJ)
            if (I /= J) then
               do B = 1, 3
                  M1 = NNEG(mod(B + 0, 3) + 1, J) ! Endnode 1 of edge B
                  M2 = NNEG(mod(B + 1, 3) + 1, J) ! Endnode 2 of edge B

                  if (N1 == M2 .and. N2 == M1) then
                     IMAP_NEIGHEDG(1, A, I) = J
                     IMAP_NEIGHEDG(2, A, I) = B

                     IMAP_NEIGHEDG(1, B, J) = I
                     IMAP_NEIGHEDG(2, B, J) = A

                     goto 2001
                  end if
               end do
            end if
         end do
2001     continue
      end do
   end do

   !C STEP 5-DG-6:
   !C--If there are any elements shared by three different sub-domains,
   !C  Re-construct Local-to-Global Element Map: IMAP_EL_LG(:,PE)
   !C  Rule: Add non-local elements adjacent to the elements
   !C        shared by three different sub-domains
   !C        to the local element list
   !C
   do I = 1, NPROC
      NELPI = NELP(I)
      do J = 1, NELPI
         IELG = IMAP_EL_LG(J, I)

         N1 = NNEG(1, IELG)
         N2 = NNEG(2, IELG)
         N3 = NNEG(3, IELG)

         PE1 = IMAP_NOD_GL(1, N1)
         PE2 = IMAP_NOD_GL(1, N2)
         PE3 = IMAP_NOD_GL(1, N3)

         if ((PE1 /= PE2) .and. (PE2 /= PE3) .and. (PE3 /= PE1)) then
            if (PE1 == I) then
               K = 1
            else if (PE2 == I) then
               K = 2
            else if (PE3 == I) then
               K = 3
            end if

            INEIELG = IMAP_NEIGHEDG(1, K, IELG)

            !Cmm   Do NOT add the neighboring element if the current
            ! element is on the global boundary.
            if (ineielg /= 0) then

               !Cmm   Also do NOT add the neighboring element if it
               ! is already a member of the current subdomain.
               !c                  if (I == 4 .or. i == 8) then
               !c                     write(*,*) "NELP(I) = ", NELP(I)
               !c                     write(*,*) "NELPI = ", NELPI
               !c                  endif
               do L = 1, NELP(I)
                  if (IMAP_EL_LG(L, I) == INEIELG) then
                     !Cmm   As much as I cringe at using GO TO's,
                     ! I feel one is necessary here
                     !     since we need to break out of a nested loop
                     ! and give up on element J.
                     !if (i == 4 .or. i == 8) then
                     !   write(*,*) "Found a match: L = ", L
                     !   write(*,*) "INEIELG = ", INEIELG
                     !endif
                     go to 520
                  end if
               end do

               if (NELP(I) >= MNEP) then
                  stop 'FATAL ERROR (STEP 5-DG-6 (2))'
               else
                  NELP(I) = NELP(I) + 1
                  IMAP_EL_LG(NELP(I), I) = INEIELG
               end if
               !c               endif
            end if
         end if
520      continue
      end do
   end do
   !C--

   !c     cmm added test here
   !c      write(*,*) "NELP(",4,") = ", NELP(4)
   !c      write(*,*) "NELP(",8,") = ", NELP(8)

   !C STEP 6:
   !C--Using Local-to-Global Element map
   !C  Construct Local-to-Global Node map:  IMAP_NOD_LG(I,PE)
   !C  and reconstruct Global-to-Local map for resident nodes
   !C
   do I = 1, NPROC
      ITOT = 0
      do J = 1, NELP(I)
         IEL = IMAP_EL_LG(J, I)
         do L = 1, 3
            ITOT = ITOT + 1
            ITVECT(ITOT) = NNEG(L, IEL)
         end do
      end do
      ITEMP = ITOT
      if (ITOT > VTMAX) stop 'step6 decomp'
      call SORT(ITEMP, ITVECT) ! Sort and remove multiple occurrences
      ITOT = 1
      IMAP_NOD_LG(1, I) = ITVECT(1)
      if (IMAP_NOD_GL(1, ITVECT(1)) == I) then
         IMAP_NOD_GL(2, ITVECT(1)) = 1
      end if
      do J = 2, ITEMP
         if (ITVECT(J) /= IMAP_NOD_LG(ITOT, I)) then
            ITOT = ITOT + 1
            IMAP_NOD_LG(ITOT, I) = ITVECT(J)
            if (IMAP_NOD_GL(1, ITVECT(J)) == I) then
               IMAP_NOD_GL(2, ITVECT(J)) = ITOT
            end if
         end if
      end do
      NNODP(I) = ITOT
      if (NNODP(I) > MNPP) stop 'NNODP > MNPP'
   end do
   !c     print *, "Number of Nodes Assigned to PEs"
   !c     DO I=1, NPROC
   !c        print *, I-1, NNODP(I)
   !c        DO J=1, NNODP(I)
   !c           print *, J,IMAP_NOD_LG(J,I)
   !c        ENDDO
   !c     ENDDO

   !C STEP 7:
   !C--Construct Local Element Connectivity Table for each PE: NNEP(3,I,PE)
   !C
   do I = 1, NPROC
      ITEMP = NNODP(I)
      do J = 1, NNODP(I)
         ITVECT(J) = IMAP_NOD_LG(J, I)
      end do
      do J = 1, NELP(I)
         IELG = IMAP_EL_LG(J, I)
         do L = 1, 3
            IG1 = NNEG(L, IELG)
            call LOCATE(ITVECT, ITEMP, IG1, K)
            if (K <= 0) then
               if (IMAP_NOD_LG(K + 1, I) == IG1) then
                  NNEP(L, J, I) = K + 1
               else
                  stop 'ERROR IN IMAP_NOD_LG'
               end if
            elseif (K >= NNODP(I)) then
               if (IMAP_NOD_LG(K, I) == IG1) then
                  NNEP(L, J, I) = K
               else
                  stop 'ERROR IN IMAP_NOD_LG'
               end if
            else
               if (IMAP_NOD_LG(K, I) == IG1) then
                  NNEP(L, J, I) = K
               elseif (IMAP_NOD_LG(K + 1, I) == IG1) then
                  NNEP(L, J, I) = K + 1
               else
                  stop 'ERROR IN IMAP_NOD_LG'
               end if
            end if
         end do
      end do
   end do

   !C STEP 8:
   !C--Compute the number of communicating PEs and their
   !C  list for each PE:  NUM_COMM_PE(PE) and COMM_PE_NUM(IPE,PE)
   !C
   do I = 1, NPROC
      NUM_COMM_PE(I) = 0
      ITEMP = 0
      do J = 1, NNODP(I)
         INDEX = IMAP_NOD_LG(J, I)
         IPR = IMAP_NOD_GL(1, INDEX)
         if (IPR /= I) then
            ITEMP = ITEMP + 1
            ITVECT(ITEMP) = IPR
         end if
      end do
      if (ITEMP == 0) then
         NUM_COMM_PE(I) = 0
      else
         if (ITEMP > VTMAX) stop 'step8 decomp'
         call SORT(ITEMP, ITVECT)
         COMM_PE_NUM(1, I) = ITVECT(1)
         ITOT = 1
         do J = 1, ITEMP
            if (ITVECT(J) /= COMM_PE_NUM(ITOT, I)) then
               ITOT = ITOT + 1
               COMM_PE_NUM(ITOT, I) = ITVECT(J)
            end if
         end do
         NUM_COMM_PE(I) = ITOT
         if (NUM_COMM_PE(I) > MNPROC) stop 'NUM_COMM_PE>MNPROC'
      end if
   end do

   !C STEP 9:
   !C--Construct a Global-to-Local node mapping: IMAP_NOD_GL2(*,J)
   !C  This is not a function, but is rather a relation
   !C  It works for both resident and ghost nodes

   do I = 1, NNODG
      ITOTPROC(I) = 0
   end do
   do I = 1, NPROC
      do J = 1, NNODP(I)
         INDEX = IMAP_NOD_LG(J, I)
         ITOTPROC(INDEX) = ITOTPROC(INDEX) + 1
         if (ITOTPROC(INDEX) > MNPROC) then
            write (6, *) 'Some nodes belong to more processors', &
               ' than MNPROC'
            stop
         end if
         ITEMP = (ITOTPROC(INDEX) - 1)*2 + 1
         IMAP_NOD_GL2(ITEMP, INDEX) = I
         IMAP_NOD_GL2(ITEMP + 1, INDEX) = J
      end do
   end do
   !c     print *, "Global Nodes assigned to more than one PE"
   !c     do J=1, NNODG
   !c        if (ITOTPROC(J).GT.1) print *, J, ITOTPROC(J)
   !c     enddo

   !Csb-DG2
   !C STEP 9-DG-1:
   !C--Compute NIEL_SEND(:) and NIEL_RECV(:) to prepare for ALLOC_DG2
   !C
   NIEL_SEND(:) = 0
   NIEL_RECV(:) = 0

   do I = 1, NPROC
      do J = 1, NELP(I)
         IELG = IMAP_EL_LG(J, I)

         N1 = NNEG(1, IELG)
         N2 = NNEG(2, IELG)
         N3 = NNEG(3, IELG)

         PE1 = IMAP_NOD_GL(1, N1)
         PE2 = IMAP_NOD_GL(1, N2)
         PE3 = IMAP_NOD_GL(1, N3)

         ! Element J in sub-domain I should receive information from
         ! a neighboring sub-domain (i.e. element J is a ghost element),
         ! if two or three nodes of the element belong to a sub-domain
         ! other than sub-domain I.
         ! However, there is an exception.  If the elemental nodes belong
         ! to three different sub-domains and one of the nodes belong to
         ! sub-domain I, this element is not a ghost element because
         ! it is immersed in sub-domain I since an adjacent element has
         ! been added at STEP 5-DG-3.

         if ((PE1 == I .and. PE2 == I) .or. &
             (PE2 == I .and. PE3 == I) .or. &
             (PE3 == I .and. PE1 == I)) then
            ! A resident element; Do nothing
         else if ((PE1 /= PE2 .and. PE2 /= PE3 .and. PE3 /= PE1) .and. &
                  (PE1 == I .or. PE2 == I .or. PE3 == I)) then
            ! An element shared by three different sub-domains; Do nothing
         else if (PE1 /= PE2 .and. PE2 /= PE3 .and. PE3 /= PE1) then
            ! An element added at STEP 5-DG-3

            ! Ghost element J can receive information from any of the
            ! neighboring sub-domains if the elemental nodes belong
            ! to three different sub-domains and none of the sub-domains
            ! is sub-domain I.

            NIEL_RECV(I) = NIEL_RECV(I) + 1
            NIEL_SEND(PE1) = NIEL_SEND(PE1) + 1

         else if ((PE1 /= I .and. PE2 /= I) .or. &
                  (PE2 /= I .and. PE3 /= I) .or. &
                  (PE3 /= I .and. PE1 /= I)) then
            ! An element which comes here should have one node which belongs
            ! to sub-domain I and two nodes which belong to one sub-domain
            ! different from sub-domain I

            ! Ghost element J can receive information from a neighboring
            ! sub-domain to which two of the elemental nodes belong to the
            ! sub-domain.

            NIEL_RECV(I) = NIEL_RECV(I) + 1

            if (PE1 == PE2) then
               NIEL_SEND(PE1) = NIEL_SEND(PE1) + 1
            else if (PE2 == PE3) then
               NIEL_SEND(PE2) = NIEL_SEND(PE2) + 1
            else if (PE3 == PE1) then
               NIEL_SEND(PE3) = NIEL_SEND(PE3) + 1
            else
               stop 'STEP 9-DG-1 You should not see this. (1)'
            end if
         else
            stop 'STEP 9-DG-1 You should not see this. (2)'
         end if
      end do
   end do

   !C STEP 9-DG-2:
   !C--Alloc memory 2
   !C
   call ALLOC_DG2()

   !C STEP 9-DG-3:
   !C--Decompose sub-domain interface edges
   !C
   ! Zero out again
   NIEL_RECV(:) = 0
   NIEL_SEND(:) = 0

   do I = 1, NPROC
      do J = 1, NELP(I)
         IELG = IMAP_EL_LG(J, I)

         N1 = NNEG(1, IELG)
         N2 = NNEG(2, IELG)
         N3 = NNEG(3, IELG)

         PE1 = IMAP_NOD_GL(1, N1)
         PE2 = IMAP_NOD_GL(1, N2)
         PE3 = IMAP_NOD_GL(1, N3)

         ! Element J in sub-domain I should receive information from
         ! a neighboring sub-domain (i.e. element J is a ghost element),
         ! if two or three nodes of the element belong to a sub-domain
         ! other than sub-domain I.
         ! However, there is an exception.  If the elemental nodes belong
         ! to three different sub-domains and one of the nodes belong to
         ! sub-domain I, this element is not a ghost element because
         ! it is immersed in sub-domain I since an adjacent element has
         ! been added at STEP 5-DG-3.

         if ((PE1 == I .and. PE2 == I) .or. &
             (PE2 == I .and. PE3 == I) .or. &
             (PE3 == I .and. PE1 == I)) then
            ! A resident element; Do nothing
         else if ((PE1 /= PE2 .and. PE2 /= PE3 .and. PE3 /= PE1) .and. &
                  (PE1 == I .or. PE2 == I .or. PE3 == I)) then
            ! An element shared by three different sub-domains; Do nothing
         else if (PE1 /= PE2 .and. PE2 /= PE3 .and. PE3 /= PE1) then
            ! An element added at STEP 5-DG-3

            ! Ghost element J can receive information from any of the
            ! neighboring sub-domains if the elemental nodes belong
            ! to three different sub-domains and none of the sub-domains
            ! is sub-domain I.

            NIEL_RECV(I) = NIEL_RECV(I) + 1
            IEL_RECV(1, NIEL_RECV(I), I) = PE1
            IEL_RECV(2, NIEL_RECV(I), I) = J

            NIEL_SEND(PE1) = NIEL_SEND(PE1) + 1
            IEL_SEND(1, NIEL_SEND(PE1), PE1) = I
            IEL_SEND(2, NIEL_SEND(PE1), PE1) = IELG

         else if ((PE1 /= I .and. PE2 /= I) .or. &
                  (PE2 /= I .and. PE3 /= I) .or. &
                  (PE3 /= I .and. PE1 /= I)) then
            ! An element which comes here should have one node which belongs
            ! to sub-domain I and two nodes which belong to one sub-domain
            ! different from sub-domain I

            ! Ghost element J can receive information from a neighboring
            ! sub-domain to which two of the elemental nodes belong to the
            ! sub-domain.

            if (PE1 == PE2) then
               K = PE1
            else if (PE2 == PE3) then
               K = PE2
            else if (PE3 == PE1) then
               K = PE3
            else
               stop 'STEP 9-DG-1 You should not see this. (1)'
            end if

            NIEL_RECV(I) = NIEL_RECV(I) + 1
            IEL_RECV(1, NIEL_RECV(I), I) = K
            IEL_RECV(2, NIEL_RECV(I), I) = J

            NIEL_SEND(K) = NIEL_SEND(K) + 1
            IEL_SEND(1, NIEL_SEND(K), K) = I
            IEL_SEND(2, NIEL_SEND(K), K) = IELG

         else
            stop 'STEP 9-DG-1 You should not see this. (2)'
         end if
      end do
   end do

   !C STEP 9-DG-4:
   !C--Replace the temporarily assigned global element IDs with local ones
   !C
   do I = 1, NPROC
      do J = 1, NIEL_SEND(I)
         IELG = IEL_SEND(2, J, I)

         do K = 1, NELP(I)
            if (IMAP_EL_LG(K, I) == IELG) then
               IEL_SEND(2, J, I) = K
               cycle
            end if
         end do
      end do
   end do

   !C STEP DG-7:
   !C--Compute NUM_COMM_PE_SEND and NUM_COMM_PE_RECV
   !C
   allocate (COMM_PE_FLAG(MNPROC))

   NUM_COMM_PE_RECV(:) = 0
   NUM_COMM_PE_SEND(:) = 0

   do I = 1, MNPROC
      ! RECV
      COMM_PE_FLAG(:) = 0

      do J = 1, NIEL_RECV(I)
         COMM_PE_FLAG(IEL_RECV(1, J, I)) = 1
      end do

      do J = 1, MNPROC
         if (COMM_PE_FLAG(J) == 1) then
            NUM_COMM_PE_RECV(I) = NUM_COMM_PE_RECV(I) + 1
         end if
      end do

      ! SEND
      COMM_PE_FLAG(:) = 0

      do J = 1, NIEL_SEND(I)
         COMM_PE_FLAG(IEL_SEND(1, J, I)) = 1
      end do

      do J = 1, MNPROC
         if (COMM_PE_FLAG(J) == 1) then
            NUM_COMM_PE_SEND(I) = NUM_COMM_PE_SEND(I) + 1
         end if
      end do

   end do

   deallocate (COMM_PE_FLAG)
   !C--

   !C STEP 10:
   !C--Print Summary of Decomposition
   !C
   print *, "Decomposition Data"
   print *, "DOMAIN  RES_NODES  GHOST_NODES  TOT_NODES  ELEMENTS"
   print *, "------  ---------  -----------  ---------  --------"
   GLOBAL = "GLOBAL"
   write (*, 90) GLOBAL, NNODG, NELG
   do I = 1, NPROC
      PE(1:8) = 'pe000000'
      call IWRITE(PE, 3, 8, I - 1)
      write (6, 92) PE, NOD_RES_TOT(I), NNODP(I) - NOD_RES_TOT(I), &
         NNODP(I), NELP(I)
   end do
90 format(1x, A6, 25x, I9, 2x, I9)
92 format(1x, A8, 1x, I9, 2x, I9, 4x, I9, 2x, I9)
   !C
   !RETURN
end subroutine DECOMP

subroutine DOMSIZE()

   use pre_global
   !C
   !C--------------------------------------------------------------------------C
   !C                  (  Serial Version 1.0  12/20/99 vjp )                   C
   !C                                                                          C
   !C  Takes dry run through the domain decomp to determine the max number of  C
   !C  nodes and elements assigned to any subdomain to determine MNPP and MNEP.C
   !C                                                                          C
   !C--------------------------------------------------------------------------C
   !C
   integer :: N1, N2, N3, VTMAX
   integer :: J, K, L, ITEMP
   integer :: ITOT, IEL, IPROC
   integer :: PE1, PE2, PE3
   integer :: INDEX, INDEX2
   integer :: RESNODE, NODES, NELEM, ONELEM, NLWEIR
   !Csb-
   integer :: NEL_ADDED
   !C--
   !C
   integer, allocatable :: ITVECT(:)
   integer, allocatable :: NODE_LG(:)
   integer, allocatable :: NODE_GL1(:)
   integer, allocatable :: NODE_GL2(:)
   integer, allocatable :: ELEM_LG(:)
   integer, allocatable :: LWEIR_LG(:), LWEIRD_LG(:)
   !C
   VTMAX = 24*MNP
   allocate (ITVECT(VTMAX))
   allocate (NODE_LG(MNP))
   allocate (NODE_GL1(MNP))
   allocate (NODE_GL2(MNP))
   allocate (LWEIR_LG(MNP), LWEIRD_LG(MNP))
   allocate (ELEM_LG(MNE))
   !C
   MNPP = 0
   MNEP = 0
   !C
   do IPROC = 1, NPROC
      !C
      !C   STEP 1:
      !C-- Use Partition of nodes to compute the number of Resident Nodes
      !C   to be assigned to each processor.
      !C-- Then construct Local-to-Global and Global-to-Local Node maps
      !C   for resident nodes
      !C
      !C
      do J = 1, MNP
         NODE_GL1(J) = 0
         NODE_GL2(J) = 0
         LWEIR_LG(J) = 0
         LWEIRD_LG(J) = 0
      end do
      !C
      do J = 1, MNE
         ELEM_LG(J) = 0
      end do
      !C
      RESNODE = 0
      NODES = 0
      do J = 1, NNODG
         !Csb-
         NODE_GL1(J) = PROC(J)
         if (IPROC == PROC(J)) then
            RESNODE = RESNODE + 1
            !Csb-
            !C           NODE_GL1(J) = PROC(J)
            NODE_GL2(J) = RESNODE
            NODE_LG(RESNODE) = J
         end if
      end do
      !C     DO I = 1, NNODG
      !C        print *, I, NODE_GL1(I)
      !C     ENDDO

      !C STEP 2:
      !C  Construct Local-to-Global Element Map
      !C  Add an element to the map if it has an resident node

      NELEM = 0
      do K = 1, NELG
         N1 = NNEG(1, K)
         N2 = NNEG(2, K)
         N3 = NNEG(3, K)
         PE1 = NODE_GL1(N1) ! Is any vertex a resident node?
         PE2 = NODE_GL1(N2)
         PE3 = NODE_GL1(N3)
         if ((PE1 == IPROC) .or. (PE2 == IPROC) .or. (PE3 == IPROC)) then
            NELEM = NELEM + 1
            ELEM_LG(NELEM) = K
         end if
      end do

      !C STEP 3:
      !C--Using Local-to-Global Element map
      !C  reconstruct Local-to-Global Node map
      !C  and Global-to-Local map for resident nodes
      !C
      ITOT = 0
      do J = 1, NELEM
         IEL = ELEM_LG(J)
         do L = 1, 3
            ITOT = ITOT + 1
            ITVECT(ITOT) = NNEG(L, IEL)
         end do
      end do
      ITEMP = ITOT
      if (ITOT > VTMAX) stop 'step3 decomp'
      call SORT(ITEMP, ITVECT) ! Sort and remove multiple occurrences
      ITOT = 1
      NODE_LG(1) = ITVECT(1)
      if (NODE_GL1(ITVECT(1)) == IPROC) then
         NODE_GL2(ITVECT(1)) = 1
      end if
      do J = 2, ITEMP
         if (ITVECT(J) /= NODE_LG(ITOT)) then
            ITOT = ITOT + 1
            NODE_LG(ITOT) = ITVECT(J)
            if (NODE_GL1(ITVECT(J)) == IPROC) then
               NODE_GL2(ITVECT(J)) = ITOT
            end if
         end if
      end do
      NODES = ITOT

      !C STEP 4:
      !C--If there are any global Weir-node pairs, construct
      !C  Local-to-Global Weir Node maps
      !C  Rule: if a global Weir node is assigned ( as resident or ghost node )
      !C        then make it and its dual a local Weir-node pair
      if (NWEIR > 0) then
         ITOT = 0
         ! Seizo 2008.07.14
         do J = 1, NWEIR
            INDEX = WEIR(J)
            INDEX2 = WEIRD(J)
            do K = 1, NODES
               N1 = NODE_LG(K)
               if ((N1 == INDEX) .or. (N1 == INDEX2)) then
                  ITOT = ITOT + 1
                  ITVECT(ITOT) = J
                  exit
               end if
            end do
         end do
         !
         ! Seizo DO J = 1,NODES
         ! Seizo    CALL SEARCH(WEIR,NWEIR,NODE_LG(J),INDEX)
         ! Seizo    IF (INDEX.NE.0) THEN
         ! Seizo      ITOT = ITOT+1
         ! Seizo      ITVECT(ITOT) = INDEX
         ! Seizo    ENDIF
         ! Seizo    CALL SEARCH(WEIRD,NWEIR,NODE_LG(J),INDEX2)
         ! Seizo    IF (INDEX2.NE.0) THEN
         ! Seizo      ITOT = ITOT+1
         ! Seizo      ITVECT(ITOT) = INDEX2
         ! Seizo    ENDIF
         ! Seizo ENDDO
         NLWEIR = 0
         ITEMP = ITOT
         if (ITOT > VTMAX) stop 'step4 decomp'
         if (ITEMP > 1) then
            call SORT(ITEMP, ITVECT)
            ITOT = 1
            INDEX = ITVECT(1)
            LWEIR_LG(ITOT) = WEIR(INDEX)
            LWEIRD_LG(ITOT) = WEIRD(INDEX)
            do J = 2, ITEMP
               if (ITVECT(J) /= INDEX) then
                  INDEX = ITVECT(J)
                  ITOT = ITOT + 1
                  LWEIR_LG(ITOT) = WEIR(INDEX)
                  LWEIRD_LG(ITOT) = WEIRD(INDEX)
               end if
            end do
            NLWEIR = ITOT
         end if
      else
         NLWEIR = 0
      end if
      !c     print *, "domsize: Number of WEIR node-pairs on PE",IPROC-1,
      !c    &         " = ",NLWEIR
      if (NLWEIR > NWEIR) then
         print *, "error in domsize: "
         print *, "local number of weir-pairs exceeds total"
         stop
      end if

      !C STEP 5:
      !C--If there are any global Weir-node pairs,
      !C  Re-construct Local-to-Global Element Map: IMAP_EL_LG(:,PE)
      !C  Rule:  Add an element if it has an resident node or
      !C         has the dual Weir node of a resident or ghost node

      ONELEM = NELEM ! Save NELEM for PEs with no WEIR-pairs
      NELEM = 0
      if (NLWEIR > 0) then
         do K = 1, NELG
            N1 = NNEG(1, K)
            N2 = NNEG(2, K)
            N3 = NNEG(3, K)
            PE1 = NODE_GL1(N1) ! Is any vertex a resident node?
            PE2 = NODE_GL1(N2)
            PE3 = NODE_GL1(N3) ! belong to a Weir-node pair ?
            call SEARCH3(LWEIR_LG(1), NLWEIR, N1, N2, N3, INDEX)
            call SEARCH3(LWEIRD_LG(1), NLWEIR, N1, N2, N3, INDEX2)
            if ((PE1 == IPROC) .or. (PE2 == IPROC) .or. (PE3 == IPROC) &
                .or. (INDEX /= 0) .or. (INDEX2 /= 0)) then
               NELEM = NELEM + 1
               ELEM_LG(NELEM) = K
            end if
         end do
         !C
      end if
      if (NELEM == 0) NELEM = ONELEM ! if necessary restore old nelem

      !Csb-DG2
      !C STEP 5-DG-3:
      !C--If there are any elements shared by three different sub-domains,
      !C  Re-construct Local-to-Global Element Map: IMAP_EL_LG(:,PE)
      !C  Rule: Add non-local elements adjacent to the elements
      !C        shared by three different sub-domains
      !C        to the local element list
      !C
      ONELEM = NELEM
      NEL_ADDED = 0
      do J = 1, ONELEM
         K = ELEM_LG(J)

         N1 = NNEG(1, K)
         N2 = NNEG(2, K)
         N3 = NNEG(3, K)

         PE1 = NODE_GL1(N1)
         PE2 = NODE_GL1(N2)
         PE3 = NODE_GL1(N3)

         if (K == 10352 .and. IPROC == 16) then
            print *, 'B=', IPROC, NELEM, K, PE1, PE2, PE3
         end if
         if (((PE1 /= PE2) .and. (PE2 /= PE3) .and. (PE3 /= PE1)) .and. &
             (PE1 == IPROC .or. PE2 == IPROC .or. PE3 == IPROC)) then
            !C        IF((PE1.NE.PE2).AND.(PE2.NE.PE3).AND.(PE3.NE.PE1)) THEN
            if (NELEM >= MNE) then
               stop 'FATAL ERROR (STEP 5-DG-3 DRY)'
            else
               NELEM = NELEM + 1
               ELEM_LG(NELEM) = K
               NEL_ADDED = NEL_ADDED + 1
            end if
         end if
      end do
      !C--

      !C STEP 6:
      !C Using Local-to-Global Element map, reconstruct Local-to-Global Node map
      !C
      ITOT = 0
      do J = 1, NELEM
         IEL = ELEM_LG(J)
         do L = 1, 3
            ITOT = ITOT + 1
            ITVECT(ITOT) = NNEG(L, IEL)
         end do
      end do
      ITEMP = ITOT
      if (ITOT > VTMAX) stop 'step6 decomp'
      call SORT(ITEMP, ITVECT) ! Sort and remove multiple occurrences
      ITOT = 1
      NODE_LG(1) = ITVECT(1)
      do J = 2, ITEMP
         if (ITVECT(J) /= NODE_LG(ITOT)) then
            ITOT = ITOT + 1
            NODE_LG(ITOT) = ITVECT(J)
         end if
      end do

      !Csb
      !C STEP 6-DG
      !C--Increase ITOT to take into acount NEL_ADDED
      !C
      ITOT = ITOT + NEL_ADDED*3
      !C--

      !C
      NODES = ITOT
      if (NODES > MNPP) MNPP = NODES
      !c     print *, "Number of nodes on PE",IPROC-1," = ",NODES
      if (NELEM > MNEP) MNEP = NELEM
      !c     print *, "Number of elements on PE",IPROC-1," = ",NELEM

      !1000 CONTINUE
   end do

   !C
   deallocate (ITVECT, NODE_LG, NODE_GL1, NODE_GL2, ELEM_LG, &
               LWEIR_LG, LWEIRD_LG)
   !C
   print *, " Setting MNPP = ", MNPP
   print *, " Setting MNEP = ", MNEP
   !C
   !RETURN
end subroutine DOMSIZE

subroutine SORT(N, RA)
   !IMPLICIT NONE
   integer :: N, L, IR, RRA, I, J
   integer :: RA(N)
   !C
   !C--------------------------------------------------------------------------
   !C  Sorts array RA of length N into ascending order using Heapsort algorithm.
   !C  N is input; RA is replaced on its output by its sorted rearrangement.
   !C  Ref: Numerical Recipes
   !C--------------------------------------------------------------------------
   !C
   L = N/2 + 1
   IR = N
10 continue
   if (L > 1) then
      L = L - 1
      RRA = RA(L)
   else
      RRA = RA(IR)
      RA(IR) = RA(1)
      IR = IR - 1
      if (IR == 1) then
         RA(1) = RRA
         return
      end if
   end if
   I = L
   J = L + L
20 if (J <= IR) then
      if (J < IR) then
         if (RA(J) < RA(J + 1)) J = J + 1
      end if
      if (RRA < RA(J)) then
         RA(I) = RA(J)
         I = J
         J = J + J
      else
         J = IR + 1
      end if
      GO TO 20
   end if
   RA(I) = RRA
   GO TO 10
end subroutine SORT

subroutine LOCATE(XX, N, X, J)
   !IMPLICIT NONE
   integer :: JM, JL, JU, J, N, X, XX(N)
   !C
   !C--Given an array XX of length N, and given a value X, returns a value J
   !C--such that X is between XX(J) and XX(J+1). XX must be monotonic, either
   !C--increasing or decreasing. J=0 or J=N is returned to indicate that X is
   !C--out of range.
   !C--
   !C--NUMERICAL RECIPES - The Art of Scientific Computing [FORTRAN Version]

   !C--Initialize lower and upper limits
   JL = 0
   JU = N + 1
   !C
   !C--If we are not done yet, compute a mid-point, and replace either the lower
   !C--limit or the upper limit, as appropriate.
   !C
10 if (JU - JL > 1) then
      JM = (JU + JL)/2
      if ((XX(N) > XX(1)) .eqv. (X > XX(JM))) then
         JL = JM
      else
         JU = JM
      end if
      !C--Repeat until the test condition 10 is satisfied.
      GO TO 10
   end if
   !C--Then set the output and return.
   J = JL
   return
end subroutine LOCATE

subroutine SEARCH3(MAP, LEN, N1, N2, N3, INDEX)
   integer :: MAP(*), LEN, N1, N2, N3, INDEX, IP, I
   !cvjp  rewritten 5/3/99
   INDEX = 0
   do I = 1, LEN
      IP = MAP(I)
      if (IP == N1 .or. IP == N2 .or. IP == N3) then
         INDEX = I
         goto 99
      end if
   end do
99 return
end subroutine SEARCH3

subroutine SEARCH(MAP, N, target, INDEX)
   integer :: MAP(*), N, target, INDEX, I
   !C
   INDEX = 0
   if (N == 0) goto 99
   do I = 1, N
      if (MAP(I) == target) then
         INDEX = I
         goto 99
      end if
   end do
99 return
end subroutine SEARCH
