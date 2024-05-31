!=========================================================
!     The purpose of this function is to generate the data structures
!     that will be used in the Configuration Bias Monte Carlo regrowth
!     functions.
!     The terminology used in this segment is as followed
!           Terminal: An atom that has one bond associated with it. These atoms are
!                     located at the end of a molecular chain.
!             Linker: An atom that has two bonds associated with it.  These atoms are
!                     segments found in the middle of a long chain
!                Hub: An atom that has three or more bonds assocaited with it. These atoms
!                     are where the molecule splits into multiple branches.
      subroutine CBMC_CreateTopology
         use Forcefield
         use SimParameters
         use Coords
         use CBMC_Variables
         use ParallelVar
         implicit none
         integer :: i, iBin, iType, iAtom, iBond
         integer :: AllocateStatus
         integer :: cnt, curAtom
         integer :: nTerminal, nLinker, nHub

         select case (trim(adjustl(ForceFieldName)))
         case ("Pedone")
            allocate (regrowType(1:nMolTypes), STAT=AllocateStatus)
            regrowType = 0
            return
         case default
            continue
         end select

!      Allocate the Topology arrays
         allocate (topolArray(1:nMolTypes), STAT=AllocateStatus)
         do iType = 1, nMolTypes
            allocate (topolArray(iType)%atom(1:nAtoms(iType)), STAT=AllocateStatus)
         end do
         allocate (regrowType(1:nMolTypes), STAT=AllocateStatus)
         allocate (regrowOrder(1:nMolTypes, 1:maxAtoms), STAT=AllocateStatus)
         allocate (pathArray(1:nMolTypes), STAT=AllocateStatus)
         allocate (usedByPath(1:nMolTypes, 1:maxAtoms), STAT=AllocateStatus)
         allocate (atomPathIndex(1:nMolTypes, 1:maxAtoms), STAT=AllocateStatus)
         allocate (probTypeCBMC(1:nMolTypes), STAT=AllocateStatus)

!      Search through each molecule and determine the number of bonds each atom is
!      a part of. This is done in order to classify each atom as either a Terminal atom(numBonds = 1),
!      Linker Atom(numBonds = 2), or Hub Atom (numBonds >= 3).  This information is used to classify
!      which type of CBMC algorithm will be required to regrow each molecule.
         totalDihed = 0
         do iType = 1, nMolTypes
            do iAtom = 1, nAtoms(iType)
               cnt = 0
               topolArray(iType)%atom(iAtom) = 0
               do iBond = 1, nBonds(iType)
                  if (bondArray(iType, iBond)%bondMembr(1) .eq. iAtom) then
                     cnt = cnt + 1
                     cycle
                  end if
                  if (bondArray(iType, iBond)%bondMembr(2) .eq. iAtom) then
                     cnt = cnt + 1
                     cycle
                  end if
               end do
               topolArray(iType)%atom(iAtom) = int(cnt, atomIntType)
            end do
         end do

!      Using the data that was previously calcualted classify each atom as terminal, linker, or hub and then tabulate
!      the data.
         do iType = 1, nMolTypes
            nTerminal = 0
            nHub = 0
            nLinker = 0
            do iAtom = 1, nAtoms(iType)
               if (topolArray(iType)%atom(iAtom) .le. 1) then
                  nTerminal = nTerminal + 1
               elseif (topolArray(iType)%atom(iAtom) .eq. 2) then
                  nLinker = nLinker + 1
               elseif (topolArray(iType)%atom(iAtom) .ge. 3) then
                  nHub = nHub + 1
!             Count the number of dihedral angles that are used.
                  if (topolArray(iType)%atom(iAtom) .eq. 3) then
                     totalDihed = totalDihed + 1
                  elseif (topolArray(iType)%atom(iAtom) .eq. 4) then
                     totalDihed = totalDihed + 2
                  end if
               end if
            end do
            pathArray(iType)%nTerminal = nTerminal
            pathArray(iType)%nLinker = nLinker
            pathArray(iType)%nHub = nHub
            if (nTerminal .ne. 0) then
               allocate (pathArray(iType)%termAtoms(1:nTerminal), STAT=AllocateStatus)
            else
               write (nout, *) "Molecule is Cyclic, not yet supported"
               stop "Molecule is Cyclic, not yet supported"
            end if
            if (nHub .ne. 0) then
               allocate (pathArray(iType)%hubAtoms(1:nHub), STAT=AllocateStatus)
               if (nHub .gt. 1) then
                  stop "Molecule contains Multiple Hubs, not yet supported"
               end if
            end if

!        This block stores the atom indicies for each Terminal and Hub atom for quick access
!        later in the program.
            nTerminal = 0
            nHub = 0
            nLinker = 0
            do iAtom = 1, nAtoms(iType)
               if (topolArray(iType)%atom(iAtom) .le. 1) then
                  nTerminal = nTerminal + 1
                  pathArray(iType)%termAtoms(nTerminal) = iAtom
               elseif (topolArray(iType)%atom(iAtom) .ge. 3) then
                  nHub = nHub + 1
                  pathArray(iType)%hubAtoms(nHub) = iAtom
               end if
            end do
         end do

         allocate (dihedData(1:totalDihed), STAT=AllocateStatus)
         cnt = 0
         write (35, *) "--------------------------------------------"
         write (35, *) "Dihedral Angles"
         do iType = 1, nMolTypes
            do iAtom = 1, nAtoms(iType)
               if (topolArray(iType)%atom(iAtom) .eq. 3) then
                  cnt = cnt + 1
                  dihedData(cnt)%molType = iType
                  dihedData(cnt)%hubIndx = iAtom
                  dihedData(cnt)%dihedIndx = 1
                  write (35, *) iType, iAtom, 1
               elseif (topolArray(iType)%atom(iAtom) .eq. 4) then
                  cnt = cnt + 1
                  dihedData(cnt)%molType = iType
                  dihedData(cnt)%hubIndx = iAtom
                  dihedData(cnt)%dihedIndx = 1
                  write (35, *) iType, iAtom, 1
                  cnt = cnt + 1
                  dihedData(cnt)%molType = iType
                  dihedData(cnt)%hubIndx = iAtom
                  dihedData(cnt)%dihedIndx = 2
                  write (35, *) iType, iAtom, 2
               end if
            end do
         end do

         do i = 1, totalDihed
            dihedData(i)%Prob = 1E0/nDihBins
            dihedData(i)%Hist = 0E0
            dihedData(i)%Integral(0) = dihedData(i)%Prob(0)
            do iBin = 1, nDihBins
               dihedData(i)%Integral(iBin) = dihedData(i)%Integral(iBin - 1) + dihedData(i)%Prob(iBin)
            end do
            dihedData(i)%accConst = nDihBins*diBinSize
         end do

         call CBMC_FindPathways

         usedByPath = 0
         do iType = 1, nMolTypes
            do i = 1, pathArray(iType)%nPaths
               do iAtom = 1, nAtoms(iType)
                  curAtom = pathArray(iType)%path(i, iAtom)
                  if (curAtom .eq. 0) then
                     cycle
                  end if
                  if (topolArray(iType)%atom(curAtom) .le. 2) then
                     usedByPath(iType, curAtom) = i
                     atomPathIndex(iType, curAtom) = iAtom
                  else
                     usedByPath(iType, curAtom) = -1
                     atomPathIndex(iType, curAtom) = -1
                  end if
               end do
            end do
         end do

         call DetermineRegrowthType
         call DetermineRegrowOrder

      end subroutine
!=========================================================
!      The purpose of this subroutine is to divide each molecule into chain segments.  A segment is the sequence
!      of atoms that start with either a hub or a terminal atom and end with another terminal or hub atom.
!      For instance in the molecule 2-methyl-butane [CH3-CH2-CH(CH3)-CH3] one segment would
!      be defined as CH3-CH2-CH since CH3 is a terminal atom (numBonds = 1) and CH is a Hub atom
!      (numBonds >= 3).  This information is used by the CBMC algorithm to map out which algorithms
!      need to be used in order to efficiently regrow each atom.  For instance a Hub atom is much more
!      difficult to regrow than a Linker (numBonds = 2) since there are multiple bending angles associated
!      with this regrowth.
      subroutine CBMC_FindPathways
         use Forcefield
         use SimParameters
         use Coords
         use CBMC_Variables
         implicit none
         logical :: atomUsed(1:maxAtoms)
         integer :: i, cnt
         integer :: iType, iAtom, iPath
         integer :: nPaths, curBonds
         integer, allocatable :: tempArray(:, :)
         integer :: atmPrev, atmCur, atmNext
         integer :: AllocateStatus

         allocate (tempArray(1:50, 1:maxAtoms))

         do iType = 1, nMolTypes
!         This first segment finds all the segments that begin at a terminal atom
            tempArray = 0
            atomUsed = .false.
            pathArray(iType)%nPaths = 0
            do iAtom = nAtoms(iType) + 1, maxAtoms
               atomUsed(iAtom) = .true.
            end do

            if (nAtoms(iType) .eq. 1) then
               pathArray(iType)%nPaths = 1
               nPaths = pathArray(iType)%nPaths
               tempArray(1, 1:1) = 1
               goto 56
            end if
            do i = 1, pathArray(iType)%nTerminal
               atmPrev = 0
               atmCur = pathArray(iType)%termAtoms(i)
               tempArray(i, 1) = atmCur
               atomUsed(atmCur) = .true.
               curBonds = 2
               cnt = 1
               do while (curBonds .eq. 2)
                  call getNextAtom(iType, atmPrev, atmCur, atmNext)
                  atmPrev = atmCur
                  atmCur = atmNext
                  cnt = cnt + 1
                  tempArray(i, cnt) = atmNext
                  atomUsed(atmNext) = .true.
                  curBonds = topolArray(iType)%atom(atmNext)
               end do
!            If curBonds = 1 we have gone from a Terminal atom to another Terminal atom
!            implying we have a perfectly straight chain molecule.  Thus only one segment is
!            present.  If curBonds >=3 we have reached a hub atom which ends the current
!            path and the hub will now.
               if (curBonds .ne. 2) then
                  pathArray(iType)%nPaths = pathArray(iType)%nPaths + 1
                  if (all(atomUsed .eqv. .true.)) then
                     goto 56
                  end if
               end if
            end do
!         If all atoms have been accounted for then no further calculations are required

!        --------------------------------------------------------------------
!            This section will create pathways between two hubs.  Not yet completed.
            stop "Hub to Hub Regrowth not currently supported"
!        --------------------------------------------------------------------
56          continue

            nPaths = pathArray(iType)%nPaths
            allocate (pathArray(iType)%path(1:nPaths, 1:nAtoms(iType)), STAT=AllocateStatus)
            pathArray(iType)%path = 0

            do i = 1, nPaths
               pathArray(iType)%path(i, 1:nAtoms(iType)) = tempArray(i, 1:nAtoms(iType))

            end do
         end do

         do iType = 1, nMolTypes
            nPaths = pathArray(iType)%nPaths
            allocate (pathArray(iType)%pathMax(1:nPaths), STAT=AllocateStatus)

            do iPath = 1, nPaths
               pathArray(iType)%pathMax(iPath) = nAtoms(iType)
               do iAtom = 1, nAtoms(iType)
                  if (pathArray(iType)%path(iPath, iAtom) .eq. 0) then
                     pathArray(iType)%pathMax(iPath) = iAtom - 1
                     exit
                  end if
               end do
            end do
         end do

         write (35, *) "------------------------------------------------------------------"
         do iType = 1, nMolTypes
            write (35, *) "iType:", iType
            write (35, *) "Number of Paths:", pathArray(iType)%nPaths
            do iPath = 1, pathArray(iType)%nPaths
               write (35, *) pathArray(iType)%path(iPath, 1:pathArray(iType)%pathMax(iPath))
            end do
            write (35, *)
         end do
         write (35, *) "------------------------------------------------------------------"

      end subroutine
!=========================================================
      subroutine getNextAtom(iType, atmPrev, atmCur, atmNext)
         use Forcefield
         use SimParameters
         use Coords
         implicit none

         integer, intent(in) :: iType, atmPrev, atmCur
         integer, intent(out) :: atmNext

         logical :: atmFound
         integer :: iBond

         do iBond = 1, nBonds(iType)
            atmFound = .false.
            if (any(bondArray(iType, iBond)%bondMembr .eq. atmCur)) then
               atmFound = .true.
            end if
            if (atmFound) then
               if (all(bondArray(iType, iBond)%bondMembr .ne. atmPrev)) then
                  if (bondArray(iType, iBond)%bondMembr(1) .eq. atmCur) then
                     atmNext = bondArray(iType, iBond)%bondMembr(2)
!                bondUsed = iBond
                     return
                  else
                     atmNext = bondArray(iType, iBond)%bondMembr(1)
!                bondUsed = iBond
                     return
                  end if
               end if
            end if
         end do

         write (*, *) "Error in getNextAtom subroutine!!"
         write (*, *) "Dead End! Unable to find next atom in the chain"
         write (*, *) "Mol Type:", iType
         write (*, *) "Previous Atom:", atmPrev
         write (*, *) "Current Atom:", atmCur
         stop

      end subroutine
!===============================================================================
!     This subroutine determines the type of regrowth function to be used in
!     a CBMC or AVBMC regrowth.  Each molecule is assigned an integer based on
!     the molecule's flexbility and size.  The integers correspond to these molecule types
!           0: Completely Ridgid Molecule
!           1: Simple Flexible Molecule
!           2: Straight Chain Molecule
!           3: Branched Chain Molecule
!
      subroutine DetermineRegrowthType
         use Forcefield
         use SimParameters
         use Coords
         use CBMC_Variables
         implicit none

         logical :: bondRidgid, bendRidgid
         integer :: iType, i
         integer :: xType
         real(dp) :: xForceConst, norm

         regrowType = -1

         do iType = 1, nMolTypes
            if (nAtoms(iType) .eq. 1) then
               regrowType(iType) = 0
               cycle
            end if
            bondRidgid = .true.
            do i = 1, nBonds(iType)
               xType = bondArray(iType, i)%bondType
               xForceConst = bondData(xType)%k_eq
               if (xForceConst .ne. 0d0) then
                  bondRidgid = .false.
                  exit
               end if
            end do
            bendRidgid = .true.
            do i = 1, nAngles(iType)
               xType = bendArray(iType, i)%bendType
               xForceConst = bendData(xType)%k_eq
               if (xForceConst .ne. 0d0) then
                  bendRidgid = .false.
                  exit
               end if
            end do
            if ((bondRidgid .eqv. .true.) .and. (bendRidgid .eqv. .true.)) then
               if (nTorsional(iType) .eq. 0) then
                  regrowType(iType) = 0
               end if
            else
               if (nTorsional(iType) .eq. 0) then
                  regrowType(iType) = 1
               else
                  if (pathArray(iType)%nHub .eq. 0) then
                     regrowType(iType) = 2
                  else
                     regrowType(iType) = 3
                  end if
               end if
            end if
            write (35, *) "Regrow Type:", regrowType(iType)
         end do

         probTypeCBMC = 0d0
         do iType = 1, nMolTypes
            if (regrowType(iType) .ge. 1) then
               probTypeCBMC(iType) = 1d0
            end if
         end do
         norm = sum(probTypeCBMC)
         if (any(probTypeCBMC .ne. 0d0)) then
            do iType = 1, nMolTypes
               probTypeCBMC(iType) = probTypeCBMC(iType)/norm
            end do
         end if

      end subroutine
!===============================================================================
!
      subroutine DetermineRegrowOrder
         use Forcefield
         use SimParameters
         use Coords
         use CBMC_Variables
         implicit none

         logical :: bondRidgid, bendRidgid
         integer :: iType, iAtom, Atom1Loc, nextAtom
         integer :: cnt

         regrowOrder = 0

         do iType = 1, nMolTypes
            select case (regrowType(iType))
            case (0)
               continue

            case (1)
               continue

            case (2)
               regrowOrder(iType, 1) = 1
               do iAtom = 1, pathArray(iType)%pathMax(1)
                  if (pathArray(iType)%path(1, iAtom) .eq. 1) then
                     Atom1Loc = iAtom
                     exit
                  end if
               end do

               cnt = 1
               if (dble(Atom1Loc) .ge. dble(pathArray(iType)%pathMax(1) + 1)/2d0) then
                  do iAtom = Atom1Loc + 1, pathArray(iType)%pathMax(1)
                     cnt = cnt + 1
                     regrowOrder(iType, cnt) = pathArray(iType)%path(1, iAtom)
                  end do
                  do iAtom = 1, Atom1Loc - 1
                     cnt = cnt + 1
                     regrowOrder(iType, cnt) = pathArray(iType)%path(1, iAtom)
                  end do
               else
                  do iAtom = 1, Atom1Loc - 1
                     cnt = cnt + 1
                     regrowOrder(iType, cnt) = pathArray(iType)%path(1, iAtom)
                  end do
                  do iAtom = Atom1Loc + 1, pathArray(iType)%pathMax(1)
                     cnt = cnt + 1
                     regrowOrder(iType, cnt) = pathArray(iType)%path(1, iAtom)
                  end do
               end if
            case default
               write (*, *) iType, regrowType(iType)
               write (*, *) "Error in DetermineRegrowOrder function!"
               write (*, *) "Invalid Regrow Type!"
               stop

            end select

         end do

      end subroutine
!=========================================================

