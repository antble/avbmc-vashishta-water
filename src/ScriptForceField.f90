!========================================================
module ForceFieldInput
   use VarPrecision
   use Units

   integer, parameter :: maxLineLen = 500
   abstract interface
      subroutine CommonSub
         use VarPrecision
         implicit none
      end subroutine
   end interface

   abstract interface
      subroutine interSub(lineStore)
         implicit none
         character(len=*), intent(in) :: lineStore(:)
      end subroutine
   end interface

   private
   logical :: fieldTypeSet = .false.

   procedure(CommonSub), pointer :: commonFunction => NULL()
   procedure(CommonSub), pointer :: FFSpecificFlags => NULL()

   procedure(interSub), pointer :: interFunction => NULL()
   procedure(interSub), pointer :: interFunction2 => NULL()
   procedure(interSub), pointer :: interFunction3 => NULL()
   
   real(dp) :: convEng = 1d0
   real(dp) :: convDist = 1d0
   real(dp) :: convAng = 1d0

   real(dp), parameter :: coulombConst = 1.671009770E5_dp

   public :: SetForcefieldType, ScriptForcefield, fieldTypeSet

!========================================================
contains
!========================================================
   subroutine SetForcefieldType(potenType)
      use ForceField, only: ForceFieldName
      use EnergyPointers
      use SwapBoundary
      use ParallelVar, only: nout
!      use UmbrellaSamplingNew, only: ReadInput_Umbrella
      use WHAM_Functions
      implicit none
      character(len=10) :: potenType
      character(len=10) :: lenDefault, distDefault, angDefaults

      lenDefault = "kb"
      distDefault = "ang"
      angDefaults = "deg"

      convEng = FindEngUnit(lenDefault)
      convDist = FindLengthUnit(distDefault)
      convAng = FindAngularUnit(angDefaults)

      select case (trim(adjustl(potenType)))
      case ("lj_q")
         write (nout, *) "Forcefield Type: Standard Lennard-Jones w/ Eletrostatic"
         ForceFieldName = "LJ_Q"
         Detailed_ECalc => Detailed_EnergyCalc_LJ_Q
         Shift_ECalc => Shift_EnergyCalc_LJ_Q
         SwapIn_ECalc => SwapIn_EnergyCalc_LJ_Q
         SwapOut_ECalc => SwapOut_EnergyCalc_LJ_Q
         Rosen_Mol_New => Rosen_BoltzWeight_Molecule_New
         Rosen_Mol_Old => Rosen_BoltzWeight_Molecule_Old
         Quick_Nei_ECalc => QuickNei_ECalc_Inter_LJ_Q
         boundaryFunction => Bound_MaxMin
         commonFunction => Allocate_LJ_Q
         FFSpecificFlags => LJ_SetFlags
         interFunction => Read_LJ_Q
         fieldTypeSet = .true.
      case ("pedone")
         write (nout, *) "Forcefield Type: Pedone"
         ForceFieldName = "Pedone"
         Detailed_ECalc => Detailed_EnergyCalc_Pedone
         Shift_ECalc => Shift_EnergyCalc_Pedone
         SwapIn_ECalc => SwapIn_EnergyCalc_Pedone
         SwapOut_ECalc => SwapOut_EnergyCalc_Pedone
         Rosen_Mol_New => Rosen_BoltzWeight_Pedone_New
         Rosen_Mol_Old => Rosen_BoltzWeight_Pedone_Old
         Quick_Nei_ECalc => QuickNei_ECalc_Inter_Pedone
         boundaryFunction => Bound_PedoneChargeBalance
         commonFunction => Allocate_Pedone
         FFSpecificFlags => Pedone_SetFlags
         interFunction => Read_Pedone
         fieldTypeSet = .true.
      case ("tersoff")
         write (nout, *) "Forcefield Type: Tersoff"
         ForceFieldName = "Tersoff"
         Detailed_ECalc => Detailed_EnergyCalc_Tersoff
         Shift_ECalc => Shift_EnergyCalc_Tersoff
         SwapIn_ECalc => SwapIn_EnergyCalc_Tersoff
         SwapOut_ECalc => SwapOut_EnergyCalc_Tersoff
         Rosen_Mol_New => Rosen_Tersoff_Molecule_New
         Rosen_Mol_Old => Rosen_Tersoff_Molecule_Old
         Quick_Nei_ECalc => QuickNei_ECalc_Inter_Tersoff
         boundaryFunction => Bound_MaxMin
         commonFunction => Allocate_Tersoff
         FFSpecificFlags => Tersoff_SetFlags
         interFunction => Read_Tersoff
         fieldTypeSet = .true.
      case ("vashishta")
         write (nout, *) "Forcefield Type: Vashishta"
         ForceFieldName = "Vashishta"
         Detailed_ECalc => Detailed_EnergyCalc_Vashishta
         Shift_ECalc => Shift_EnergyCalc_Vashishta
         SwapIn_ECalc => SwapIn_EnergyCalc_Vashishta
         SwapOut_ECalc => SwapOut_EnergyCalc_Vashishta
         Rosen_Mol_New => Rosen_Vashishta_Molecule_New
         ! Rosen_Mol_Old => Rosen_Vashishta_Molecule_Old
         Rosen_Mol_Old_V => Rosen_Vashishta_Molecule_Old
         ! Quick_Nei_ECalc => QuickNei_ECalc_Inter_Vashishta
         boundaryFunction => Bound_MaxMin
         commonFunction => Allocate_Vashishta
         FFSpecificFlags => Vashishta_SetFlags
         interFunction => Read_Vashishta
         interFunction2 => Read_Vashishta2
         interFunction3 => Read_Vashishta3
         fieldTypeSet = .true.
      case ("custompairwise")
      case default
         stop "Unknown potential type given in forcefield input"
      end select

   end subroutine

!========================================================
   subroutine ScriptForcefield(lineStore)
      use CBMC_Variables
      use Coords
      use CoordinateFunctions
      use EnergyTables
      use ForceField
      use SimParameters
      use Units
      use VarPrecision
      implicit none
      character(len=maxLineLen), intent(in) :: lineStore(:)

      character(len=25) :: dummy, command, command2, stringValue
      logical :: logicValue
      integer :: i, j, iLine, nLines, lineBuffer, lineStat
      integer :: intValue, AllocateStat
      integer :: nAtomsMax, nBondsMax, nAnglesMax
      integer :: nTorsMax, nImpropMax, nNonBondMax
      real(dp) :: realValue

      if (.not. fieldTypeSet) then
         write (*, *) "ERROR! Forcefield input has been called before the forcefield type has been set!"
         stop
      end if

      if (nMolTypes .eq. 0) then
         write (*, *) "ERROR! Forcefield input has been called before the number of molecule types have been defined"
         stop
      end if

      nAtomTypes = 0
      nBondTypes = 0
      nAngleTypes = 0
      nTorsionalTypes = 0

      call FindMax(lineStore, "atoms", nAtomsMax)
      call FindMax(lineStore, "nonbonded", nNonBondMax)
      call FindMax(lineStore, "bonds", nBondsMax)
      call FindMax(lineStore, "angles", nAnglesMax)
      call FindMax(lineStore, "torsional", nTorsMax)
      ALLOCATE (atomArray(1:nMolTypes, 1:nAtomsMax), STAT=AllocateStat)
      ALLOCATE (nonBondArray(1:nMolTypes, 1:nNonBondMax), STAT=AllocateStat)
      ALLOCATE (bondArray(1:nMolTypes, 1:nBondsMax), STAT=AllocateStat)
      ALLOCATE (bendArray(1:nMolTypes, 1:nAnglesMax), STAT=AllocateStat)
      ALLOCATE (torsArray(1:nMolTypes, 1:nTorsMax), STAT=AllocateStat)

      do i = 1, maxRosenTrial
         allocate (rosenTrial(i)%x(1:nAtomsMax))
         allocate (rosenTrial(i)%y(1:nAtomsMax))
         allocate (rosenTrial(i)%z(1:nAtomsMax))
      end do

      lineBuffer = 0
      nLines = size(lineStore)
      do iLine = 1, nLines
         if (lineBuffer .gt. 0) then
            lineBuffer = lineBuffer - 1
            cycle
         end if
         call GetCommand(lineStore(iLine), command, lineStat)
         if (lineStat .gt. 0) then
            cycle
         end if
         call LowerCaseLine(command)
!        write(*,*) lineStore(iLine)
         select case (adjustl(trim(command)))
         case ("define")
            call FindCommandBlock(iLine, lineStore, "end_define", lineBuffer)
            call DefineForcefield(lineStore(iLine:iLine + lineBuffer))
         case ("create")
            call FindCommandBlock(iLine, lineStore, "end_create", lineBuffer)
            call CreateForcefield(lineStore(iLine:iLine + lineBuffer))
         case ("set")
            call SetForcefieldParam(lineStore(iLine))
         case default
            write (*, *) "ERROR! Invalid command in Forcefield file on line:", iLine
            write (*, *) lineStore(iLine)
            stop "INPUT ERROR!"
         end select

      end do

      vmdAtoms = 0
      do i = 1, nMolTypes
         vmdAtoms = vmdAtoms + (NMAX(i)*nAtoms(i))
      end do

      totalMass = 0d0
      do i = 1, nMolTypes
         do j = 1, nAtoms(i)
            totalMass(i) = totalMass(i) + atomData(atomArray(i, j))%mass
         end do
      end do

      call AllocateCoordinateArrays
      call CreateJointArray
      call FFSpecificFlags

   end subroutine
!========================================================
   subroutine DefineForcefield(lineStore)
      use Coords
      use ForceField
      use SimParameters
      use Units
      use VarPrecision
      implicit none
      character(len=maxLineLen), intent(in) :: lineStore(:)
      character(len=30) :: dummy, defType, stringValue
      logical :: logicValue
      integer :: i, j, iLine, nLines, nParam
      integer :: intValue, AllocateStat
      real(dp) :: realValue

      nLines = size(lineStore)
      read (lineStore(1), *) dummy, defType
      call LowerCaseLine(defType)
      ! print *,"DEFINE: " , trim(adjustl(defType))
      
      select case (trim(adjustl(defType)))
      case ("atomtypes")
         if (allocated(atomData)) then
            write (*, *) "ERROR! AtomTypes already defined in the forcefield file"
            stop
         end if
         read (lineStore(1), *) dummy, defType, intValue
         nAtomTypes = intValue
         ALLOCATE (atomData(1:nAtomTypes), STAT=AllocateStat)
         call commonFunction
         call interFunction(lineStore)
      case ("paramtypes_2body")
         call interFunction2(lineStore)
      case ("paramtypes_3body")
         if (allocated(atomDataV)) then
            write (*, *) "ERROR! paramtypes_3body already defined in the forcefield file"
            stop
         end if
         read (lineStore(1), *) dummy, defType, intValue
         nAtomTypes2 = intValue
         ALLOCATE (atomDataV(1:nAtomTypes2), STAT=AllocateStat)
         call interFunction3(lineStore)
      case ("bondtypes")
         read (lineStore(1), *) dummy, defType, intValue
         nBondTypes = intValue
         ALLOCATE (bondData(1:nBondTypes), STAT=AllocateStat)
         i = 0
         do iLine = 2, nLines - 1
            i = i + 1
            read (lineStore(iLine), *) bondData(i)%bondName, bondData(i)%k_eq, bondData(i)%r_eq
            bondData(i)%k_eq = bondData(i)%k_eq*convEng
            bondData(i)%r_eq = bondData(i)%r_eq*convDist
            if (echoInput) then
               write (35, *) bondData(i)%bondName, bondData(i)%k_eq, bondData(i)%r_eq
            end if
         end do
      case ("angletypes")
         read (lineStore(1), *) dummy, defType, intValue
         nAngleTypes = intValue
         ALLOCATE (bendData(1:nAngleTypes), STAT=AllocateStat)
         i = 0
         do iLine = 2, nLines - 1
            i = i + 1
            read (lineStore(iLine), *) bendData(i)%angleName, bendData(i)%k_eq, bendData(i)%ang_eq
            if (echoInput) then
               write (35, *) bendData(i)%angleName, bendData(i)%k_eq, bendData(i)%ang_eq
            end if
            bendData(i)%k_eq = bendData(i)%k_eq*convEng
            bendData(i)%ang_eq = bendData(i)%ang_eq*convAng
         end do
      case ("torsiontypes")
         read (lineStore(1), *) dummy, defType, intValue
         nTorsionalTypes = intValue
         ALLOCATE (torsData(1:nTorsionalTypes), STAT=AllocateStat)
         i = 0
         do iLine = 2, nLines - 1
            i = i + 1
            read (lineStore(iLine), *) torsData(i)%torsName, nParam
            allocate (torsData(i)%a(1:nParam), STAT=AllocateStat)
            read (lineStore(iLine), *) torsData(i)%torsName, torsData(i)%nPara, (torsData(i)%a(j), j=1, nParam)
            if (echoInput) then
               write (35, *) torsData(i)%torsName, torsData(i)%nPara, (torsData(i)%a(j), j=1, nParam)
            end if
            do j = 1, nParam
               torsData(i)%a(j) = torsData(i)%a(j)*convEng
            end do
         end do
      case ("rmin")
         if (nAtomTypes .eq. 0) then
            write (*, *) "ERROR! RMin is called before the number of atom types"
            write (*, *) "have been defined"
            stop
         end if
         call SetRMin(lineStore)
      case default
         write (*, *) "ERROR! Invalid command in Forcefield file on line:"
         write (*, *) lineStore(1)
         stop "INPUT ERROR!"
      end select

   end subroutine
!========================================================
   subroutine SetForcefieldParam(line)
      use Coords
      use ForceField
      use SimParameters
      use Units
      use VarPrecision
      implicit none
      character(len=maxLineLen), intent(in) :: line
      character(len=30) :: dummy, deftype, stringValue
      character(len=25) :: command

      read (line, *) dummy, defType
      call LowerCaseLine(defType)
      select case (adjustl(trim(defType)))
      case ("lenunits")
         read (line, *) dummy, defType, stringValue
         convDist = FindLengthUnit(stringValue)
      case ("engunits")
         read (line, *) dummy, defType, stringValue
         convEng = FindEngUnit(stringValue)
      case ("angunits")
         read (line, *) dummy, defType, stringValue
         convAng = FindAngularUnit(stringValue)
      end select

   end subroutine
!========================================================
   subroutine CreateForcefield(lineStore)
      use Coords
      use ForceField
      use SimParameters
      use Units
      use VarPrecision
      implicit none
      character(len=*), intent(in) :: lineStore(:)
      character(len=30) :: dummy, labelField, stringValue
      character(len=25) :: command
      logical :: logicValue
      integer :: i, j, iUnit, iLine, jLine, nLines, nMol
      integer :: intValue, AllocateStat, lineStat, lineBuffer
      real(dp) :: realValue

      nLines = size(lineStore)
      read (lineStore(1), *) dummy, nMol

      do iLine = 2, nLines - 1
         call GetCommand(lineStore(iLine), command, lineStat)
         if (lineStat .gt. 0) then
            cycle
         end if
!        write(*,*) iLine, lineStore(iLine)
         call LowerCaseLine(command)
         select case (adjustl(trim(command)))
         case ("atoms")
            call FindCommandBlock(iLine, lineStore, "end_atoms", lineBuffer)
            read (lineStore(iLine), *) dummy, intValue
            if (lineBuffer - 1 .ne. intValue) then
               write (*, *) "Error in Forcefield Def! The Number of Atoms specified is different than the number of lines given."
               write (*, *) "Number of atoms specified:", intValue
               write (*, *) "Number of lines in command block:", lineBuffer - 1
               stop
            end if
            nAtoms(nMol) = intValue
            iUnit = 0
            ! atomArray contains the atomic composition of molecule #nMol
            do jLine = iLine + 1, iLine + lineBuffer - 1
               iUnit = iUnit + 1
               read (lineStore(jLine), *) atomArray(nMol, iUnit)
               ! print *, "atomArray", atomArray(nMol, iUnit), iUnit, jLine, nMol
            end do
         case ("bonds")
            call FindCommandBlock(iLine, lineStore, "end_bonds", lineBuffer)
            read (lineStore(iLine), *) dummy, intValue
            nBonds(nMol) = intValue
            iUnit = 0
            do jLine = iLine + 1, iLine + lineBuffer - 1
               iUnit = iUnit + 1
               read (lineStore(jLine), *) bondArray(nMol, iUnit)%bondType, (bondArray(nMol, iUnit)%bondMembr(j), j=1, 2)
            end do
         case ("angles")
            call FindCommandBlock(iLine, lineStore, "end_angles", lineBuffer)
            read (lineStore(iLine), *) dummy, intValue
            nAngles(nMol) = intValue
            iUnit = 0
            do jLine = iLine + 1, iLine + lineBuffer - 1
               iUnit = iUnit + 1
               read (lineStore(jLine), *) bendArray(nMol, iUnit)%bendType, (bendArray(nMol, iUnit)%bendMembr(j), j=1, 3)
            end do
         case ("torsional")
            call FindCommandBlock(iLine, lineStore, "end_torsional", lineBuffer)
            read (lineStore(iLine), *) dummy, intValue
            nTorsional(nMol) = intValue
            iUnit = 0
            do jLine = iLine + 1, iLine + lineBuffer - 1
               iUnit = iUnit + 1
               read (lineStore(jLine), *) torsArray(nMol, iUnit)%torsType, (torsArray(nMol, iUnit)%torsMembr(j), j=1, 4)
            end do
         case ("nonbonded")
            call FindCommandBlock(iLine, lineStore, "end_nonbonded", lineBuffer)
            read (lineStore(iLine), *) dummy, intValue
            nIntraNonBond(nMol) = intValue
            iUnit = 0
            do jLine = iLine + 1, iLine + lineBuffer - 1
               iUnit = iUnit + 1
               read (lineStore(jLine), *) (nonBondArray(nMol, iUnit)%nonMembr(j), j=1, 2)
            end do
         end select
      end do

   end subroutine
!========================================================
   subroutine FindCommandBlock(iLine, lineStore, endCommand, lineBuffer)
      implicit none
      integer, intent(in) :: iLine
      character(len=*), intent(in) :: lineStore(:)
      character(len=*), intent(in) :: endCommand
      integer, intent(out) :: lineBuffer
      logical :: found
      integer :: i, lineStat, nLines
      character(len=35) :: dummy

      dummy = " "
      nLines = size(lineStore)
      found = .false.
      do i = iLine + 1, nLines
         call GetCommand(lineStore(i), dummy, lineStat)
         if (trim(adjustl(dummy)) .eq. trim(adjustl(endCommand))) then
            lineBuffer = i - iLine
            found = .true.
            exit
         end if
      end do

      if (.not. found) then
         write (*, *) "ERROR! A command block was opened in the input script, but no closing END statement found!"
         write (*, *) lineStore(iLine)
         stop
      end if

   end subroutine
!========================================================
!     This subrotuine searches a given input line for the first command.
   subroutine SetRMin(lineStore)
      use VarPrecision
      use ForceField
      use ForceFieldFunctions
      use SimParameters
      implicit none
      character(len=maxLineLen), intent(in) :: lineStore(:)

      logical :: custom
      integer :: i, j, iLine, nLines
      integer :: indx1, indx2
      character(len=30) :: dummy, dummy2, mixingRule
      procedure(MixRule), pointer :: rmin_func => null()
      real(dp) :: curVal

      nLines = size(lineStore)
      read (lineStore(1), *) dummy, dummy2, mixingRule

      custom = .false.
      select case (trim(adjustl(mixingRule)))
      case ("geometric")
         rmin_func => GeoMean_MixingFunc
      case ("average")
         rmin_func => Mean_MixingFunc
      case ("min")
         rmin_func => Min_MixingFunc
      case ("max")
         rmin_func => Max_MixingFunc
      case ("custom")
         custom = .true.
      case default
         stop "Error! Invalid R_Min Mixing Rule Type"
      end select

      r_min = 0E0_dp
      r_min_tab = 0E0_dp
      if (custom) then
         do iLine = 2, nLines - 1
            read (lineStore(iLine), *) indx1, indx2, curVal
            if (indx1 .gt. nAtomTypes) then
               cycle
            end if
            if (indx2 .gt. nAtomTypes) then
               cycle
            end if
            if (echoInput) then
               write (35, *) indx1, indx2, curVal
               flush (35)
            end if
            r_min_tab(indx1, indx2) = curVal
            r_min_tab(indx2, indx1) = curVal
         end do
         do i = 1, nAtomTypes
            do j = 1, nAtomTypes
               r_min_tab(i, j) = r_min_tab(i, j)**2
            end do
         end do

      else

         do iLine = 2, nLines - 1
            read (lineStore(iLine), *) indx1, curVal
            r_min(indx1) = curVal
         end do
         do i = 1, nAtomTypes
            do j = i, nAtomTypes
               r_min_tab(i, j) = rmin_func(r_min(i), r_min(j))**2
               r_min_tab(j, i) = r_min_tab(i, j)
            end do
         end do

      end if

      write (35, *) "Rmin Table:"
      do i = 1, nAtomTypes
         write (35, *) (sqrt(r_min_tab(i, j)), j=1, nAtomTypes)
      end do

   end subroutine
!========================================================
!     This subrotuine searches a given input line for the first command.
   subroutine GetCommand(line, command, lineStat)
      use VarPrecision
      implicit none
      character(len=*), intent(in) :: line
      character(len=25), intent(out) :: command
      integer, intent(out) :: lineStat
      integer :: i, sizeLine, lowerLim, upperLim

      sizeLine = len(line)
      lineStat = 0
      i = 1
      command = " "
!      Find the first non-blank character in the string
      do while (i .le. sizeLine)
         if (ichar(line(i:i)) .ne. ichar(' ')) then
            !If a non-blank character is found, check first to see if it is the comment character.
            if (ichar(line(i:i)) .eq. ichar('#')) then
               lineStat = 1
               return
            else
               exit
            end if
         end if
         i = i + 1
      end do
!      If no characters are found the line is empty,
      if (i .ge. sizeLine) then
         lineStat = 1
         return
      end if
      lowerLim = i

!
      do while (i .le. sizeLine)
         if (line(i:i) .eq. " ") then
            exit
         end if
         i = i + 1
      end do
      upperLim = i

      command = line(lowerLim:upperLim)

   end subroutine
!================================================================
   subroutine Allocate_Vashishta
      use SimParameters
      use ForceField
      use ForceFieldPara_Vashishta
      use AcceptRates
      implicit none
      integer :: AllocateStatus
      integer :: i, j 
      ! print *, "nAtomTypes", nAtomTypes
      ALLOCATE (r_min(1:nAtomTypes), STAT=AllocateStatus)
      ALLOCATE (r_min_sq(1:nAtomTypes), STAT=AllocateStatus)
      ALLOCATE (r_min_tab(1:nAtomTypes, 1:nAtomTypes), STAT=AllocateStatus)

      ! 2-body parameters 
      ALLOCATE (H_tab(1:nAtomTypes, 1:nAtomTypes), STAT=AllocateStatus)
      ALLOCATE (eta_tab(1:nAtomTypes, 1:nAtomTypes), STAT=AllocateStatus)
      ALLOCATE (r4s_tab(1:nAtomTypes, 1:nAtomTypes), STAT=AllocateStatus)
      ALLOCATE (D_tab(1:nAtomTypes, 1:nAtomTypes), STAT=AllocateStatus)
      ALLOCATE (W_tab(1:nAtomTypes, 1:nAtomTypes), STAT=AllocateStatus)
      ! ALLOCATE (r1s_tab(1:nAtomTypes, 1:nAtomTypes), STAT=AllocateStatus)
      ! ALLOCATE (rc_tab(1:nAtomTypes, 1:nAtomTypes), STAT=AllocateStatus)

      ALLOCATE (ep_tab(1:nAtomTypes, 1:nAtomTypes), STAT=AllocateStatus)
      ALLOCATE (sig_tab(1:nAtomTypes, 1:nAtomTypes), STAT=AllocateStatus)
      ALLOCATE (q_tab(1:nAtomTypes, 1:nAtomTypes), STAT=AllocateStatus)

      ep_tab = 0d0
      sig_tab = 0d0
      q_tab = 0d0
      r_min = 0d0
      r_min_sq = 0d0
      r_min_tab = 0d0

      ! check content of the tab 
      ! do i = 1, nAtomTypes
      !    do j = i, nAtomTypes 
      !       print *, "-->", i, j, H_tab(i, j)
      !    end do 
      ! end do 

      ALLOCATE (totalMass(1:nMolTypes), STAT=AllocateStatus)

      ALLOCATE (nAtoms(1:nMolTypes), STAT=AllocateStatus)
      ALLOCATE (nIntraNonBond(1:nMolTypes), STAT=AllocateStatus)
      ALLOCATE (nBonds(1:nMolTypes), STAT=AllocateStatus)
      ALLOCATE (nAngles(1:nMolTypes), STAT=AllocateStatus)
      ALLOCATE (nTorsional(1:nMolTypes), STAT=AllocateStatus)
      ALLOCATE (nImproper(1:nMolTypes), STAT=AllocateStatus)

      nAtoms = 0
      nIntraNonBond = 0
      nBonds = 0
      nAngles = 0
      nTorsional = 0
      nImproper = 0

      ALLOCATE (acptTrans(1:nMolTypes), STAT=AllocateStatus)
      ALLOCATE (atmpTrans(1:nMolTypes), STAT=AllocateStatus)
      ALLOCATE (acptRot(1:nMolTypes), STAT=AllocateStatus)
      ALLOCATE (atmpRot(1:nMolTypes), STAT=AllocateStatus)
      ALLOCATE (acptSwapIn(1:nMolTypes), STAT=AllocateStatus)
      ALLOCATE (acptSwapIn(1:nMolTypes), STAT=AllocateStatus)
      ALLOCATE (atmpSwapIn(1:nMolTypes), STAT=AllocateStatus)
      ALLOCATE (acptSwapOut(1:nMolTypes), STAT=AllocateStatus)
      ALLOCATE (atmpSwapOut(1:nMolTypes), STAT=AllocateStatus)

      ALLOCATE (max_dist(1:nMolTypes), STAT=AllocateStatus)
      ALLOCATE (max_dist_single(1:nMolTypes), STAT=AllocateStatus)
      ALLOCATE (max_rot(1:nMolTypes), STAT=AllocateStatus)

      IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"

   end subroutine

!================================================================
   subroutine Allocate_LJ_Q
      use SimParameters
      use ForceField
      use ForceFieldPara_LJ_Q
      use AcceptRates
      implicit none
      integer :: AllocateStatus

      ALLOCATE (r_min(1:nAtomTypes), STAT=AllocateStatus)
      ALLOCATE (r_min_sq(1:nAtomTypes), STAT=AllocateStatus)
      ALLOCATE (r_min_tab(1:nAtomTypes, 1:nAtomTypes), STAT=AllocateStatus)

      ALLOCATE (ep_tab(1:nAtomTypes, 1:nAtomTypes), STAT=AllocateStatus)
      ALLOCATE (sig_tab(1:nAtomTypes, 1:nAtomTypes), STAT=AllocateStatus)
      ALLOCATE (q_tab(1:nAtomTypes, 1:nAtomTypes), STAT=AllocateStatus)

      ep_tab = 0d0
      sig_tab = 0d0
      q_tab = 0d0
      r_min = 0d0
      r_min_sq = 0d0
      r_min_tab = 0d0

      ALLOCATE (totalMass(1:nMolTypes), STAT=AllocateStatus)

      ALLOCATE (nAtoms(1:nMolTypes), STAT=AllocateStatus)
      ALLOCATE (nIntraNonBond(1:nMolTypes), STAT=AllocateStatus)
      ALLOCATE (nBonds(1:nMolTypes), STAT=AllocateStatus)
      ALLOCATE (nAngles(1:nMolTypes), STAT=AllocateStatus)
      ALLOCATE (nTorsional(1:nMolTypes), STAT=AllocateStatus)
      ALLOCATE (nImproper(1:nMolTypes), STAT=AllocateStatus)

      nAtoms = 0
      nIntraNonBond = 0
      nBonds = 0
      nAngles = 0
      nTorsional = 0
      nImproper = 0

      ALLOCATE (acptTrans(1:nMolTypes), STAT=AllocateStatus)
      ALLOCATE (atmpTrans(1:nMolTypes), STAT=AllocateStatus)
      ALLOCATE (acptRot(1:nMolTypes), STAT=AllocateStatus)
      ALLOCATE (atmpRot(1:nMolTypes), STAT=AllocateStatus)
      ALLOCATE (acptSwapIn(1:nMolTypes), STAT=AllocateStatus)
      ALLOCATE (acptSwapIn(1:nMolTypes), STAT=AllocateStatus)
      ALLOCATE (atmpSwapIn(1:nMolTypes), STAT=AllocateStatus)
      ALLOCATE (acptSwapOut(1:nMolTypes), STAT=AllocateStatus)
      ALLOCATE (atmpSwapOut(1:nMolTypes), STAT=AllocateStatus)

      ALLOCATE (max_dist(1:nMolTypes), STAT=AllocateStatus)
      ALLOCATE (max_dist_single(1:nMolTypes), STAT=AllocateStatus)
      ALLOCATE (max_rot(1:nMolTypes), STAT=AllocateStatus)

      IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"

   end subroutine

!================================================================
   subroutine Allocate_Pedone
      use SimParameters
      use ForceField
      use ForceFieldPara_Pedone
      use AcceptRates
      implicit none
      integer :: AllocateStatus

      ALLOCATE (pedoneData(1:nAtomTypes), STAT=AllocateStatus)
      ALLOCATE (atomData(1:nAtomTypes), STAT=AllocateStatus)
      ALLOCATE (nAtoms(1:nMolTypes), STAT=AllocateStatus)

      nAtoms = 1

      ALLOCATE (r_min(1:nAtomTypes), STAT=AllocateStatus)
      ALLOCATE (r_min_sq(1:nAtomTypes), STAT=AllocateStatus)
      ALLOCATE (r_min_tab(1:nAtomTypes, 1:nAtomTypes), STAT=AllocateStatus)
      ALLOCATE (bornRad(1:nAtomTypes), STAT=AllocateStatus)

      ALLOCATE (alpha_Tab(1:nAtomTypes, 1:nAtomTypes), STAT=AllocateStatus)
      ALLOCATE (D_Tab(1:nAtomTypes, 1:nAtomTypes), STAT=AllocateStatus)
      ALLOCATE (repul_tab(1:nAtomTypes, 1:nAtomTypes), STAT=AllocateStatus)
      ALLOCATE (rEq_tab(1:nAtomTypes, 1:nAtomTypes), STAT=AllocateStatus)
      ALLOCATE (q_tab(1:nAtomTypes, 1:nAtomTypes), STAT=AllocateStatus)

      ALLOCATE (totalMass(1:nAtomTypes), STAT=AllocateStatus)

      repul_tab = 0d0
      D_Tab = 0d0
      alpha_Tab = 0d0
      rEq_tab = 0d0
      q_tab = 0d0
      r_min_tab = 0d0

      ALLOCATE (acptTrans(1:nMolTypes), STAT=AllocateStatus)
      ALLOCATE (atmpTrans(1:nMolTypes), STAT=AllocateStatus)
      acptTrans = 0E0_dp
      atmpTrans = 1E-40_dp

      ALLOCATE (acptRot(1:nMolTypes), STAT=AllocateStatus)
      ALLOCATE (atmpRot(1:nMolTypes), STAT=AllocateStatus)
      acptRot = 0E0_dp
      atmpRot = 1E-40_dp

      ALLOCATE (acptSwapIn(1:nMolTypes), STAT=AllocateStatus)
      ALLOCATE (atmpSwapIn(1:nMolTypes), STAT=AllocateStatus)
      acptSwapIn = 0E0_dp
      atmpSwapIn = 1E-40_dp

      ALLOCATE (acptSwapOut(1:nMolTypes), STAT=AllocateStatus)
      ALLOCATE (atmpSwapOut(1:nMolTypes), STAT=AllocateStatus)
      acptSwapOut = 0E0_dp
      atmpSwapOut = 1E-40_dp

      ALLOCATE (acptInSize(1:maxMol), STAT=AllocateStatus)
      ALLOCATE (atmpInSize(1:maxMol), STAT=AllocateStatus)

      acptInSize = 0E0_dp
      atmpInSize = 1E-100_dp

      ALLOCATE (max_dist(1:nMolTypes), STAT=AllocateStatus)
      ALLOCATE (max_dist_single(1:nMolTypes), STAT=AllocateStatus)
      ALLOCATE (max_rot(1:nMolTypes), STAT=AllocateStatus)

      IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"

   end subroutine

!================================================================
   subroutine Allocate_Tersoff
      use SimParameters
      use ForceField
      use ForceFieldPara_Tersoff
      use AcceptRates
      implicit none
      integer :: AllocateStatus

      ALLOCATE (tersoffData(1:nAtomTypes), STAT=AllocateStatus)
      ALLOCATE (atomData(1:nAtomTypes), STAT=AllocateStatus)
      ALLOCATE (nAtoms(1:nMolTypes), STAT=AllocateStatus)

      nAtoms = 1

      ALLOCATE (r_min(1:nAtomTypes), STAT=AllocateStatus)
      ALLOCATE (r_min_sq(1:nAtomTypes), STAT=AllocateStatus)
      ALLOCATE (r_min_tab(1:nAtomTypes, 1:nAtomTypes), STAT=AllocateStatus)
      ALLOCATE (totalMass(1:nAtomTypes), STAT=AllocateStatus)

      ALLOCATE (acptTrans(1:nMolTypes), STAT=AllocateStatus)
      ALLOCATE (atmpTrans(1:nMolTypes), STAT=AllocateStatus)
      acptTrans = 0E0_dp
      atmpTrans = 1E-40_dp

      ALLOCATE (acptRot(1:nMolTypes), STAT=AllocateStatus)
      ALLOCATE (atmpRot(1:nMolTypes), STAT=AllocateStatus)
      acptRot = 0E0_dp
      atmpRot = 1E-40_dp

      ALLOCATE (acptSwapIn(1:nMolTypes), STAT=AllocateStatus)
      ALLOCATE (atmpSwapIn(1:nMolTypes), STAT=AllocateStatus)
      acptSwapIn = 0E0_dp
      atmpSwapIn = 1E-40_dp

      ALLOCATE (acptSwapOut(1:nMolTypes), STAT=AllocateStatus)
      ALLOCATE (atmpSwapOut(1:nMolTypes), STAT=AllocateStatus)
      acptSwapOut = 0E0_dp
      atmpSwapOut = 1E-40_dp

      ALLOCATE (acptInSize(1:maxMol), STAT=AllocateStatus)
      ALLOCATE (atmpInSize(1:maxMol), STAT=AllocateStatus)

      acptInSize = 0E0_dp
      atmpInSize = 1E-100_dp

      ALLOCATE (max_dist(1:nMolTypes), STAT=AllocateStatus)
      ALLOCATE (max_dist_single(1:nMolTypes), STAT=AllocateStatus)
      ALLOCATE (max_rot(1:nMolTypes), STAT=AllocateStatus)

      IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"

   end subroutine
!===================================================================================
   subroutine Read_Vashishta3(lineStore)
      use SimParameters
      use ForceField
      use ForceFieldPara_Vashishta
      implicit none
      character(len=maxLineLen), intent(in) :: lineStore(:)
      integer :: i, j, iLine, nLines

      nLines = size(lineStore)

      i = 0
      do iLine = 2, nLines - 1
         i = i + 1
         read (lineStore(iLine), *) atomDataV(i)%atmName, atomDataV(i)%B, atomDataV(i)%xi, &
                                    atomDataV(i)%r0, atomDataV(i)%C, atomDataV(i)%theta
      end do
   end subroutine
!===================================================================================
   subroutine Read_Vashishta2(lineStore)
      use SimParameters
      use ForceField
      use ForceFieldPara_Vashishta
      implicit none
      character(len=maxLineLen), intent(in) :: lineStore(:)
      integer :: i, j, iLine, nLines
      integer :: indx1, indx2 
      real(dp) :: H, eta, r4s, D, W
      character(len=20) :: atmName

      nLines = size(lineStore)

      i = 0
      do iLine = 2, nLines - 1
         i = i + 1
         read (lineStore(iLine), *) atmName, indx1, indx2, H, eta, r4s, D, W
         H_tab(indx1, indx2) = H
         eta_tab(indx1, indx2) = eta 
         r4s_tab(indx1, indx2) = r4s 
         D_tab(indx1, indx2) = D 
         W_tab(indx1, indx2) = W
         
         ! OH == HO 
         if (indx1 .ne. indx2) then
            H_tab(indx2, indx1) = H
            eta_tab(indx2, indx1) = eta 
            r4s_tab(indx2, indx1) = r4s 
            D_tab(indx2, indx1) = D 
            W_tab(indx2, indx1) = W
         end if
         
      end do
   end subroutine
!===================================================================================
   subroutine Read_Vashishta(lineStore)
      use SimParameters
      use ForceField
      use ForceFieldPara_Vashishta
      implicit none
      character(len=maxLineLen), intent(in) :: lineStore(:)
      integer :: i, j, iLine, nLines
      real(dp), parameter :: qqr2e = 14.399645000024973 ! e^2/(4*pi*eps0), eV*A 

      nLines = size(lineStore)

      i = 0
      do iLine = 2, nLines - 1
         i = i + 1
         read (lineStore(iLine), *) atomData(i)%atmName, atomData(i)%Symb, atomData(i)%q, atomData(i)%mass
      end do
!      Generate the look up tables for the inter molecular interactions
      if (echoInput) then
         write (35, *) "---------------------------------------------"
         write (35, *) "Interaction Table"
         write (35, *) " i ", " j ", " q(i) ", " q(j) ", " qqr2e "
      end if
      do i = 1, nAtomTypes
         do j = i, nAtomTypes
            q_tab(i, j) = atomData(i)%q*atomData(j)%q * qqr2e
            q_tab(j, i) = q_tab(i, j)

            if (echoInput) then
               write (35, *) i, j, atomData(i)%q, atomData(j)%q, qqr2e
            end if
         end do
      end do
      if (echoInput) then
         write (35, *) "---------------------------------------------"
         flush (35)
      end if

   end subroutine
!===================================================================================
   subroutine Read_LJ_Q(lineStore)
      use SimParameters
      use ForceField
      use ForceFieldPara_LJ_Q
      implicit none
      character(len=maxLineLen), intent(in) :: lineStore(:)
      integer :: i, j, iLine, nLines

      nLines = size(lineStore)

      i = 0
      do iLine = 2, nLines - 1
         i = i + 1
         read (lineStore(iLine), *) atomData(i)%atmName, atomData(i)%Symb, atomData(i)%ep, &
            atomData(i)%sig, atomData(i)%q, atomData(i)%mass
      end do
      ep_func => GeoMean_MixingFunc
      sig_func => Mean_MixingFunc
!      Generate the look up tables for the inter molecular interactions
      if (echoInput) then
         write (35, *) "---------------------------------------------"
         write (35, *) "Interaction Table"
         write (35, *) " i ", " j ", " eps ", " sig ", " q "
      end if
      do i = 1, nAtomTypes
         do j = i, nAtomTypes
            ep_tab(i, j) = 4d0*ep_func(atomData(i)%ep, atomData(j)%ep)
            sig_tab(i, j) = sig_func(atomData(i)%sig, atomData(j)%sig)**2
            q_tab(i, j) = atomData(i)%q*atomData(j)%q*coulombConst

            ep_tab(j, i) = ep_tab(i, j)
            sig_tab(j, i) = sig_tab(i, j)
            q_tab(j, i) = q_tab(i, j)
            if (echoInput) then
               write (35, *) i, j, ep_tab(i, j)/4d0, sqrt(sig_tab(i, j)), q_tab(i, j)
            end if
         end do
      end do
      if (echoInput) then
         write (35, *) "---------------------------------------------"
         flush (35)
      end if

   end subroutine

!===================================================================================
   subroutine Read_Pedone(lineStore)
      use SimParameters
      use ForceField
      use ForceFieldPara_Pedone
      implicit none
      character(len=maxLineLen), intent(in) :: lineStore(:)
      integer :: i, j, iLine, nLines

      nLines = size(lineStore)

      bornRad = 0E0_dp

      i = 0
      do iLine = 2, nLines - 1
         i = i + 1
!        write(*,*)  lineStore(iLine)
         read (lineStore(iLine), *) pedoneData(i)%atmName, pedoneData(i)%Symb, pedoneData(i)%repul, pedoneData(i)%rEq, &
            pedoneData(i)%alpha, pedoneData(i)%delta, pedoneData(i)%q, pedoneData(i)%mass, bornRad(i)
      end do

      do i = 1, nAtomTypes
         if (echoInput) then
            write (35, *) pedoneData(i)%atmName, pedoneData(i)%Symb, pedoneData(i)%repul, pedoneData(i)%rEq, &
               pedoneData(i)%alpha, pedoneData(i)%delta, pedoneData(i)%q, pedoneData(i)%mass, bornRad(i)
         end if
         pedoneData(i)%repul = pedoneData(i)%repul*convEng
         pedoneData(i)%delta = pedoneData(i)%delta*convEng
         pedoneData(i)%rEq = pedoneData(i)%rEq*convDist
         bornRad(i) = bornRad(i)*convDist
         atomData(i)%atmName = pedoneData(i)%atmName
         atomData(i)%Symb = pedoneData(i)%Symb
      end do

      q_tab = 0d0
      do i = 1, nAtomTypes
         do j = i, nAtomTypes
            q_tab(i, j) = pedoneData(i)%q*pedoneData(j)%q*coulombConst
            q_tab(j, i) = q_tab(i, j)
            if (echoInput) then
               write (35, *) i, j, q_tab(i, j)
            end if
         end do
      end do

      repul_tab = 0d0
      alpha_Tab = 0d0
      rEq_tab = 0d0
      D_Tab = 0d0
      do i = 1, nAtomTypes
         repul_tab(i, 1) = pedoneData(i)%repul
         alpha_Tab(i, 1) = pedoneData(i)%alpha
         rEq_tab(i, 1) = pedoneData(i)%rEq
         D_Tab(i, 1) = pedoneData(i)%delta

         repul_tab(1, i) = repul_tab(i, 1)
         alpha_Tab(1, i) = alpha_Tab(i, 1)
         rEq_tab(1, i) = rEq_tab(i, 1)
         D_Tab(1, i) = D_Tab(i, 1)

      end do

      if (echoInput) then
         do i = 1, nAtomTypes
            do j = 1, nAtomTypes
               write (35, *) i, j, repul_tab(i, j), alpha_Tab(i, j), rEq_tab(i, j), D_Tab(i, j), q_tab(i, j)
            end do
         end do
      end if

      if (any(bornRad .ne. 0E0_dp)) then
         implcSolvent = .true.
      else
         implcSolvent = .false.
      end if

   end subroutine

!===================================================================================
   subroutine Read_Tersoff(lineStore)
      use SimParameters
      use ForceField
      use ForceFieldPara_Tersoff
      implicit none
      character(len=maxLineLen), intent(in) :: lineStore(:)
      integer :: i, j, iLine, nLines

      nLines = size(lineStore)

      i = 0
      do iLine = 2, nLines - 1
         i = i + 1
!        write(*,*)  lineStore(iLine)
         read (lineStore(iLine), *) tersoffData(i)%atmName, tersoffData(i)%Symb, tersoffData(i)%A, &
            tersoffData(i)%B, tersoffData(i)%c, tersoffData(i)%d, tersoffData(i)%n, &
            tersoffData(i)%lam1, tersoffData(i)%lam2, tersoffData(i)%h, tersoffData(i)%R, &
            tersoffData(i)%D2, tersoffData(i)%beta, tersoffData(i)%mass
      end do

      do i = 1, nAtomTypes
         if (echoInput) then
            write (35, *) tersoffData(i)%atmName, tersoffData(i)%Symb, tersoffData(i)%A, &
               tersoffData(i)%B, tersoffData(i)%c, tersoffData(i)%d, tersoffData(i)%n, &
               tersoffData(i)%lam1, tersoffData(i)%lam2, tersoffData(i)%h, tersoffData(i)%R, &
               tersoffData(i)%D2, tersoffData(i)%beta, tersoffData(i)%mass
         end if
         tersoffData(i)%A = tersoffData(i)%A*convEng
         tersoffData(i)%B = tersoffData(i)%B*convEng

         tersoffData(i)%lam1 = tersoffData(i)%lam1/convDist
         tersoffData(i)%lam2 = tersoffData(i)%lam2/convDist

         tersoffData(i)%R = tersoffData(i)%R*convDist
         tersoffData(i)%D2 = tersoffData(i)%D2*convDist

         atomData(i)%atmName = tersoffData(i)%atmName
         atomData(i)%Symb = tersoffData(i)%Symb
      end do

   end subroutine
!===================================================================================
   subroutine LJ_SetFlags
      use SimParameters
      use Coords
      use ForceField
      use ForceFieldPara_LJ_Q
      use PairStorage, only: SetStorageFlags, distStorage, rPair, useDistStore
      implicit none
      integer :: iType, iMol, iAtom
      integer :: jType, jMol, jAtom
      integer :: atmType1, atmType2, globIndx1, globIndx2
      real(dp) :: q, ep

      call IntegrateBendAngleProb

      if (useDistStore) then
         call SetStorageFlags(q_tab)
         do iType = 1, nMolTypes
            do iMol = 1, NMAX(iType)
               do iAtom = 1, nAtoms(iType)
                  atmType1 = atomArray(iType, iAtom)
                  globIndx1 = MolArray(iType)%mol(iMol)%globalIndx(iAtom)
                  do jType = 1, nMolTypes
                     do jMol = 1, NMAX(jType)
                        do jAtom = 1, nAtoms(jType)
                           atmType2 = atomArray(jType, jAtom)
                           globIndx2 = MolArray(jType)%mol(jMol)%globalIndx(jAtom)
                           ep = ep_tab(atmType1, atmType2)
                           q = q_tab(atmType1, atmType2)
                           rPair(globIndx1, globIndx2)%p%usePair = .true.
                           if (abs(q) .gt. 1E-16_dp) then
                              rPair(globIndx1, globIndx2)%p%storeRValue = .true.
                           end if
                           if (ep .eq. 0E0_dp) then
                              if (q .eq. 0E0_dp) then
                                 rPair(globIndx1, globIndx2)%p%usePair = .false.

                              end if
                           end if
                        end do
                     end do
                  end do
               end do
            end do
         end do
      end if

   end subroutine
!===================================================================================
   subroutine Vashishta_SetFlags
      use SimParameters
      use Coords
      use ForceField
      use ForceFieldPara_Vashishta
      use PairStorage, only: SetStorageFlags, distStorage, rPair, useDistStore
      implicit none
      integer :: iType, iMol, iAtom
      integer :: jType, jMol, jAtom
      integer :: atmType1, atmType2, globIndx1, globIndx2
      real(dp) :: q, ep

      ! call IntegrateBendAngleProb

      ! if (useDistStore) then
      !    call SetStorageFlags(q_tab)
      !    do iType = 1, nMolTypes
      !       do iMol = 1, NMAX(iType)
      !          do iAtom = 1, nAtoms(iType)
      !             atmType1 = atomArray(iType, iAtom)
      !             globIndx1 = MolArray(iType)%mol(iMol)%globalIndx(iAtom)
      !             do jType = 1, nMolTypes
      !                do jMol = 1, NMAX(jType)
      !                   do jAtom = 1, nAtoms(jType)
      !                      atmType2 = atomArray(jType, jAtom)
      !                      globIndx2 = MolArray(jType)%mol(jMol)%globalIndx(jAtom)
      !                      ep = ep_tab(atmType1, atmType2)
      !                      q = q_tab(atmType1, atmType2)
      !                      rPair(globIndx1, globIndx2)%p%usePair = .true.
      !                      if (abs(q) .gt. 1E-16_dp) then
      !                         rPair(globIndx1, globIndx2)%p%storeRValue = .true.
      !                      end if
      !                      if (ep .eq. 0E0_dp) then
      !                         if (q .eq. 0E0_dp) then
      !                            rPair(globIndx1, globIndx2)%p%usePair = .false.

      !                         end if
      !                      end if
      !                   end do
      !                end do
      !             end do
      !          end do
      !       end do
      !    end do
      ! end if

   end subroutine

!===================================================================================
   subroutine Pedone_SetFlags
      use SimParameters
      use ForceField
      use ForceFieldPara_Pedone
      use PairStorage, only: SetStorageFlags, useDistStore
      implicit none

      if (useDistStore) then
         call SetStorageFlags(q_tab)
      end if

   end subroutine

!===================================================================================
   subroutine Tersoff_SetFlags
      use SimParameters
      use ForceField
      use ForceFieldPara_Pedone
      use PairStorage, only: TurnOnAllStorageFlags, useDistStore
      implicit none

      if (useDistStore) then
         call TurnOnAllStorageFlags
      end if

   end subroutine
!===================================================================================
   subroutine FindMax(lineStore, targetCommand, commandMax)
      use SimParameters
      use ForceField
      use ForceFieldPara_LJ_Q
      implicit none
      character(len=maxLineLen), intent(in) :: lineStore(:)
      character(len=*), intent(in) :: targetCommand
      integer, intent(out) :: commandMax

      integer :: intValue
      character(len=25) :: curCommand, dummy

      integer :: iLine, nLines, lineStat

      nLines = size(lineStore)
      commandMax = 0
      do iLine = 1, nLines
         call GetCommand(lineStore(iLine), curCommand, lineStat)
         if (trim(adjustl(curCommand)) .eq. trim(adjustl(targetCommand))) then
            read (lineStore(iLine), *) dummy, intValue
            commandMax = max(commandMax, intValue)
         end if
      end do

   end subroutine

!================================================================================
end module
