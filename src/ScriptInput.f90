!========================================================
module ScriptInput
   integer, parameter :: maxLineLen = 500
contains
!========================================================
   subroutine Script_ReadParameters(seed, screenEcho)
      use AnalysisMain, only: ScriptAnalysisInput
      use CBMC_Variables
      use Constants
      use Coords
      use CoordinateFunctions
      use EnergyTables
      use ForceField
      use ForceFieldInput, only: SetForcefieldType, ScriptForcefield, fieldTypeSet
      use MoveTypeModule, only: ScriptInput_MCMove
      use ParallelVar
      use SimParameters
      use UmbrellaSamplingNew, only: useUmbrella, ScriptInput_Umbrella
      use Units
      use WHAM_Functions
      implicit none
      logical, intent(OUT)  :: screenEcho
      integer, intent(OUT) :: seed
!      integer(kind=8), intent(OUT) :: ncycle,nmoves
      integer :: i, ii, j, nArgs
      integer :: iLine, lineStat, AllocateStat
      integer :: nLines, nForceLines, lineBuffer
      integer, allocatable :: lineNumber(:), ffLineNumber(:)
      character(len=maxLineLen), allocatable :: lineStore(:)
      character(len=maxLineLen), allocatable :: forcefieldStore(:)
      character(len=25) :: command, command2, dummy
      character(len=maxLineLen) :: command_ff, command_cf
      ! character(len=25) :: command,dummy
      ! character(len=100) :: command2 
      character(len=maxLineLen) :: fileName
      character(len=maxLineLen) :: forcefieldFile

!      Get the filename from the command line.
      nArgs = command_argument_count()
      if (nArgs > 2) then
         stop "This program only takes one argument"
      elseif (nArgs == 1) then
         call get_command_argument(1, fileName)
         call LoadFile(lineStore, nLines, lineNumber, fileName)
      elseif (nArgs == 2) then
         call get_command_argument(2, fileName)
         call LoadFile(lineStore, nLines, lineNumber, fileName)
      elseif (nArgs == 0) then
         write (*, *) "ERROR! No Input File has been given!"
         stop
      end if
!      Read in the script file

!      fileName = "ScriptTest.dat"

!      This block counts the number of lines in the input file to determine how large the lineStorage array needs to be.

      call setDefaults(seed, screenEcho)

      lineBuffer = 0
      do iLine = 1, nLines
         if (lineBuffer .gt. 0) then
            lineBuffer = lineBuffer - 1
            cycle
         end if
         lineStat = 0
!        write(*,*) trim(adjustl(lineStore(iLine)))
         call getCommand(lineStore(iLine), command, lineStat)
         call LowerCaseLine(command)
!         If line is empty or commented, move to the next line.
         if (lineStat .eq. 1) then
            cycle
         end if

         select case (trim(adjustl(command)))
         case ("analysis")
            call FindCommandBlock(iLine, lineStore, "end_analysis", lineBuffer)
            call ScriptAnalysisInput(lineStore(iLine:iLine + lineBuffer))
         case ("set")
            call setVariable(lineStore(iLine), seed, screenEcho, lineStat)
            if (lineStat .eq. -1) then
               write (*, "(A,2x,I10)") "ERROR! Unknown Variable Name on Line", lineNumber(iLine)
               write (*, *) lineStore(iLine)
               stop
            end if
         case ("movetypes")
            call FindCommandBlock(iLine, lineStore, "end_movetypes", lineBuffer)
            call ScriptInput_MCMove(lineStore(iLine:iLine + lineBuffer))
         case ("umbrella")
            call FindCommandBlock(iLine, lineStore, "end_umbrella", lineBuffer)
            call ScriptInput_Umbrella(lineStore(iLine:iLine + lineBuffer))
         case ("forcefieldfile")
            read (lineStore(iLine), *) dummy, command_ff
            forcefieldFile = trim(adjustl(command_ff))
            call LoadFile(forcefieldStore, nForceLines, ffLineNumber, forcefieldFile)
            call ScriptForcefield(forcefieldStore)
            if (allocated(forcefieldStore)) then
               deallocate (forcefieldStore)
            end if
         case ("forcefieldtype")
            read (lineStore(iLine), *) dummy, command2
            call LowerCaseLine(command2)
            call SetForcefieldType(command2)
            fieldTypeSet = .true.
         case ("boundary")
            if (nMolTypes .eq. 0) then
               write (*, *) "INPUT ERROR! Boundary is called before the number of molecule types has been assigned"
               stop
            end if
            call FindCommandBlock(iLine, lineStore, "end_boundary", lineBuffer)
            allocate (Eng_Critr(1:nMolTypes, 1:nMolTypes), STAT=AllocateStat)
            ii = 0
            do i = iLine + 1, iLine + lineBuffer - 1
               ii = ii + 1
               read (lineStore(i), *) (Eng_Critr(ii, j), j=1, nMolTypes)
            end do
         case ("biasalpha")
            if (nMolTypes .eq. 0) then
               write (*, *) "INPUT ERROR! BiasAlpha is called before the number of molecule types has been assigned"
               stop
            end if
            call FindCommandBlock(iLine, lineStore, "end_biasalpha", lineBuffer)
            allocate (biasAlpha(1:nMolTypes, 1:nMolTypes), stat=AllocateStat)
            if (lineBuffer - 1 .gt. nMolTypes) then
               write (*, *) "ERROR! User has given the biasalpha command far too many lines.  BiasAlpha takes a number of lines"
               write (*, *) "equal to the number of molecule types."
               stop
            end if
            ii = 0
            do i = iLine + 1, iLine + lineBuffer - 1
               ii = ii + 1
               read (lineStore(i), *) (biasAlpha(ii, j), j=1, nMolTypes)
               do j = 1, nMolTypes
                  biasAlpha(ii, j) = biasAlpha(ii, j)/temperature
               end do
            end do
         case ("configfile")
            read (lineStore(iLine), *) dummy, command_cf
            call ReadInitialConfiguration(command_cf)
            write(*,*) "Finished Reading Initial Configuration."
            call RecenterCoordinates
         case ("gasconfigfile")
            read (lineStore(iLine), *) dummy, command_cf
            call ReadInitialGasPhase(command_cf)
            write(*,*) "Finished Setting up Gas Phase Configuration."
         case default
            write (*, "(A,2x,I10)") "ERROR! Unknown Command on Line", lineNumber(iLine)
            write (*, *) lineStore(iLine)
            stop
         end select
      end do
!   write(*,*) "Finished Reading Input Script."
!   deallocate(lineStore)
!   call ReadInitialConfiguration
!   write(*,*) "Finished Reading Initial Configuration."
!   call RecenterCoordinates
!   call ReadInitialGasPhase
!   write(*,*) "Finished Setting up Gas Phase Configuration."
   if (useWHAM) then
      if (.not. useUmbrella) then
         write (*, *) "ERROR! The WHAM method can not be used if no Umbrella sampling variables are given!"
         stop
      end if
      nWhamItter = ceiling(dble(ncycle)/dble(intervalWHAM))
      call WHAM_Initialize
   end if

   end subroutine
!========================================================
   subroutine setVariable(line, seed, screenEcho, lineStat)
      use VarPrecision
      use SimParameters
      use CBMC_Variables
      use Coords
      use EnergyTables
      use WHAM_Module
      use AcceptRates
      use Units
      use PairStorage, only: useDistStore
      implicit none
      character(len=maxLineLen), intent(in) :: line
      logical, intent(out) :: screenEcho
      integer, intent(out) :: seed, lineStat

      character(len=30) :: dummy, command!, stringValue
!      character(len=15) :: fileName
      logical :: logicValue
      integer :: j
      integer :: intValue, AllocateStat
      real(dp) :: realValue

      lineStat = 0

      read (line, *) dummy, command
      call LowerCaseLine(command)
      select case (trim(adjustl(command)))
      case ("avbmc_distance")
         read (line, *) dummy, command, realValue
         Dist_Critr = realValue
         Dist_Critr_sq = realValue**2
      case ("calcpressure")
         read (line, *) dummy, command, logicValue
         calcPressure = logicValue
      case ("cycles")
         read (line, *) dummy, command, realValue
         nCycle = nint(realValue)
      case ("distancecriteria")
         read (line, *) dummy, command, logicValue
         distCriteria = logicValue
      case ("distancestorage")
         read (line, *) dummy, command, logicValue
         useDistStore = logicValue
      case ("r1s")
         read (line, *) dummy, command, realValue
         r1s = realValue
      case ("rc")
         read (line, *) dummy, command, realValue
         rc = realValue
      case ("gasdensity")
         if (.not. allocated(gas_dens)) then
            write (*, *) "INPUT ERROR! GasDensity is called before the number of molecule types has been assigned"
            stop
         end if
         read (line, *) dummy, command, (gas_dens(j), j=1, nMolTypes)
      case ("moves")
         read (line, *) dummy, command, realValue
         ncycle2 = nint(realValue)
      case ("moleculetypes")
         read (line, *) dummy, command, realValue
         nMolTypes = nint(realValue)
         allocate (NPART(1:nMolTypes), STAT=AllocateStat)
         allocate (NPART_new(1:nMolTypes), STAT=AllocateStat)
         allocate (NMIN(1:nMolTypes), STAT=AllocateStat)
         allocate (NMAX(1:nMolTypes), STAT=AllocateStat)
         allocate (gas_dens(1:nMolTypes), STAT=AllocateStat)
         allocate (nRosenTrials(1:nMolTypes), STAT=AllocateStat)
         NMIN = 0
         NMAX = 0
         gas_dens = 0
         nRosenTrials = 1
         allocate (biasAlpha(1:nMolTypes, 1:nMolTypes), STAT=AllocateStat)
         biasAlpha = 0E0_dp
      case ("molmin")
         if (.not. allocated(NMIN)) then
            write (*, *) "INPUT ERROR! molmin is called before the number of molecule types has been assigned"
            stop
         end if
         read (line, *) dummy, command, (NMIN(j), j=1, nMolTypes)
      case ("molmax")
         if (.not. allocated(NMAX)) then
            write (*, *) "INPUT ERROR! molmax is called before the number of molecule types has been assigned"
            stop
         end if
         read (line, *) dummy, command, (NMAX(j), j=1, nMolTypes)
         maxMol = sum(NMAX)
!          ALLOCATE (acptInSize(1:maxMol), STAT = AllocateStat)
!          ALLOCATE (atmpInSize(1:maxMol), STAT = AllocateStat)
      case ("rosentrials")
         if (.not. allocated(nRosenTrials)) then
            write (*, *) "INPUT ERROR! molmax is called before the number of molecule types has been assigned"
            stop
         end if
         read (line, *) dummy, command, (nRosenTrials(j), j=1, nMolTypes)
         maxRosenTrial = maxval(nRosenTrials)
         allocate (rosenTrial(1:maxRosenTrial))
      case ("temperature")
         read (line, *) dummy, command, realValue
         temperature = realValue
         beta = 1d0/temperature
      case ("screenecho")
         read (line, *) dummy, command, logicValue
         screenEcho = logicValue
      case ("rng_seed")
         read (line, *) dummy, command, intValue
         seed = intValue
         if (seed .lt. 0) then
            call system_clock(intValue)
            seed = mod(intValue, 10000)
         end if
      case ("softcutoff")
         read (line, *) dummy, command, realValue
         softCutoff = realValue
      case ("screen_outfreq")
         read (line, *) dummy, command, intValue
         outFreq_Screen = intValue
      case ("trajectory_outfreq")
         read (line, *) dummy, command, intValue
         outFreq_Traj = intValue
      case ("multipleinputconfig")
         read (line, *) dummy, command, logicValue
         multipleInput = logicValue
      case ("out_energyunits")
         read (line, *) dummy, command, outputEngUnits
         outputEConv = FindEngUnit(outputEngUnits)
      case ("out_distunits")
         read (line, *) dummy, command, outputLenUnits
         outputLenConv = FindLengthUnit(outputLenUnits)
      case ("usewham")
         read (line, *) dummy, command, logicValue
         useWHAM = logicValue
      case ("whamseglength")
         read (line, *) dummy, command, intValue
         intervalWHAM = intValue
      case ("whammaxiteration")
         read (line, *) dummy, command, intValue
         maxSelfConsist = intValue
      case ("whamdgreplace")
         read (line, *) dummy, command, intValue
         whamEstInterval = intValue
      case ("whamequilcycle")
         read (line, *) dummy, command, intValue
         equilInterval = intValue
      case ("whamtol")
         read (line, *) dummy, command, realValue
         tolLimit = realValue
      case default
         lineStat = -1
      end select

   end subroutine
!========================================================
   subroutine setDefaults(seed, screenEcho)
      use VarPrecision
      use SimParameters
      use CBMC_Variables
      use Coords
      use EnergyTables
      use WHAM_Module
      use AcceptRates
      use Units
      implicit none
      logical, intent(out) :: screenEcho
      integer, intent(out) :: seed

      seed = -1

      Dist_Critr = 0E0_dp
      Dist_Critr_sq = 0E0_dp
      nCycle = 0
      ncycle2 = 0
      maxMol = 0

      maxRosenTrial = 1

      temperature = 0E0_dp
      beta = 0E0_dp

      screenEcho = .true.
      softCutoff = 100E0_dp

      outFreq_Screen = 1000
      outFreq_Traj = 1000
      multipleInput = .false.

      outputEConv = FindEngUnit("kb")
      outputLenConv = FindLengthUnit("ang")
      useWHAM = .false.
      intervalWHAM = 10000
      maxSelfConsist = 100
      whamEstInterval = -1
      equilInterval = 10
      tolLimit = 1E-5_dp

   end subroutine
!========================================================
   subroutine variableSafetyCheck
      implicit none

   end subroutine

!========================================================
   subroutine LoadFile(lineArray, nLines, lineNumber, fileName)
      use ParallelVar, only: myid
      use SimParameters, only: echoInput
      implicit none
      character(len=maxLineLen), allocatable, intent(inout) :: lineArray(:)
      character(len=maxLineLen), intent(in) :: fileName
      integer, allocatable, intent(inout) :: lineNumber(:)

      character(len=maxLineLen), allocatable :: rawLines(:)
      character(len=25) :: command
      integer, intent(out) :: nLines
      integer :: i, iLine, nRawLines, lineStat, AllocateStat, InOutStat

      open (unit=54, file=trim(adjustl(fileName)), status='OLD', iostat=InOutStat)
      ! print *, "FORCEFIELD NAME:", trim(adjustl(fileName))
      if (InOutStat .gt. 0) then
         if (myid .eq. 0) then
            write (*, *) "The file specified in either the input script or on the command line could not be opened."
            write (*, *) "Please double check to ensure the file exists and the name in the input is accurate."
            write (*, *) "  Attempted to open file:", trim(adjustl(fileName))
         end if
         stop
      end if

!      This block counts the number of lines in the input file to determine how large the lineStorage array needs to be.
      nRawLines = 0
      do iLine = 1, nint(1d7)
         read (54, *, iostat=lineStat)
         if (lineStat .lt. 0) then
            exit
         end if
         nRawLines = nRawLines + 1
      end do
      rewind (54)

!      Read in the file line by line
      allocate (rawLines(1:nRawLines), stat=AllocateStat)

      do iLine = 1, nRawLines
         read (54, "(A)") rawLines(iLine)
         if (echoInput) then
            write (35, *) rawLines(iLine)
         end if
!        call LowerCaseLine(rawLines(iLine))
!        write(*,*) rawLines(iLine)
      end do
      close (54)

      nLines = 0
      do iLine = 1, nRawLines
         lineStat = 0
         call getCommand(rawLines(iLine), command, lineStat)
!         If line is empty or commented, move to the next line.
         if (lineStat .eq. 0) then
            nLines = nLines + 1
         end if
      end do

      allocate (lineArray(1:nLines))
      allocate (lineNumber(1:nLines))
      i = 0
      do iLine = 1, nRawLines
         lineStat = 0
         call getCommand(rawLines(iLine), command, lineStat)
         if (lineStat .eq. 0) then
            i = i + 1
            lineArray(i) = rawLines(iLine)
            lineNumber(i) = iLine
!          write(*,*) lineNumber(i), lineArray(i)
         end if
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
!========================================================
   subroutine FindCommandBlock(iLine, lineStore, endCommand, lineBuffer)
      implicit none
      integer, intent(in) :: iLine
      character(len=maxLineLen), intent(in) :: lineStore(:)
      character(len=*), intent(in) :: endCommand
      integer, intent(out) :: lineBuffer
      logical :: found
      integer :: i, lineStat, nLines
      character(len=35) :: command

      command = " "
      nLines = size(lineStore)
      found = .false.
      do i = iLine + 1, nLines
         call GetCommand(lineStore(i), command, lineStat)
         call LowerCaseLine(command)
!        write(*,*)  dummy
         if (trim(adjustl(command)) .eq. trim(adjustl(endCommand))) then
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
end module
