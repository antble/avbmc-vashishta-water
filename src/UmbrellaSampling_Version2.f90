!==========================================================================
module UmbrellaSamplingNew
   use VarPrecision
   use SimpleDistPair
   use UmbrellaTypes
   implicit none
   private

   integer, parameter :: maxLineLen = 500

   logical :: useUmbrella = .false.
   logical :: UScreenOut
   logical :: energyAnalytics = .true.
   integer :: nBiasVariables = 0
   integer :: curUIndx, umbrellaLimit
   integer, allocatable :: curVarIndx
   integer, allocatable :: binIndx(:)
   integer, allocatable :: varMax(:), varMin(:)
   integer, allocatable :: binMax(:), binMin(:)
   integer, allocatable :: indexCoeff(:)
   integer, allocatable :: UArray(:)
   real(dp), allocatable :: UBias(:)
   real(dp), allocatable :: UHist(:)
   real(dp), allocatable :: UHistTotal(:)
   real(dp), allocatable :: UBinSize(:)
   real(dp), allocatable :: varValues(:)
   character(len=50):: inputFile
   character(len=10), allocatable :: outputFormat(:)
   character(len=100) :: screenFormat
   type(BiasVariablePointer), allocatable :: biasvar(:)
   type(BiasVariablePointer), allocatable :: biasvarnew(:)

   integer :: nDispFunc, nSwapInFunc, nSwapOutFunc
   type(DispUmbrellaArray), allocatable :: DispUmbrella(:)
   type(SwapInUmbrellaArray), allocatable :: SwapInUmbrella(:)
   type(SwapOutUmbrellaArray), allocatable :: SwapOutUmbrella(:)

   integer, parameter :: E_Bins = 1000
   real(dp), parameter :: dE = -100000/E_Bins
   real(dp), allocatable :: U_EAvg(:)
   real(dp), allocatable :: U_EHist(:, :)

!    public :: ReadInput_Umbrella
   public :: AllocateUmbrellaVariables, AllocateUmbrellaArray, UmbrellaHistAdd
   public :: useUmbrella, OutputUmbrellaHist, GetUmbrellaBias_Disp, findVarValues, getBiasIndex
   public :: nBiasVariables, umbrellaLimit, UBias, UHist, UBinSize, outputFormat, curUIndx
   public :: DispUmbrella, SwapInUmbrella, SwapOutUmbrella, biasVar, biasVarnew
   public :: GetUmbrellaBias_SwapIn, GetUmbrellaBias_SwapOut, ScreenOutputUmbrella, screenFormat
   public :: CheckInitialValues, energyAnalytics, OutputUmbrellaAnalytics
   public :: ScriptInput_Umbrella, GetUmbrellaBias_Temperature
   public :: inputFile
!==========================================================================================
contains
!==========================================================================================
   subroutine AllocateUmbrellaVariables
      implicit none
      integer :: AllocateStatus
      integer :: i

      allocate (binMin(1:nBiasVariables), STAT=AllocateStatus)
      allocate (binMax(1:nBiasVariables), STAT=AllocateStatus)
      allocate (UBinSize(1:nBiasVariables), STAT=AllocateStatus)

      allocate (biasvar(1:nBiasVariables), STAT=AllocateStatus)
      do i = 1, nBiasVariables
         biasvar(i)%varType = 0
         biasvar(i)%intVar => null()
         biasvar(i)%realVar => null()
      end do

      allocate (biasvarnew(1:nBiasVariables), STAT=AllocateStatus)
      do i = 1, nBiasVariables
         biasvarnew(i)%varType = 0
         biasvarnew(i)%intVar => null()
         biasvarnew(i)%realVar => null()
      end do

      allocate (binIndx(1:nBiasVariables), STAT=AllocateStatus)
      allocate (indexCoeff(1:nBiasVariables), STAT=AllocateStatus)
      allocate (outputFormat(1:nBiasVariables), STAT=AllocateStatus)
      do i = 1, nBiasVariables
         outputFormat(i) = " "
      end do

      allocate (DispUmbrella(1:nDispFunc), STAT=AllocateStatus)
      do i = 1, nDispFunc
         DispUmbrella(i)%func => null()
      end do

      allocate (SwapInUmbrella(1:nSwapInFunc), STAT=AllocateStatus)
      do i = 1, nSwapInFunc
         SwapInUmbrella(i)%func => null()
      end do

      allocate (SwapOutUmbrella(1:nSwapOutFunc), STAT=AllocateStatus)
      do i = 1, nSwapOutFunc
         SwapOutUmbrella(i)%func => null()
      end do

      IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
   end subroutine
!==========================================================================================
   subroutine ScriptInput_Umbrella(inputLines)
      use AnalysisMain, only: internalIndx, loadUmbArray, nAnalysisVar
      use MiscelaniousVars
      use SimpleDistPair, only: nDistPair, pairArrayIndx, CalcDistPairs_New, UmbrellaVar_DistPair
      use SimParameters, only: NMAX, NMIN, NPART, NPart_New, nMolTypes, maxMol, echoInput, NTotal, NTotal_New, temperature, TempNew
      use Q6Functions, only: q6ArrayIndx, CalcQ6_Disp, CalcQ6_SwapIn, CalcQ6_SwapOut
      use ParallelVar, only: nout
      use WHAM_Module, only: refBin, refSizeNumbers
      use Q6Functions, only: UmbrellaVar_Q6
      implicit none
      character(len=maxLineLen), intent(in) :: inputLines(:)
      integer :: nLines
      integer :: iUmbrella, AllocateStatus
      integer :: indxVar, stat
      integer :: iDisp, iSwapIn, iSwapOut
      real(dp) :: binSize, valMax, valMin
      real(dp), allocatable :: refVals(:)

      character(len=30) :: labelField
      character(len=30) :: umbrellaName

!      read(fileUnit, *) labelField, nBiasVariables
      nLines = size(inputLines)
      nBiasVariables = nLines - 3
      if (nBiasVariables .lt. 0) then
         write (*, *) "ERROR! The user has specified an invalid number of Umbrella Sampling Variables"
         write (*, *) labelField, nBiasVariables
         stop
      end if
      if (nBiasVariables .eq. 0) then
         useUmbrella = .false.
         write (nout, *) "Umbrella Sampling Used? = ", useUmbrella
         return
      else
         useUmbrella = .true.
         write (nout, *) "Umbrella Sampling Used? = ", useUmbrella
         write (nout, *) "Number of Umbrella Variables:", nBiasVariables
      end if

      allocate (refVals(1:nBiasVariables))

!     Read the file name for the initial Umbrealla Sampling Bias
      read (inputLines(1), *) labelField, inputFile

!     This first block performs an initial run through the input script to determine
!     how large the various pointer arrays must be to satisfy the various needs
!     of each

      nDispFunc = 0
      nSwapInFunc = 0
      nSwapOutFunc = 0
      do iUmbrella = 1, nBiasVariables
         read (inputLines(iUmbrella + 2), *) umbrellaName
         select case (trim(adjustl(umbrellaName)))
         case ("clustersize")
            continue
         case ("totalclustersize")
            continue
         case ("temperature")
            continue
         case ("pairdist")
            nDispFunc = nDispFunc + 1
            nSwapInFunc = nSwapInFunc + 1
            nSwapOutFunc = nSwapOutFunc + 1
         case ("q6")
            nDispFunc = nDispFunc + 1
            nSwapInFunc = nSwapInFunc + 1
            nSwapOutFunc = nSwapOutFunc + 1
         case ("analysisvar")
            nDispFunc = nDispFunc + 1
            nSwapInFunc = nSwapInFunc + 1
            nSwapOutFunc = nSwapOutFunc + 1
         case default
            write (*, *) "ERROR! Invalid variable type specified in input file"
            write (*, *) umbrellaName
            stop
         end select
      end do

      call AllocateUmbrellaVariables

!      In the next block all pointers used by the simulation are assigned to the user defined variables. And
!      parameters such as the largest and smallest allowed values for each variable are assigned.
!

      iDisp = 0
      iSwapIn = 0
      iSwapOut = 0
      do iUmbrella = 1, nBiasVariables
         read (inputLines(iUmbrella + 2), *) umbrellaName
!        write(*,*) inputLines(iUmbrella+2)
         select case (trim(adjustl(umbrellaName)))
         case ("clustersize")
            read (inputLines(iUmbrella + 2), *) labelField, indxVar
            if (indxVar .le. 0) then
               write (*, *) "Error! An invalid molecule type has been chosen!"
               write (*, *) "Defined Mol Types:", nMolTypes
               write (*, *) "Chosen Mol Type:", indxVar
               stop
            end if
            if (indxVar .gt. nMolTypes) then
               write (*, *) "Error! An invalid molecule type has been chosen!"
               write (*, *) "Defined Mol Types:", nMolTypes
               write (*, *) "Chosen Mol Type:", indxVar
               stop
            end if
            biasvar(iUmbrella)%varType = 1
            biasvar(iUmbrella)%intVar => NPart(indxVar)
            biasvarnew(iUmbrella)%varType = 1
            biasvarnew(iUmbrella)%intVar => NPart_New(indxVar)
            binMax(iUmbrella) = NMAX(indxVar)
            binMin(iUmbrella) = NMin(indxVar)
            UBinSize(iUmbrella) = 1E0_dp
            outputFormat(iUmbrella) = "2x,F5.1,"

!          iSwapIn = iSwapIn + 1
!          iSwapOut = iSwapOut + 1
         case ("totalclustersize")
            biasvar(iUmbrella)%varType = 1
            biasvar(iUmbrella)%intVar => NTotal
            biasvarnew(iUmbrella)%varType = 1
            biasvarnew(iUmbrella)%intVar => NTotal_New
            binMax(iUmbrella) = maxMol
            binMin(iUmbrella) = sum(NMin)
            UBinSize(iUmbrella) = 1E0_dp
            outputFormat(iUmbrella) = "2x,F5.1,"
         case ("temperature")
            read (inputLines(iUmbrella + 2), *) labelField, valMin, valMax, binSize
            biasvar(iUmbrella)%varType = 2
            biasvar(iUmbrella)%realVar => temperature
            biasvarnew(iUmbrella)%varType = 2
            biasvarnew(iUmbrella)%realVar => TempNew
            UBinSize(iUmbrella) = binSize
            binMin(iUmbrella) = nint(valMin/binSize)
            binMax(iUmbrella) = nint(valMax/binSize)
            outputFormat(iUmbrella) = "2x,F12.6,"
         case ("pairdist")
            indxVar = 0
            read (inputLines(iUmbrella + 2), *) labelField, indxVar, valMin, valMax, binSize
            if (nDistPair .le. 0) then
               write (*, *) "Error! An invalid distance variable has been chosen!"
               write (*, *) "Defined Distance Pairs:", nDistPair
               write (*, *) "Chosen Distance Pair:", indxVar
               stop
            end if
            if (indxVar .gt. nDistPair) then
               write (*, *) "Error! An invalid distance variable has been chosen!"
               write (*, *) "Defined Distance Pairs:", nDistPair
               write (*, *) "Chosen Distance Pair:", indxVar
               stop
            end if
            biasvar(iUmbrella)%varType = 2
            biasvar(iUmbrella)%realVar => miscCoord(pairArrayIndx(indxVar))
            biasvarnew(iUmbrella)%varType = 2
            biasvarnew(iUmbrella)%realVar => miscCoord_New(pairArrayIndx(indxVar))
            UBinSize(iUmbrella) = binSize
            binMin(iUmbrella) = nint(valMin/binSize)
            binMax(iUmbrella) = nint(valMax/binSize)
            outputFormat(iUmbrella) = "2x,F12.6,"

            iDisp = iDisp + 1
            DispUmbrella(iDisp)%func => CalcDistPairs_New
            iSwapIn = iSwapIn + 1
            SwapInUmbrella(iSwapIn)%func => CalcDistPairs_SwapIn
            iSwapOut = iSwapOut + 1
            SwapOutUmbrella(iSwapOut)%func => CalcDistPairs_SwapOut

         case ("q6")
            read (inputLines(iUmbrella + 2), *) labelField, valMin, valMax, binSize
            biasvar(iUmbrella)%varType = 2
            biasvar(iUmbrella)%realVar => miscCoord(q6ArrayIndx)
            biasvarnew(iUmbrella)%varType = 2
            biasvarnew(iUmbrella)%realVar => miscCoord_New(q6ArrayIndx)
            UBinSize(iUmbrella) = binSize
            binMin(iUmbrella) = nint(valMin/binSize)
            binMax(iUmbrella) = nint(valMax/binSize)
!          outputFormat(iUmbrella) = "2x,F12.6,"

            iDisp = iDisp + 1
            DispUmbrella(iDisp)%func => CalcQ6_Disp
            iSwapIn = iSwapIn + 1
            SwapInUmbrella(iSwapIn)%func => CalcQ6_SwapIn
            iSwapOut = iSwapOut + 1
            SwapOutUmbrella(iSwapOut)%func => CalcQ6_SwapOut
!          write(*,*) labelField, valMin, valMax, binSize
         case ("analysisvar")
            read (inputLines(iUmbrella + 2), *) labelField, indxVar, valMin, valMax, binSize
            if (indxVar .gt. nAnalysisVar) then
               write (*, *) "ERROR! The analysis index is above the maximum index of", nAnalysisVar
               write (*, *) inputLines(iUmbrella)
               stop
            end if
!          write(*,*) "BLAH!"
            call loadUmbArray(indxVar)%func(iUmbrella, indxVar, biasVar, biasVarNew, outputFormat, &
                                            iDisp, DispUmbrella, iSwapIn, SwapInUmbrella, iSwapOut, SwapOutUmbrella)
            UBinSize(iUmbrella) = binSize
            binMin(iUmbrella) = nint(valMin/binSize)
            binMax(iUmbrella) = nint(valMax/binSize)
            write (*, *) valMin, valMax, binSize
            write (*, *) binMin(iUmbrella), binMax(iUmbrella)
         case default
            write (*, *) "ERROR! Invalid variable type specified in input file"
            write (*, *) umbrellaName
            stop
         end select
      end do
      write (screenFormat, *) "(", (trim(adjustl(outputFormat(iUmbrella))), iUmbrella=1, nBiasVariables), "1x)"
      allocate (varValues(1:nBiasVariables))
      allocate (UArray(1:nBiasVariables))
      call AllocateUmbrellaArray
      call ReadInitialBias

!      Use the input to specify the reference bin
      allocate (refSizeNumbers(1:nBiasVariables), STAT=AllocateStatus)
!      write(*, *) trim(adjustl(inputLines(2)))
      read (inputLines(2), *) labelField, (refVals(iUmbrella), iUmbrella=1, nBiasVariables)
      call getUIndexArray(refVals, refBin, stat)
      if (stat .ne. 0) then
         if (stat .eq. 1) then
          write (*, *) "ERROR! Reference bin is above the largest allowed value for one or more of the umbrella sampling variables!"
            write (*, *) trim(adjustl(inputLines(2)))
            stop
         elseif (stat .eq. -1) then
         write (*, *) "ERROR! Reference bin is below the smallest allowed value for one or more of the umbrella sampling variables!"
            write (*, *) trim(adjustl(inputLines(2)))
            stop
         else
            write (*, *) "ERROR! An invalid value for the umbrella sampling reference bin was given!"
            write (*, *) trim(adjustl(inputLines(2)))
            stop
         end if
      end if
      do iUmbrella = 1, nBiasVariables
         refSizeNumbers(iUmbrella) = refVals(iUmbrella)
      end do
!      write(*,*) refVals

      deallocate (refVals)

   end subroutine
!==========================================================================================
   subroutine AllocateUmbrellaArray
      use ParallelVar
      implicit none
      integer :: i, j
      integer :: AllocateStatus

      indexCoeff(1) = 1
      do i = 2, nBiasVariables
         indexCoeff(i) = 1
         do j = 1, i - 1
            indexCoeff(i) = indexCoeff(i) + indexCoeff(j)*(binMax(j) - binMin(j))
         end do
      end do

      umbrellaLimit = 1
      do i = 1, nBiasVariables
         umbrellaLimit = umbrellaLimit + indexCoeff(i)*(int(binMax(i), 4) - int(binMin(i), 4))
      end do

      write (nout, *) "Number of Umbrella Bins:", umbrellaLimit

      allocate (UBias(1:umbrellaLimit + 1), STAT=AllocateStatus)
      allocate (UHist(1:umbrellaLimit + 1), STAT=AllocateStatus)

      UBias = 0E0_dp
      UHist = 0E0_dp

      if (energyAnalytics) then
         allocate (U_EAvg(1:umbrellaLimit + 1), STAT=AllocateStatus)
         allocate (U_EHist(1:umbrellaLimit + 1, 0:E_Bins), STAT=AllocateStatus)
         allocate (UHistTotal(1:umbrellaLimit + 1), STAT=AllocateStatus)
         U_EAvg = 0E0_dp
         U_EHist = 0E0_dp
         UHistTotal = 0E0_dp
      end if

      IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
   end subroutine
!==========================================================================================
   subroutine ReadInitialBias
      implicit none
      integer :: AllocateStatus
      integer :: j, iInput, iBias, inStat, biasIndx
      real(dp), allocatable :: varValue(:)
      real(dp) :: curBias

!    open(unit=80, file="TempIn.txt")
      open (unit=80, file=trim(adjustl(inputFile)))
      allocate (varValue(1:nBiasVariables), STAT=AllocateStatus)
      IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"

      UBias = 0E0_dp
      do iInput = 1, nint(1d7)
         read (80, *, IOSTAT=inStat) (varValue(j), j=1, nBiasVariables), curBias

         if (inStat .lt. 0) then
            exit
         end if

         call getUIndexArray(varValue, biasIndx, inStat)
!      write(*, *) (varValue(j), j=1,nBiasVariables), curBias
!      write(*,*) varValue, biasIndx
         if (inStat .eq. 1) then
            cycle
         end if

         UBias(biasIndx) = curBias
      end do

      deallocate (varValue)

      close (80)

   end subroutine
!==========================================================================================
   subroutine UmbrellaHistAdd(E_T)
      implicit none
      real(dp), intent(in) :: E_T

      integer :: bin

      curUIndx = getBiasIndex()
      UHist(curUIndx) = UHist(curUIndx) + 1E0_dp
      ! write(*,*) curUIndx, UHist(curUIndx)

      if (energyAnalytics) then
         U_EAvg(curUIndx) = U_EAvg(curUIndx) + E_T
         UHistTotal(curUIndx) = UHistTotal(curUIndx) + 1E0_dp
         bin = floor(E_T/dE)

         if ((bin .ge. 0) .and. (bin .lt. E_Bins)) then
            U_EHist(curUIndx, bin) = U_EHist(curUIndx, bin) + 1d0
         else
            U_EHist(curUIndx, E_Bins) = U_EHist(curUIndx, E_Bins) + 1d0
         end if
         ! write(*,*) curUIndx, UHist(curUIndx), U_EHist(curUIndx, bin)
      end if

   end subroutine
!==========================================================================================
   subroutine GetUmbrellaBias_Disp(disp, biasDiff, rejMove)
      use CoordinateTypes
      use SimParameters, only: NPART, NPART_new, NTotal, NTotal_New, temperature, TempNew
      implicit none
      logical, intent(out) :: rejMove
      type(Displacement), intent(in) :: disp(:)
      real(dp), intent(out) :: biasDiff
      integer :: iDispFunc, newUIndx, sizeDisp
      real(dp) :: biasOld, biasNew

      if (.not. useUmbrella) then
         biasDiff = 0E0_dp
         return
      end if

      NPART_New = NPART
      NTotal_New = NTotal
      TempNew = temperature
      rejMove = .false.
      curUIndx = getBiasIndex()
      biasOld = UBias(curUIndx)
      if (nDispFunc .ne. 0) then
         sizeDisp = size(disp)
         do iDispFunc = 1, nDispFunc
            call DispUmbrella(iDispFunc)%func(disp(1:sizeDisp))
         end do
      end if

      call getNewBiasIndex(newUIndx, rejMove)
      if (rejMove) then
         write (*, *) "BLAH!"
         return
      end if
      biasNew = UBias(newUIndx)
      biasDiff = biasNew - biasOld
   end subroutine

!==========================================================================================
   subroutine GetUmbrellaBias_SwapIn(biasDiff, rejMove)
      use SimParameters, only: NPART, NPART_new, NTotal, NTotal_New, temperature, TempNew
      implicit none
      logical, intent(out) :: rejMove
      real(dp), intent(out) :: biasDiff
      integer :: iSwapFunc, newUIndx, sizeDisp
      real(dp) :: biasOld, biasNew

      if (.not. useUmbrella) then
         biasDiff = 0E0_dp
         return
      end if

      TempNew = temperature
      rejMove = .false.
      curUIndx = getBiasIndex()
      biasOld = UBias(curUIndx)
      if (nSwapInFunc .ne. 0) then
         do iSwapFunc = 1, nSwapInFunc
            call SwapInUmbrella(iSwapFunc)%func
         end do
      end if

      call getNewBiasIndex(newUIndx, rejMove)
      if (rejMove) then
         return
      end if
      biasNew = UBias(newUIndx)
      biasDiff = biasNew - biasOld

!     write(2,*) biasNew, biasOld, biasDiff

   end subroutine
!==========================================================================================
   subroutine GetUmbrellaBias_SwapOut(nType, nMol, biasDiff, rejMove)
      use SimParameters, only: NPART, NPART_new, NTotal, NTotal_New
      implicit none
      logical, intent(out) :: rejMove
      integer, intent(in) :: nType, nMol
      real(dp), intent(out) :: biasDiff
      integer :: iSwapFunc, newUIndx, sizeDisp
      real(dp) :: biasOld, biasNew

      if (.not. useUmbrella) then
         biasDiff = 0E0_dp
         return
      end if

      rejMove = .false.
      curUIndx = getBiasIndex()
      biasOld = UBias(curUIndx)
      if (nSwapOutFunc .ne. 0) then
         do iSwapFunc = 1, nSwapInFunc
            call SwapOutUmbrella(iSwapFunc)%func(nType, nMol)
         end do
      end if

      call getNewBiasIndex(newUIndx, rejMove)
      if (rejMove) then
         return
      end if
      biasNew = UBias(newUIndx)
      biasDiff = biasNew - biasOld

   end subroutine
!==========================================================================================
   subroutine GetUmbrellaBias_Temperature(biasDiff, rejMove)
      use SimParameters, only: NPART, NPART_new, NTotal, NTotal_New, temperature, TempNew
      implicit none
      logical, intent(out) :: rejMove
      real(dp), intent(out) :: biasDiff
      integer :: iUmbrella, newUIndx
      real(dp) :: biasNew, biasOld

      if (.not. useUmbrella) then
         biasDiff = 0E0_dp
         return
      end if

      rejMove = .false.
      curUIndx = getBiasIndex()
      biasOld = UBias(curUIndx)

      NPART_New = NPART
      NTotal_New = NTotal

!     do iUmbrella = 1, nBiasVariables
!       if(biasvarnew(iUmbrella) % varType .eq. 1) then
!         biasvarnew(iUmbrella) % intVar = biasvar(iUmbrella) % intVar
!       elseif(biasvarnew(iUmbrella) % varType .eq. 2) then
!         biasvarnew(iUmbrella) % realVar = biasvar(iUmbrella) % realVar
!       endif
!     enddo

      call getNewBiasIndex(newUIndx, rejMove)
      if (rejMove) then
         return
      end if
      biasNew = UBias(newUIndx)
      biasDiff = biasNew - biasOld

   end subroutine
!==========================================================================================
   subroutine ScreenOutputUmbrella
      use ParallelVar
      implicit none
      integer :: iBias, j
      ! print *, "This the screen output umbrella working!"
      do iBias = 1, nBiasVariables
         if (biasvar(iBias)%varType .eq. 1) then
            varValues(iBias) = real(biasvar(iBias)%intVar, dp)
         elseif (biasvar(iBias)%varType .eq. 2) then
            varValues(iBias) = biasvar(iBias)%realVar
         end if
      end do
      ! print *, "nBiasVariable: ",nBiasVariables, varValues
      write (nout, screenFormat) (varValues(j), j=1, nBiasVariables)

   end subroutine
!==========================================================================================
   subroutine OutputUmbrellaHist
      implicit none
      integer :: iUmbrella, iBias, iBin
!    integer, allocatable :: UArray(:)
      character(len=100) :: outputString

!    allocate(UArray(1:nBiasVariables))
!    allocate(varValues(1:nBiasVariables))

      write (outputString, *) "(", (trim(outputFormat(iBias)), iBias=1, nBiasVariables), "2x, F18.1)"
      open (unit=60, file="TemporaryHist.txt")
      do iUmbrella = 1, umbrellaLimit
         call findVarValues(iUmbrella, UArray)
         do iBias = 1, nBiasVariables
            varValues(iBias) = real(UArray(iBias), dp)*UBinSize(iBias)
         end do
         if (UHist(iUmbrella) .ne. 0E0_dp) then
            write (60, outputString) (varValues(iBias), iBias=1, nBiasVariables), UHist(iUmbrella)
         end if
      end do
      flush (60)
      close (60)

!    deallocate(UArray)
!    deallocate(varValues)

   end subroutine

!==========================================================================================
   subroutine OutputUmbrellaAnalytics
      use ParallelVar
!    use MPI
      implicit none
      include 'mpif.h'
      integer :: iUmbrella, iBias, iBin
      integer :: arraySize
      real(dp), allocatable :: TempHist(:), TempAvg(:), Temp2D(:, :)
      character(len=100) :: outputString

      if (.not. useUmbrella) then
         return
      end if

!    if(myid .eq. 0) then
      allocate (TempHist(1:umbrellaLimit + 1))
      allocate (TempAvg(1:umbrellaLimit + 1))
      allocate (Temp2D(1:umbrellaLimit + 1, 0:E_Bins))
      TempHist = 0E0_dp
      TempAvg = 0E0_dp
      Temp2D = 0E0_dp
!    endif
      call MPI_BARRIER(MPI_COMM_WORLD, ierror)
      arraySize = size(U_EAvg)
!    write(*,*) arraySize, size(TempHist)
      call MPI_REDUCE(U_EAvg, TempAvg, arraySize, &
                      MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierror)

      call MPI_REDUCE(UHistTotal, TempHist, arraySize, &
                      MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierror)
!    write(*,*) myid, ierror
!    do iUmbrella = 1, umbrellaLimit
!      write(*,*) iUmbrella, U_EAvg(iUmbrella), TempHist(iUmbrella)
!    enddo

      if (myid .eq. 0) then
         do iUmbrella = 1, umbrellaLimit
            U_EAvg(iUmbrella) = TempAvg(iUmbrella)
            UHistTotal(iUmbrella) = TempHist(iUmbrella)
         end do
      end if

      call MPI_BARRIER(MPI_COMM_WORLD, ierror)
      arraySize = size(U_EHist)
      call MPI_REDUCE(U_EHist, Temp2D, arraySize, &
                      MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierror)
!    write(*,*) myid, ierror

      if (myid .eq. 0) then
         do iUmbrella = 1, umbrellaLimit
            do iBin = 0, E_Bins
               U_EHist(iUmbrella, iBin) = Temp2D(iUmbrella, iBin)
            end do
         end do
      end if

      if (myid .eq. 0) then
         open (unit=60, file="Umbrella_AvgE.txt")
         do iUmbrella = 1, umbrellaLimit
            if (UHistTotal(iUmbrella) .ne. 0E0_dp) then
               call findVarValues(iUmbrella, UArray)
               do iBias = 1, nBiasVariables
                  varValues(iBias) = real(UArray(iBias), dp)*UBinSize(iBias)
               end do
               write (60, *) (varValues(iBias), iBias=1, nBiasVariables), U_EAvg(iUmbrella)/UHistTotal(iUmbrella)
            end if
         end do
         close (60)

         open (unit=60, file="Umbrella_DensityStates.txt")
         do iUmbrella = 1, umbrellaLimit
            if (UHistTotal(iUmbrella) .ne. 0E0_dp) then
               call findVarValues(iUmbrella, UArray)
               do iBias = 1, nBiasVariables
                  varValues(iBias) = real(UArray(iBias), dp)*UBinSize(iBias)
               end do
               write (60, *) (varValues(iBias), iBias=1, nBiasVariables)
               do iBin = 0, E_Bins
                  if (U_EHist(iUmbrella, iBin) .ne. 0d0) then
                     write (60, *) iBin*dE, U_EHist(iUmbrella, iBin)
                  end if
               end do
            end if
         end do
         close (60)
      end if

      deallocate (TempHist)
      deallocate (TempAvg)
      deallocate (Temp2D)

   end subroutine
!==========================================================================
   function getBiasIndex() result(biasIndx)
      integer :: biasIndx
      integer :: iBias

      !  print *, "WORKING IN THE getBiasIndex()", nBiasVariables
      do iBias = 1, nBiasVariables
         if (biasvar(iBias)%varType .eq. 1) then
            !  print *, "intVar", biasvar(iBias) % intVar
            binIndx(iBias) = biasvar(iBias)%intVar
         elseif (biasvar(iBias)%varType .eq. 2) then
            binIndx(iBias) = floor(biasvar(iBias)%realVar/UBinSize(iBias))
         end if
      end do

      biasIndx = 1
      do iBias = 1, nBiasVariables
!       write(*,*) biasIndx,indexCoeff(iBias), binIndx(iBias), binMin(iBias)
         biasIndx = biasIndx + indexCoeff(iBias)*(binIndx(iBias) - binMin(iBias))
      end do

!     write(*,*) biasIndx

   end function
!==========================================================================
   subroutine getNewBiasIndex(biasIndx, rejMove)
      logical, intent(out) :: rejMove
      integer, intent(out) :: biasIndx
      integer :: iBias

      rejMove = .false.
      do iBias = 1, nBiasVariables
         if (biasvarnew(iBias)%varType .eq. 1) then
            binIndx(iBias) = biasvarnew(iBias)%intVar
         elseif (biasvar(iBias)%varType .eq. 2) then
            binIndx(iBias) = floor(biasvarnew(iBias)%realVar/UBinSize(iBias))
         end if
         if (binIndx(iBias) .lt. binMin(iBias)) then
            rejMove = .true.
            return
         end if
         if (binIndx(iBias) .gt. binMax(iBias)) then
            rejMove = .true.
            return
         end if

      end do

      biasIndx = 1
      do iBias = 1, nBiasVariables
         biasIndx = biasIndx + indexCoeff(iBias)*(binIndx(iBias) - binMin(iBias))
      end do
      if (biasIndx .gt. umbrellaLimit) then
         rejMove = .true.
         return
      end if

   end subroutine
!==========================================================================
   subroutine getUIndexArray(varArray, biasIndx, stat)
      real(dp), intent(in) :: varArray(:)
      integer, intent(out) :: biasIndx, stat
      integer :: iBias

      stat = 0
      do iBias = 1, nBiasVariables
!       binIndx(iBias) = floor( varArray(iBias) / UBinSize(iBias) + 1E-8 )
         binIndx(iBias) = nint(varArray(iBias)/UBinSize(iBias))
         if (binIndx(iBias) .gt. binMax(iBias)) then
            stat = 1
            return
         end if
         if (binIndx(iBias) .lt. binMin(iBias)) then
            stat = -1
            return
         end if
      end do

      biasIndx = 1
      do iBias = 1, nBiasVariables
         biasIndx = biasIndx + indexCoeff(iBias)*(binIndx(iBias) - binMin(iBias))
      end do

   end subroutine
!===========================================================================
   subroutine findVarValues(UIndx, UArray)
      implicit none
      integer, intent(in) :: UIndx
      integer, intent(inout) :: UArray(:)
      integer :: i, iBias
      integer :: remainder, curVal

      remainder = UIndx - 1
      do i = 1, nBiasVariables
         iBias = nBiasVariables - i + 1
         curVal = int(real(remainder, dp)/real(indexCoeff(iBias), dp))
         UArray(iBias) = curVal + binMin(iBias)
         remainder = remainder - curVal*indexCoeff(iBias)
      end do

   end subroutine
!==========================================================================================
   subroutine CheckInitialValues
      implicit none
      integer :: iBias

      do iBias = 1, nBiasVariables
         if (biasvar(iBias)%varType .eq. 1) then
            binIndx(iBias) = nint(biasvar(iBias)%intVar/UBinSize(iBias))
         elseif (biasvar(iBias)%varType .eq. 2) then
            binIndx(iBias) = floor(biasvar(iBias)%realVar/UBinSize(iBias))
         end if

         if (binIndx(iBias) .gt. binMax(iBias)) then
            write (*, *) "The initital system state is above the upper bounds"
            write (*, *) "specified by the umbrella sampling input"
            write (*, *) "Umbrella Variable:", iBias
            write (*, *) "Bin Index:", binIndx(iBias)
            write (*, *) "Largest allowed bin:", binMax(iBias)
            stop
         end if
         if (binIndx(iBias) .lt. binMin(iBias)) then
            write (*, *) "The initital system state is below the lower bounds"
            write (*, *) "specified by the umbrella sampling input"
            write (*, *) "Umbrella Variable:", iBias
            write (*, *) "Bin Index:", binIndx(iBias)
            write (*, *) "Smallest allowed bin:", binMin(iBias)
            stop
         end if
      end do

   end subroutine
!==========================================================================
end module
!==========================================================================

