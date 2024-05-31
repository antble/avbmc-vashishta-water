
!======================================================
module CoordinateTypes
   use VarPrecision

   type FloatPointer
      real(dp), pointer :: pnt
   end type

   type SimpleMolCoords
      real(dp), allocatable :: x(:), y(:), z(:)
   end type

   type SimpleAtomCoords
      real(dp) :: x, y, z
   end type

   type GlobalAtomIndex
      integer(kind=atomIntType) :: nType, nMol, nAtom, atmType
   end type

   type Molecule
      integer(kind=atomIntType) :: indx
      integer, allocatable :: globalIndx(:)
      real(dp), allocatable :: x(:), y(:), z(:)
   end type

   type MolArrayType
      type(Molecule), allocatable :: mol(:)
   end type

   type MolPointers
      integer(kind=atomIntType) :: molType, molIndx
      type(FloatPointer), allocatable :: x(:), y(:), z(:)
   end type

   type Displacement
      integer(kind=atomIntType) :: molType, atmIndx, molIndx
      real(dp):: x_new, y_new, z_new
      real(dp), pointer :: x_old, y_old, z_old
      ! new data structure for the single atom translation
      real(dp), pointer :: x1_new, y1_new, z1_new
      real(dp) :: x1_old, y1_old, z1_old
   end type

   type TrialCoordinates
      integer(kind=atomIntType) :: molType, molIndx
      real(dp), allocatable :: x(:), y(:), z(:)
   end type

end module

!==============================================================
module AcceptRates
   use VarPrecision
   real(dp), allocatable :: acptTrans(:), acptRot(:)
   real(dp), allocatable :: atmpTrans(:), atmpRot(:)

   real(dp), allocatable :: acptSwapIn(:), atmpSwapIn(:)
   real(dp), allocatable :: acptSwapOut(:), atmpSwapOut(:)

   real(dp), allocatable :: acptInSize(:), atmpInSize(:)

   real(dp) :: angGen_accpt, angGen_atmp
   real(dp) :: dihedGen_accpt, dihedGen_atmp
   real(dp) :: distGen_accpt, distGen_atmp
   real(dp) :: clusterCritRej

   real(dp) :: reversibility_count 
end module
!==============================================================
module AVBMC_RejectionVar
   use CoordinateTypes

   real(dp) :: totalRej
   real(dp) :: ovrlapRej, dbalRej, critriaRej, boundaryRej

   real(dp) :: totalRej_out
   real(dp) :: boundaryRej_out, dbalRej_out, critriaRej_out

   real(dp) :: OHrej_out

end module
!======================================================
module ParallelVar
   integer :: myid, p_size, ierror, tag, nout
end module
!==============================================================
module Coords
   use CoordinateTypes

   type(MolArrayType), allocatable, target :: MolArray(:)
!      type(MolPointers), allocatable :: JointArray(:)
   type(GlobalAtomIndex), allocatable :: atomIndicies(:)

   type(TrialCoordinates) :: newMol, newMol2
   type(SimpleMolCoords), allocatable :: rosenTrial(:)

   logical, allocatable :: NeighborList(:, :)
   integer, allocatable :: NumNei(:)
   integer, allocatable :: typeList(:)
   integer, allocatable :: subIndxList(:)

   ! real(dp), allocatable :: PairList(:)

   type(SimpleMolCoords), allocatable :: gasConfig(:)

end module

!==============================================================
module CBMC_Variables
   use Constants
   use VarPrecision

   type BondNumber
      integer(kind=atomIntType), allocatable :: atom(:)
   end type

   type Pathing
      integer :: nTerminal, nLinker, nHub
      integer(kind=atomIntType), allocatable :: termAtoms(:)
      integer(kind=atomIntType), allocatable :: hubAtoms(:)
      integer :: nPaths
      integer(kind=atomIntType), allocatable :: path(:, :)
      integer(kind=atomIntType), allocatable :: pathMax(:)
   end type

   integer, parameter :: nDihBins = 2000
   real(dp), parameter :: diBinSize = two_pi/real(nDihBins, dp)
   real(dp) :: startProb = 0.05E0
   type DihedralAngle
      integer :: molType, hubIndx, dihedIndx
      integer :: startBin
      real(dp) :: accConst
      real(dp) :: Hist(0:nDihBins)
      real(dp) :: Prob(0:nDihBins)
      real(dp) :: Integral(0:nDihBins)
   end type

   integer :: maxRosenTrial
   integer, allocatable :: nRosenTrials(:)
   integer, allocatable :: regrowType(:)
   integer, allocatable :: regrowOrder(:, :)
   real(dp), allocatable :: probTypeCBMC(:)
   type(BondNumber), allocatable :: topolArray(:)
   type(Pathing), allocatable :: pathArray(:)

   integer, allocatable :: usedByPath(:, :)
   integer, allocatable :: atomPathIndex(:, :)

   integer :: totalDihed
   type(DihedralAngle), allocatable :: dihedData(:)

   integer, parameter :: maxBranches = 6
   integer, parameter :: nRosenTwoBranch = 5
   integer, parameter :: nRosenThreeBranch = 1
   integer, parameter :: nRosenTorsion = 100

   type GrowthSchedule
      integer :: GrowthSteps
      integer, allocatable :: GrowFrom(:), GrowPrev(:), GrowNum(:), TorNum(:)
      integer, allocatable :: GrowList(:, :), TorList(:, :)
   end type

   type(GrowthSchedule), allocatable :: SwapGrowOrder(:)

end module

!==============================================================

!      module CBMC_Pointers
!
!      interface
!        subroutine growFunc(nType)
!          integer, intent(in) :: nType
!        end subroutine
!      end interface
!
!      type GrowPointer
!        procedure(growFunc), pointer :: pnt
!      end type
!
!      type(GrowPointer), allocatable :: growthFunc(:)
!
!      end module
!==============================================================
module EnergyTables
   use VarPrecision

   integer, allocatable :: neiCount(:)
   real(dp), allocatable :: ETable(:)
   real(dp), allocatable :: NeiETable(:)
   real(dp), allocatable :: biasAlpha(:, :)

   real(dp) :: E_Inter_T, E_Inter_Diff
   real(dp) :: E_NBond_T, E_NBond_Diff
   real(dp) :: E_Stretch_T, E_Strch_Diff
   real(dp) :: E_Bend_T, E_Bend_Diff
   real(dp) :: E_Torsion_T, E_Tors_Diff
end module

!================================================================
module IndexingFunctions
contains

!     ----------------------------------------------------
   subroutine Get_MolIndex(nMove, NPart, molType, molIndx)
      implicit none
      integer, intent(in) :: nMove, NPart(:)
      integer, intent(inout) :: molIndx, molType
      integer :: i, sizeN, curLimit

      sizeN = size(NPART)
!       curLimit = NPART(1)
      curLimit = 0
      molIndx = 1

      molType = 1
      do i = 1, sizeN
         curLimit = curLimit + NPART(i)
         if (nMove .le. curLimit) then
            molIndx = nMove - curLimit + NPART(i)
            molType = i
            return
         end if
      end do

      write (35, *) "NotFound", nMove, molType, molIndx
   end subroutine
!     ----------------------------------------------------
   pure subroutine Get_SubIndex(nIndx, nType, NMAX)
      implicit none
      integer, intent(in) :: nIndx, NMAX(:)
      integer, intent(out) :: nType
      integer :: i, sizeN, curLimit

      sizeN = size(NMAX)
      curLimit = 0
      do i = 1, nType - 1
         curLimit = curLimit + NMAX(i)
      end do
      nType = nIndx - curLimit

   end subroutine

!     ----------------------------------------------------
   pure integer function Get_MolType(Indx, NMAX)
      implicit none
      integer, intent(in) :: Indx
      integer, intent(in) :: NMAX(:)
      integer :: i, sizeN, curLimit

      sizeN = size(NMAX)
      curLimit = 0
      Get_MolType = 0
      do i = 1, sizeN
         curLimit = curLimit + NMAX(i)
         if (Indx .le. curLimit) then
            Get_MolType = i
            return
         end if
      end do

   end function
!     ----------------------------------------------------

end module

!==============================================================
module SimParameters
   use Constants, only : pi
   use VarPrecision
   
   ! logical, parameter :: echoInput = .false.
   logical, parameter :: echoInput = .true.
   logical, parameter :: useScriptInput = .true.
   logical :: prevMoveAccepted

   logical :: distCriteria = .false.
   logical :: multipleInput = .false.
   logical :: calcPressure = .false.
   integer(kind=8) :: ncycle, ncycle2
   integer(kind=atomIntType) :: nMolTypes = 1
   integer, allocatable :: NMin(:), NMax(:)
   integer, allocatable, target :: NPart(:)
   integer, allocatable, target :: NPart_New(:)
   logical, allocatable :: isActive(:)
   integer, target :: NTotal, NTotal_New
   integer :: maxMol
   integer :: maxAtoms, vmdAtoms
   real(dp) :: avbmc_vol

   integer :: umbrellaLimit
   real(dp), allocatable :: NHist(:), NBias(:), E_Avg(:)
   real(dp), allocatable :: gas_dens(:)
   real(dp), allocatable :: max_dist(:), max_dist_single(:)
   real(dp), allocatable :: max_rot(:)

   real(dp) :: beta
   real(dp), target :: temperature, tempNew
   real(dp) :: global_r_min
   real(dp) :: softCutoff
   real(dp) :: Pressure, P_Diff
   real(dp), allocatable :: P_Avg(:)
   ! Vashishta parameter required 
   real(dp) :: r1s 
   real(dp) :: rc 

   ! Gas-phase geometry of water
   ! real(dp) :: rb_max = 1.4  ! angstrom
   ! real(dp) :: rb_min = 0.5 
   ! real(dp) :: theta_max = pi*114.52/180 ! radian
   ! real(dp) :: theta_min = pi*94.52/180 

!! ORIGINAL PARAMETERS USED IN TRAINING THE WATER MODEL
   ! real(dp) :: rb_max = 1.2 !0.9572 
   ! real(dp) :: rb_min = 0.9572 
   ! real(dp) :: theta_max = pi*104.52/180
   ! real(dp) :: theta_min = pi*97/180!pi*104.52/180

! REVISED PARAMETERS | Gas-phase distribution 
   real(dp) :: rb_max = 1.10  !0.9572 
   real(dp) :: rb_min = 0.92
   real(dp) :: theta_max = pi*100/180
   real(dp) :: theta_min = pi*92/180
   

!     Cluster Criteria block
   real(dp) :: Dist_Critr
   real(dp) :: Dist_Critr_sq
!      real(dp), allocatable :: Dist_Critr(:,:)
!      real(dp), allocatable :: Dist_Critr_sq(:,:)
   real(dp) :: NeighRej
   real(dp), allocatable :: Eng_Critr(:, :)
   real(dp) :: ECritMax

   character(len=10) :: outputEngUnits, outputLenUnits
   real(dp) :: outputEConv, outputLenConv

   integer :: outFreq_GCD = 1000
   integer :: outFreq_Screen = 1000
   integer :: outFreq_Traj = 1000

end module
!================================================================
module UmbrellaFunctions

contains

   pure integer function getBiasIndex(NPart, NMAX)
      implicit none
      integer, intent(in) :: NPart(:), NMAX(:)
      integer :: i, sizeN
      integer :: curIndx, maxIndx

      sizeN = size(NPART)
!       curIndx = 0
      curIndx = NPart(sizeN)
      maxIndx = NMAX(sizeN) + 1
      do i = 1, sizeN - 1
         curIndx = curIndx + maxIndx*NPart(sizeN - i)
         maxIndx = maxIndx*(NMAX(sizeN - i) + 1)
      end do
      getBiasIndex = curIndx + 1

   end function
!     ---------------------------------------------------------------
   pure integer function getNewBiasIndex(NPart, NMAX, increment)
      implicit none
      integer, intent(in) :: NPart(:), NMAX(:), increment(:)
      integer :: i, sizeN
      integer :: curIndx, maxIndx

      sizeN = size(NPART)
!       curIndx = 0
      curIndx = NPart(sizeN) + increment(sizeN)
      maxIndx = NMAX(sizeN) + 1
      do i = 1, sizeN - 1
         curIndx = curIndx + maxIndx*(NPart(sizeN - i) + increment(sizeN - i))
         maxIndx = maxIndx*(NMAX(sizeN - i) + 1)
      end do
      getNewBiasIndex = curIndx + 1

   end function

!     ---------------------------------------------------------------
end module
!======================================================
module WHAM_Module
   use VarPrecision

   logical, parameter :: WHAM_ExtensiveOutput = .false.

!      integer, allocatable :: refSizeNumbers(:)
   real(dp), allocatable :: refSizeNumbers(:)
   integer :: refBin
   logical :: useWHAM = .false.
   integer :: intervalWHAM
   integer :: maxSelfConsist
   real(dp) :: tolLimit
   integer :: equilInterval
   integer :: whamEstInterval
!      real(dp) :: unsampledIncrmnt = 1.5d0

   integer :: nWhamItter, nCurWhamItter

   real(dp), allocatable :: WHAM_Numerator(:)
   real(dp), allocatable :: WHAM_Denominator(:, :)
   real(dp), allocatable :: HistStorage(:)
   real(dp), allocatable :: FreeEnergyEst(:)
   real(dp), allocatable :: BiasStorage(:, :)
   real(dp), allocatable :: NewBias(:)
   real(dp), allocatable :: ProbArray(:)
   real(dp), allocatable :: TempHist(:)

end module


module DetailedBalance

end module 
!================================================================
