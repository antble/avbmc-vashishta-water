!****************************************************************************************

!****************************************************************************************
module E_Interface_Vashishta
   use CoordinateTypes
!=============================================================================
contains
!=============================================================================
   subroutine Detailed_EnergyCalc_Vashishta(E_T, rejMove)
      use InterEnergy_Vashishta, only: Detailed_ECalc_Inter
      use DistanceCriteria, only: Detailed_DistanceCriteria
      use EnergyCriteria, only: Detailed_EnergyCriteria

      use SimParameters, only: maxMol, distCriteria

      implicit none

      logical, intent(inout) :: rejMove
      real(dp), intent(inout) :: E_T
      
      integer :: i, j
      real(dp) :: PairList(1:maxMol, 1:maxMol)

      E_T = 0E0
      
      call Detailed_ECalc_Inter(E_T, PairList)
      if (distCriteria) then
         call Detailed_DistanceCriteria(PairList, rejMove)
      else
         call Detailed_EnergyCriteria(PairList, rejMove)
      end if

   end subroutine
!=============================================================================
!     This function contains the energy and cluster criteria functions for any move
!     where a molecule is moved within a cluster.  This function takes a set of displacement
!     vectors (the "disp" variable) and returns the change in energy for both the Intra- and
!     Inter-molecular components.  This function can be used for for moves that any number of
!     atoms in a given molecule, but it can not be used if more than one molecule changes
!     in a given move.
   subroutine Shift_EnergyCalc_Vashishta(E_Inter, E_Intra, disp, PairList, dETable, useIntra, rejMove, update)
      use SimParameters, only: distCriteria, beta, softcutoff, NTotal, outputEConv 
      use CBMC_Variables
      use Coords
      use DistanceCriteria
      use EnergyCriteria
      use EnergyTables
      use InterEnergy_Vashishta
      implicit none

      logical, intent(in), optional :: update ! useInter, 
      logical, intent(in) :: useIntra(1:4)
      type(Displacement), intent(in) :: disp(:)
      logical, intent(inout) :: rejMove
      real(dp), intent(inout) :: dETable(:)
      real(dp), intent(out) :: E_Intra, E_Inter
      real(dp), intent(InOut) :: PairList(:)

      logical :: interSwitch
      integer :: nIndx, nDisp
      real(dp) :: E_NonBond, E_Stretch, E_Bend
      real(dp) :: E_Torsion, E_Improper

      nDisp = size(disp)
      rejMove = .false.
      E_Inter = 0E0
      E_Intra = 0E0
      E_NonBond = 0E0
      dETable = 0E0

      E_Inter_Diff = 0E0
      E_NBond_Diff = 0E0

      ! if (present(useInter)) then
      !    interSwitch = useInter
      ! else
      !    interSwitch = .true.
      ! end if

!     Begin by calculating the intermolecular potential. If any atoms overlap the move will be rejected
!     immediately. If the calling function has specfied that useInter is .false. then this part will be
!     skipped and only the intra molecular potential will be calculated.

      dETable = 0E0 
      PairList = 0E0 
      call Shift_ECalc_Inter(E_Inter, disp, PairList, dETable, rejMove, update)

      if (rejMove) then
         return
      end if
      E_Inter_Diff = E_Inter

      ! if (E_Inter*beta .gt. softCutOff) then
      !    print *, "SOFTCUT-OFF", E_Inter*beta
      !    rejMove = .true.
      !    stop
      !    return
      ! end if

!        Using the data collected from the intermolecular function, check to see that the new position
!        satisfies the cluster criteria.
      nIndx = MolArray(disp(1)%molType)%mol(disp(1)%molIndx)%indx

      if (any(disp(1:nDisp)%atmIndx .eq. 1)) then
         call Shift_DistanceCriteria(PairList, nIndx, rejMove)
      end if

      if (rejMove) then
         return
      end if
   end subroutine
!=============================================================================
   subroutine SwapIn_EnergyCalc_Vashishta(E_Inter, E_Intra, PairList, dETable, rejMove, useInter)
      use InterEnergy_Vashishta
      use SimParameters, only: outputEConv 
      use EnergyCriteria
      use EnergyTables
      use Coords
      use CBMC_Variables
      implicit none

      logical, intent(out) :: rejMove
      logical, intent(in), optional :: useInter
      real(dp), intent(out) :: E_Inter, E_Intra
      real(dp), intent(inout) :: PairList(:), dETable(:)

      logical :: interSwitch
      real(dp) :: E_NonBond, E_Stretch, E_Bend
      real(dp) :: E_Torsion, E_Improper

      rejMove = .false.

      E_Inter = 0E0
      E_Intra = 0E0
      E_NonBond = 0E0
      PairList = 0E0
      dETable = 0E0

      E_Inter_Diff = 0E0
      E_NBond_Diff = 0E0

      call NewMol_ECalc_Inter(E_Inter, PairList, dETable, rejMove, E_Intra)
      if (rejMove) then
         return
      end if

      E_Inter_Diff = E_Inter
   end subroutine
!=============================================================================
!     This function contains the energy calculations that are used when a molecule
!     has been selected for removal.
   subroutine SwapOut_EnergyCalc_Vashishta(E_Inter, E_Intra, nType, nMol, dETable, useInter)
      use InterEnergy_Vashishta
      use SimParameters, only: outputEConv 
      use EnergyCriteria
      use EnergyTables
      use Coords
      use CBMC_Variables
      implicit none

      logical, intent(in), optional :: useInter
      real(dp), intent(out) :: E_Inter, E_Intra
      integer, intent(in) :: nType, nMol
      real(dp), intent(inout) :: dETable(:)

      logical :: interSwitch
      real(dp) :: E_NonBond, E_Stretch, E_Bend
      real(dp) :: E_Torsion, E_Improper

      E_Inter = 0E0
      E_Intra = 0E0

      E_Inter_Diff = 0E0
      dETable = 0E0

      call Mol_ECalc_Inter(nType, nMol, dETable, E_Inter, E_Intra)
      E_Inter = -E_Inter 
      E_Inter_Diff = E_Inter
   end subroutine

!=============================================================================
!     This function contains the energy calculations that are used when a molecule
!     has been selected for removal.
   subroutine Update_SubEnergies_Vashishta
      use EnergyTables
      implicit none

      E_Inter_T = E_Inter_T + E_Inter_Diff
      E_Inter_Diff = 0E0
   end subroutine
!=============================================================================

end module

