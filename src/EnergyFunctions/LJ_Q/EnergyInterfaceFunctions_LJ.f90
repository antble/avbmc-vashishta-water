!****************************************************************************************

!****************************************************************************************
module E_Interface_LJ_Q
   use CoordinateTypes
!=============================================================================
contains
!=============================================================================
   subroutine Detailed_EnergyCalc_LJ_Q(E_T, rejMove)
      use InterEnergy_LJ_Electro, only: Detailed_ECalc_Inter
      use IntraEnergy_LJ_Electro, only: Detailed_ECalc_IntraNonBonded
      use BondStretchFunctions, only: Detailed_ECalc_BondStretch
      use BendingFunctions, only: Detailed_ECalc_Bending
      use TorsionalFunctions, only: Detailed_ECalc_Torsional
      use ImproperAngleFunctions

      use DistanceCriteria, only: Detailed_DistanceCriteria
      use EnergyCriteria, only: Detailed_EnergyCriteria

      use SimParameters, only: maxMol, distCriteria
      use E_Interface_LJ_Q_Diststore, only: Detailed_EnergyCalc_LJ_Q_DStore
      use PairStorage, only: useDistStore

      implicit none

      logical, intent(inout) :: rejMove
      real(dp), intent(inout) :: E_T
      integer :: i, j
      real(dp) :: PairList(1:maxMol, 1:maxMol)

      if (useDistStore) then
         call Detailed_EnergyCalc_LJ_Q_DStore(E_T, rejMove)
         return
      end if

      E_T = 0E0
      call Detailed_ECalc_Inter(E_T, PairList)

      if (distCriteria) then
         call Detailed_DistanceCriteria(PairList, rejMove)
      else
         call Detailed_EnergyCriteria(PairList, rejMove)
      end if

      call Detailed_ECalc_IntraNonBonded(E_T)
      call Detailed_ECalc_BondStretch(E_T)
      call Detailed_ECalc_Bending(E_T)
      call Detailed_ECalc_Torsional(E_T)
!      call Detailed_ECalc_Improper(E_T)

   end subroutine
!=============================================================================
!     This function contains the energy and cluster criteria functions for any move
!     where a molecule is moved within a cluster.  This function takes a set of displacement
!     vectors (the "disp" variable) and returns the change in energy for both the Intra- and
!     Inter-molecular components.  This function can be used for for moves that any number of
!     atoms in a given molecule, but it can not be used if more than one molecule changes
!     in a given move.
   subroutine Shift_EnergyCalc_LJ_Q(E_Inter, E_Intra, disp, PairList, dETable, useIntra, rejMove, useInter)
      use SimParameters, only: distCriteria, beta, softcutoff, NTotal
      use BendingFunctions
      use BondStretchFunctions
      use CBMC_Variables
      use Coords
      use DistanceCriteria
      use EnergyCriteria
      use EnergyTables
      use ImproperAngleFunctions
      use InterEnergy_LJ_Electro
      use IntraEnergy_LJ_Electro
      use TorsionalFunctions

      use E_Interface_LJ_Q_Diststore, only: Shift_EnergyCalc_LJ_Q_DStore
      use PairStorage, only: useDistStore
      implicit none

      logical, intent(in), optional :: useInter
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

      if (useDistStore) then
         call Shift_EnergyCalc_LJ_Q_DStore(E_Inter, E_Intra, disp, PairList, dETable, useIntra, rejMove, useInter)
         return
      end if

      nDisp = size(disp)
      rejMove = .false.
      E_Inter = 0E0
      E_Intra = 0E0
      E_NonBond = 0E0
      E_Stretch = 0E0
      E_Bend = 0E0
      E_Torsion = 0E0
      E_Improper = 0E0
      dETable = 0E0

      E_Inter_Diff = 0E0
      E_NBond_Diff = 0E0
      E_Strch_Diff = 0E0
      E_Bend_Diff = 0E0
      E_Tors_Diff = 0E0

      if (present(useInter)) then
         interSwitch = useInter
      else
         interSwitch = .true.
      end if

!     Begin by calculating the intermolecular potential. If any atoms overlap the move will be rejected
!     immediately. If the calling function has specfied that useInter is .false. then this part will be
!     skipped and only the intra molecular potential will be calculated.
      if (interSwitch) then
         if (NTotal .gt. 1) then
            dETable = 0E0
            PairList = 0E0
            call Shift_ECalc_Inter(E_Inter, disp, PairList, dETable, rejMove)

            if (rejMove) then
               return
            end if
            E_Inter_Diff = E_Inter
            if (E_Inter*beta .gt. softCutOff) then
               rejMove = .true.
               return
            end if
         end if

!        Using the data collected from the intermolecular function, check to see that the new position
!        satisfies the cluster criteria.
         nIndx = MolArray(disp(1)%molType)%mol(disp(1)%molIndx)%indx
         if (distCriteria) then
            if (any(disp(1:nDisp)%atmIndx .eq. 1)) then
               call Shift_DistanceCriteria(PairList, nIndx, rejMove)
            end if
         else
            call Shift_EnergyCriteria(PairList, nIndx, rejMove)
         end if
         if (rejMove) then
            return
         end if
      end if

!      This block contains the calculations for all Intramolecular interactions.  For moves
!      that do not change the internal configuration of a molecule (Translation for example)
!      the variable "useIntra" is set to false which skips these calculations.
      if (regrowType(disp(1)%molType) .ne. 0) then
         if (useIntra(1)) then
            call Shift_ECalc_IntraNonBonded(E_NonBond, disp)
            E_NBond_Diff = E_NonBond
         end if
         if (useIntra(2)) then
            call Shift_ECalc_BondStretch(E_Stretch, disp)
            E_Strch_Diff = E_Stretch
         end if
         if (useIntra(3)) then
            call Shift_ECalc_Bending(E_Bend, disp)
            E_Bend_Diff = E_Bend
         end if
         if (useIntra(4)) then
            call Shift_ECalc_Torsional(E_Torsion, disp)
            E_Tors_Diff = E_Torsion
         end if
!           call Shift_ECalc_Improper(E_Improper, disp)
         E_Intra = E_NonBond + E_Stretch + E_Bend + E_Torsion + E_Improper
      end if

   end subroutine
!=============================================================================
   subroutine SwapIn_EnergyCalc_LJ_Q(E_Inter, E_Intra, PairList, dETable, rejMove, useInter)
      use InterEnergy_LJ_Electro
      use IntraEnergy_LJ_Electro
      use BondStretchFunctions
      use BendingFunctions
      use TorsionalFunctions
      use ImproperAngleFunctions
      use EnergyCriteria
      use EnergyTables
      use Coords
      use CBMC_Variables

      use E_Interface_LJ_Q_Diststore, only: SwapIn_EnergyCalc_LJ_Q_DStore
      use PairStorage, only: useDistStore
      implicit none

      logical, intent(out) :: rejMove
      logical, intent(in), optional :: useInter
      real(dp), intent(out) :: E_Inter, E_Intra
      real(dp), intent(inout) :: PairList(:), dETable(:)

      logical :: interSwitch
      real(dp) :: E_NonBond, E_Stretch, E_Bend
      real(dp) :: E_Torsion, E_Improper

      if (useDistStore) then
         call SwapIn_EnergyCalc_LJ_Q_DStore(E_Inter, E_Intra, PairList, dETable, rejMove, useInter)
         return
      end if

      rejMove = .false.

      if (present(useInter)) then
         interSwitch = useInter
      else
         interSwitch = .true.
      end if

      E_Inter = 0E0
      E_Intra = 0E0
      E_NonBond = 0E0
      E_Stretch = 0E0
      E_Bend = 0E0
      E_Torsion = 0E0
      E_Improper = 0E0
      PairList = 0E0
      dETable = 0E0

      E_Inter_Diff = 0E0
      E_NBond_Diff = 0E0
      E_Strch_Diff = 0E0
      E_Bend_Diff = 0E0
      E_Tors_Diff = 0E0

      if (interSwitch) then
         call NewMol_ECalc_Inter(E_Inter, PairList, dETable, rejMove)
         if (rejMove) then
            return
         end if
      end if

!     This block contains the calculations for all Intramolecular interactions.
      if (regrowType(newMol%molType) .ne. 0) then
         call NewMol_ECalc_IntraNonBonded(E_NonBond)
         call NewMol_ECalc_BondStretch(E_Stretch)
         call NewMol_ECalc_Bending(E_Bend)
         call NewMol_ECalc_Torsional(E_Torsion)
!        call NewMol_ECalc_Improper(E_Improper)
         E_Intra = E_NonBond + E_Stretch + E_Bend + E_Torsion + E_Improper
         E_NBond_Diff = E_NonBond
         E_Strch_Diff = E_Stretch
         E_Bend_Diff = E_Bend
         E_Tors_Diff = E_Torsion
      end if

      if (interSwitch) then
         E_Inter_Diff = E_Inter
      end if

   end subroutine
!=============================================================================
!     This function contains the energy calculations that are used when a molecule
!     has been selected for removal.
   subroutine SwapOut_EnergyCalc_LJ_Q(E_Inter, E_Intra, nType, nMol, dETable, useInter)
      use InterEnergy_LJ_Electro
      use IntraEnergy_LJ_Electro
      use BondStretchFunctions
      use BendingFunctions
      use TorsionalFunctions
      use ImproperAngleFunctions
      use EnergyCriteria
      use EnergyTables
      use Coords
      use CBMC_Variables

      use E_Interface_LJ_Q_Diststore, only: SwapOut_EnergyCalc_LJ_Q_DStore
      use PairStorage, only: useDistStore
      implicit none

      logical, intent(in), optional :: useInter
      real(dp), intent(out) :: E_Inter, E_Intra
      integer, intent(in) :: nType, nMol
      real(dp), intent(inout) :: dETable(:)

      logical :: interSwitch
      real(dp) :: E_NonBond, E_Stretch, E_Bend
      real(dp) :: E_Torsion, E_Improper

      if (useDistStore) then
         call SwapOut_EnergyCalc_LJ_Q_DStore(E_Inter, E_Intra, nType, nMol, dETable, useInter)
         return
      end if

      E_Inter = 0E0
      E_Intra = 0E0
      E_NonBond = 0E0
      E_Stretch = 0E0
      E_Bend = 0E0
      E_Torsion = 0E0
      E_Improper = 0E0

      E_Inter_Diff = 0E0
      E_NBond_Diff = 0E0
      E_Strch_Diff = 0E0
      E_Bend_Diff = 0E0
      E_Tors_Diff = 0E0
      dETable = 0E0

      if (present(useInter)) then
         interSwitch = useInter
      else
         interSwitch = .true.
      end if

      if (interSwitch) then
         call Mol_ECalc_Inter(nType, nMol, dETable, E_Inter)
         E_Inter = -E_Inter
      end if

!     This block contains the calculations for all Intramolecular interactions.
      if (regrowType(nType) .ne. 0) then
         call Mol_ECalc_IntraNonBonded(nType, nMol, E_NonBond)
         call Mol_ECalc_BondStretch(nType, nMol, E_Stretch)
         call Mol_ECalc_Bending(nType, nMol, E_Bend)
         call Mol_ECalc_Torsional(nType, nMol, E_Torsion)
!        call Mol_ECalc_Improper(nType, nMol, E_Improper)
         E_Intra = E_NonBond + E_Stretch + E_Bend + E_Torsion + E_Improper
         E_Intra = -E_Intra
         E_NBond_Diff = -E_NonBond
         E_Strch_Diff = -E_Stretch
         E_Bend_Diff = -E_Bend
         E_Tors_Diff = -E_Torsion
      end if

      if (interSwitch) then
         E_Inter_Diff = E_Inter
      end if
   end subroutine

!=============================================================================
!     This function contains the energy calculations that are used when a molecule
!     has been selected for removal.
   subroutine Update_SubEnergies
      use EnergyTables
      implicit none

      E_Inter_T = E_Inter_T + E_Inter_Diff
      E_NBond_T = E_NBond_T + E_NBond_Diff
      E_Stretch_T = E_Stretch_T + E_Strch_Diff
      E_Bend_T = E_Bend_T + E_Bend_Diff
      E_Torsion_T = E_Torsion_T + E_Tors_Diff

      E_Inter_Diff = 0E0
      E_NBond_Diff = 0E0
      E_Strch_Diff = 0E0
      E_Bend_Diff = 0E0
      E_Tors_Diff = 0E0
      ! print *, "Default Updated Energies:", E_Inter_T
   end subroutine
!=============================================================================

end module

