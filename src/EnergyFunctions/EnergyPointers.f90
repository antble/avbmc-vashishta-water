!================================================================================================================
!    This file contains the Energy Function Pointer Module.  This module exists to allows
!    the end user to freely switch between different types of forcefield potentials via the inputfile.
!    This is done by
!================================================================================================================
module EnergyPointers
   use E_Interface_LJ_Q
   use E_Interface_Pedone
   use E_Interface_Tersoff
   use E_Interface_Vashishta
   use Rosenbluth_Functions_LJ_Q
   use Rosenbluth_Functions_Vashishta_Q
   use Rosenbluth_Functions_Pedone
   use Rosenbluth_Functions_Tersoff
   use InterEnergy_LJ_Electro, only: QuickNei_ECalc_Inter_LJ_Q
   use InterEnergy_Pedone, only: QuickNei_ECalc_Inter_Pedone
   use InterEnergy_Tersoff, only: QuickNei_ECalc_Inter_Tersoff
   ! use InterEnergy_Vashishta_Electro, only: QuickNei_ECalc_Inter_Vashishta

   interface
      subroutine DetailedInterface(E_T, rejMove)
         use VarPrecision
         implicit none
         logical, intent(inout) :: rejMove
         real(dp), intent(inout) :: E_T
      end subroutine
   end interface

   interface
      subroutine ShiftInterface(E_Inter, E_Intra, disp, PairList, dETable, useIntra, rejMove, useInter)
         use CoordinateTypes
         use VarPrecision
         implicit none
         logical, intent(in), optional :: useInter
         logical, intent(in) :: useIntra(1:4)
         type(Displacement), intent(in) :: disp(:)
         logical, intent(inout) :: rejMove
         real(dp), intent(inout) :: dETable(:)
         real(dp), intent(out) :: E_Intra, E_Inter
         real(dp), intent(InOut) :: PairList(:)
      end subroutine
   end interface

   interface
      subroutine SwapInInterface(E_Inter, E_Intra, PairList, dETable, rejMove, useInter)
         use VarPrecision
         implicit none
         logical, intent(out) :: rejMove
         logical, intent(in), optional :: useInter
         real(dp), intent(out) :: E_Inter, E_Intra
         real(dp), intent(inout) :: PairList(:), dETable(:)
      end subroutine
   end interface

   interface
      subroutine SwapOutInterface(E_Inter, E_Intra, nType, nMol, dETable, useInter)
         use VarPrecision
         implicit none
         logical, intent(in), optional :: useInter
         real(dp), intent(out) :: E_Inter, E_Intra
         integer, intent(in) :: nType, nMol
         real(dp), intent(inout) :: dETable(:)
      end subroutine
   end interface

   interface
      subroutine RosenMolNewInterface(nRosen, nType, included, E_Trial, overlap)
         use VarPrecision
         implicit none
         logical, intent(in) :: included(:)
         integer, intent(in) :: nType, nRosen

         logical, intent(out) :: overlap
         real(dp), intent(out) :: E_Trial
      end subroutine
   end interface

   ! an interface for Vashishta FF
   interface
      subroutine RosenMolOldInterface_V(mol_x, mol_y, mol_z, nType, included, E_Trial, nMol, X)
         use VarPrecision
         implicit none
         logical, intent(in) :: included(:)
         integer, intent(in) :: nType, nMol, X
         real(dp), intent(in) :: mol_x(:), mol_y(:), mol_z(:)
         real(dp), intent(out) :: E_Trial
      end subroutine
   end interface

   interface
      subroutine RosenMolOldInterface(mol_x, mol_y, mol_z, nType, included, E_Trial)
         use VarPrecision
         implicit none
         logical, intent(in) :: included(:)
         integer, intent(in) :: nType
         real(dp), intent(in) :: mol_x(:), mol_y(:), mol_z(:)
         real(dp), intent(out) :: E_Trial
      end subroutine
   end interface

   interface
      subroutine QuickInterface(jType, jMol, rejMove)
         implicit none
         integer, intent(in) :: jType, jMol
         logical, intent(out) :: rejMove
      end subroutine
   end interface

   procedure(DetailedInterface), pointer :: Detailed_ECalc => NULL()
   procedure(ShiftInterface), pointer  :: Shift_ECalc => NULL()
   procedure(SwapInInterface), pointer :: SwapIn_ECalc => NULL()
   procedure(SwapOutInterface), pointer :: SwapOut_ECalc => NULL()

   procedure(RosenMolNewInterface), pointer :: Rosen_Mol_New => NULL()
   procedure(RosenMolOldInterface), pointer :: Rosen_Mol_Old => NULL()
   procedure(RosenMolOldInterface_V), pointer :: Rosen_Mol_Old_V => NULL()

!     procedure(RosenAtomNewInterface), pointer :: Rosen_Atom_New => NULL()
!     procedure(RosenAtomOldInterface), pointer :: Rosen_Atom_Old => NULL()
   procedure(QuickInterface), pointer :: Quick_Nei_ECalc => NULL()

end module
!================================================================================================================


