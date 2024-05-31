!=====================================================================================
module UmbrellaTypes
   use VarPrecision

   type BiasVariablePointer
      integer :: varType
      integer, pointer :: intVar
      real(dp), pointer :: realVar
   end type

   abstract interface
      subroutine UDispFunc(disp)
         use CoordinateTypes
         implicit none
         type(Displacement), intent(in) :: disp(:)
      end subroutine
   end interface

   abstract interface
      subroutine USwapOutFunc(nType, nMol)
         use CoordinateTypes
         implicit none
         integer, intent(in) :: nType, nMol
      end subroutine
   end interface

   type DispUmbrellaArray
      procedure(UDispFunc), pointer, nopass :: func
   end type

   type SwapInUmbrellaArray
      procedure(), pointer, nopass :: func
   end type

   type SwapOutUmbrellaArray
      procedure(USwapOutFunc), pointer, nopass :: func
   end type

end module
!=====================================================================================
