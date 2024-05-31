module DistanceGen
!========================================================================
   interface
      real(dp) function distPointer

      end function
   end interface

   procedure(distPointer), pointer :: genAVMBCDistance => null()

!------------------------------------------------------------------------
contains

   real(dp) function UniformDist
   end function
!========================================================================
end module
