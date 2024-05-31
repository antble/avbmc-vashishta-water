!======================================================
!     This module
module VarPrecision
   ! use iso_fortran_env, only:  real64
   ! use, intrinsic :: iso_fortran_env, dp=>real64
   implicit none
   integer, parameter :: dp =  kind(0.0d0)
   ! integer, parameter :: qp = selected_real_kind(33, 4931)
   ! integer, parameter :: dp = real64
!        integer, parameter :: dp = 8
   integer, parameter :: atomIntType = kind(0)
end module
!======================================================
