!=============================================================
module RandomTools
contains
!=======================================================
   subroutine SelectBin_IntTable(integralTable, nBin)
      use VarPrecision
      implicit none
      real(dp), intent(in) :: integralTable(:)
      integer, intent(out) :: nBin
      integer :: i, iBin, binMin, binMax, binMid
      integer :: lowIndx, highIndx
      real(dp) ::  grnd, ranNum

!      binMin = LBOUND(integralTable,1)
      binMax = size(integralTable) + binMin - 1
      ranNum = grnd()

      do i = 1, 7
         binMid = nint(ranNum*binMax + (1E0_dp - ranNum)*binMin)
         if (integralTable(binMid) .gt. ranNum) then
            binMax = binMid
         elseif (integralTable(binMid) .lt. ranNum) then
            binMin = binMid
         end if
      end do

      nBin = binMin
      do while (nBin .lt. binMax)
         if (integralTable(nBin) .gt. ranNum) then
            exit
         end if
         nBin = nBin + 1
      end do
!      nBin = nBin - 1

!      write(3,*) ranNum, nBin, integralTable(nBin), integralTable(nBin-1), integralTable(nBin+1)

   end subroutine
!=============================================================
end module
!=============================================================
