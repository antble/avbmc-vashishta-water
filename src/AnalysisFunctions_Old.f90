!======================================================================================
      subroutine Calc_BendAngle(iType, iBend)
         use ParallelVar
         use SimParameters
         use ForceField
         use Coords
         use ForceFieldFunctions
         use Constants
         use Histogram
         implicit none
         integer, intent(in) :: iType, iBend
         integer :: iMol
         integer :: bin
         integer :: memb1, memb2, memb3
         real(dp) :: rx12, ry12, rz12, r12
         real(dp) :: rx23, ry23, rz23, r23
         real(dp) :: Angle

         do iMol = 1, NPART(iType)
            memb1 = bendArray(iType, iBend)%bendMembr(1)
            memb2 = bendArray(iType, iBend)%bendMembr(2)
            memb3 = bendArray(iType, iBend)%bendMembr(3)

            rx12 = MolArray(iType)%mol(iMol)%x(memb1) - MolArray(iType)%mol(iMol)%x(memb2)
            ry12 = MolArray(iType)%mol(iMol)%y(memb1) - MolArray(iType)%mol(iMol)%y(memb2)
            rz12 = MolArray(iType)%mol(iMol)%z(memb1) - MolArray(iType)%mol(iMol)%z(memb2)
            r12 = sqrt(rx12**2 + ry12**2 + rz12**2)

            rx23 = MolArray(iType)%mol(iMol)%x(memb3) - MolArray(iType)%mol(iMol)%x(memb2)
            ry23 = MolArray(iType)%mol(iMol)%y(memb3) - MolArray(iType)%mol(iMol)%y(memb2)
            rz23 = MolArray(iType)%mol(iMol)%z(memb3) - MolArray(iType)%mol(iMol)%z(memb2)
            r23 = sqrt(rx23**2 + ry23**2 + rz23**2)

            Angle = rx12*rx23 + ry12*ry23 + rz12*rz23
            Angle = Angle/(r12*r23)
            if (abs(Angle) .gt. 1d0) then
               Angle = sign(real(1, dp), Angle)
            end if
            Angle = acos(Angle)
            bin = floor(Angle*d_ang)
            HistAngle(bin) = HistAngle(bin) + 1d0
         end do

      end subroutine
!======================================================================================
      subroutine Calc_BondStretch(iType, iBond)
         use ParallelVar
         use SimParameters
         use ForceField
         use Coords
         use ForceFieldFunctions
         use Histogram
         implicit none

         integer, intent(in) :: iType, iBond

         integer :: iMol, bin
         integer :: memb1, memb2
         real(dp) :: rx, ry, rz, r

         do iMol = 1, NPART(iType)
            memb1 = bondArray(iType, iBond)%bondMembr(1)
            memb2 = bondArray(iType, iBond)%bondMembr(2)

            rx = MolArray(iType)%mol(iMol)%x(memb1) - MolArray(iType)%mol(iMol)%x(memb2)
            ry = MolArray(iType)%mol(iMol)%y(memb1) - MolArray(iType)%mol(iMol)%y(memb2)
            rz = MolArray(iType)%mol(iMol)%z(memb1) - MolArray(iType)%mol(iMol)%z(memb2)
            r = sqrt(rx**2 + ry**2 + rz**2)
            bin = floor(r*d_r)
            HistDist(bin) = HistDist(bin) + 1d0
         end do

      end subroutine

!======================================================================================
