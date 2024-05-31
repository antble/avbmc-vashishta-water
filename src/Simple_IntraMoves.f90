!===========================================================================================
      subroutine Simple_Intra(E_T, acc_x, atmp_x)
         use SimParameters
         implicit none
         real(dp), intent(inout) :: E_T, acc_x, atmp_x
         real(dp) :: ran_num, grnd

         atmp_x = atmp_x + 1d0
         ran_num = grnd()
         if (ran_num .lt. 1d0/2d0) then
            call Simple_Bond_Stretch(E_T, acc_x)
         else
            call Simple_Angle_Bend(E_T, acc_x)
         end if

      end subroutine

!===========================================================================================
      subroutine Simple_Bond_Stretch(E_T, acc_x)
         use SimParameters
         use IndexingFunctions
         use Coords
         use Forcefield
         use E_Interface
         use EnergyTables
         use EnergyCriteria
         use DistanceCriteria
         use CBMC_Variables
         implicit none

         real(dp), intent(inout) :: E_T, acc_x

         logical, parameter :: useIntra(1:4) = [.false., .true., .false., .false.]

         logical rejMove
         integer :: nType, nMol, nIndx, nMove, nBondMove
         integer :: nBondType
         integer :: mem1, mem2
         real(dp) :: grnd
         real(dp) :: dx, dy, dz
         real(dp) :: r_new, r_old
         real(dp) :: E_Inter, E_Intra
         type(displacement) :: disp(1:1)
         real(dp) :: PairList(1:maxMol)
         real(dp) :: dETable(1:maxMol)
         real(dp) :: k_bond, r_eq, Prob

         rejMove = .false.
!     Randomly Select a Particle from the cluster and obtain its molecule type and index
         nMove = floor(NTotal*grnd() + 1d0)
         call Get_MolIndex(nMove, NPart, nType, nMol)
         if (regrowType(nType) .eq. 0) return
         do while (regrowType(nType) .ne. 1)
            nMove = floor(NTotal*grnd() + 1d0)
            call Get_MolIndex(nMove, NPart, nType, nMol)
         end do

         nIndx = MolArray(nType)%mol(nMol)%indx

!      Choose a Random Bond
         nBondMove = floor(nBonds(nType)*grnd() + 1d0)
         nBondType = bondArray(nType, nBondMove)%bondType

         if (topolArray(nType)%atom(bondArray(nType, nBondMove)%bondMembr(1)) .eq. 1) then
            mem1 = bondArray(nType, nBondMove)%bondMembr(2)
            mem2 = bondArray(nType, nBondMove)%bondMembr(1)
         else
            mem1 = bondArray(nType, nBondMove)%bondMembr(1)
            mem2 = bondArray(nType, nBondMove)%bondMembr(2)
         end if
         k_bond = bondData(nBondType)%k_eq
         r_eq = bondData(nBondType)%r_eq
!      Generate a random bond displacement
         call GenerateBondLength(r_new, k_bond, r_eq, Prob)

         dx = MolArray(nType)%mol(nMol)%x(mem2) - MolArray(nType)%mol(nMol)%x(mem1)
         dy = MolArray(nType)%mol(nMol)%y(mem2) - MolArray(nType)%mol(nMol)%y(mem1)
         dz = MolArray(nType)%mol(nMol)%z(mem2) - MolArray(nType)%mol(nMol)%z(mem1)
         r_old = dsqrt(dx*dx + dy*dy + dz*dz)

!     Construct the Displacement Vectors for each atom in the molecule that was chosen.'
         disp(1)%molType = int(nType, 2)
         disp(1)%molIndx = int(nMol, 2)
         disp(1)%atmIndx = int(mem2, 2)

         disp(1)%x_old => MolArray(nType)%mol(nMol)%x(mem2)
         disp(1)%y_old => MolArray(nType)%mol(nMol)%y(mem2)
         disp(1)%z_old => MolArray(nType)%mol(nMol)%z(mem2)

         disp(1)%x_new = dx*(r_new/r_old) + MolArray(nType)%mol(nMol)%x(mem1)
         disp(1)%y_new = dy*(r_new/r_old) + MolArray(nType)%mol(nMol)%y(mem1)
         disp(1)%z_new = dz*(r_new/r_old) + MolArray(nType)%mol(nMol)%z(mem1)

!     Calculate the Energy Associated with this move and ensure it conforms to
!     the cluster criteria.
         E_Inter = 0d0
         E_Intra = 0d0
         call Shift_EnergyCalc(E_Inter, E_Intra, disp, PairList, dETable, useIntra, rejMove)
         if (rejMove) return
!     Calculate Acceptance and determine if the move is accepted or not
         if (E_Inter .le. 0d0) then
            disp(1)%x_old = disp(1)%x_new
            disp(1)%y_old = disp(1)%y_new
            disp(1)%z_old = disp(1)%z_new
            E_T = E_T + E_Inter + E_Intra
            ETable = ETable + dETable
            acc_x = acc_x + 1d0
            if (distCriteria) then
               call NeighborUpdate_Distance(PairList, nIndx)
            else
               call NeighborUpdate(PairList, nIndx)
            end if
         elseif (exp(-beta*E_Inter) .gt. grnd()) then
            disp(1)%x_old = disp(1)%x_new
            disp(1)%y_old = disp(1)%y_new
            disp(1)%z_old = disp(1)%z_new
            E_T = E_T + E_Inter + E_Intra
            ETable = ETable + dETable
            acc_x = acc_x + 1d0
            if (distCriteria) then
               call NeighborUpdate_Distance(PairList, nIndx)
            else
               call NeighborUpdate(PairList, nIndx)
            end if
         end if

      end subroutine
!===========================================================================================
      subroutine Simple_Angle_Bend(E_T, acc_x)
         use SimParameters, only: beta, maxMol, NTotal, NPart, distCriteria
         use IndexingFunctions, only: Get_MolIndex
         use Coords, only: MolArray
         use Forcefield
         use E_Interface
         use EnergyTables
         use EnergyCriteria
         use DistanceCriteria
         use CBMC_Variables
         implicit none

         real(dp), intent(inout) :: E_T, acc_x

         logical, parameter :: useIntra(1:4) = [.false., .false., .true., .false.]

         logical rejMove
         integer :: nType, nMol, nIndx, nMove, nBend
         integer :: nBendType
         integer :: mem1, mem2, mem3
         real(dp) :: grnd
         real(dp) :: r1, r2
         real(dp) :: angle_new
         real(dp) :: E_Inter, E_Intra
         type(displacement) :: disp(1:1)
         real(dp) :: PairList(1:maxMol)
         real(dp) :: dETable(1:maxMol)
         real(dp) :: k_bend, ang_eq, Prob
         type(SimpleAtomCoords) :: v1, v2, v3

         rejMove = .false.
!     Randomly Select a Particle from the cluster and obtain its molecule type and index
         nMove = floor(NTotal*grnd() + 1d0)
         call Get_MolIndex(nMove, NPart, nType, nMol)
         do while (regrowType(nType) .ne. 1)
            nMove = floor(NTotal*grnd() + 1d0)
            call Get_MolIndex(nMove, NPart, nType, nMol)
         end do

         nIndx = MolArray(nType)%mol(nMol)%indx

!      Choose a Random Bond
         nBend = floor(nAngles(nType)*grnd() + 1d0)
         nBendType = bendArray(nType, nBend)%bendType

         k_bend = bendData(nBendType)%k_eq
         ang_eq = bendData(nBendType)%ang_eq
!      Generate a random bond displacement

         mem1 = bendArray(nType, nBend)%bendMembr(1)
         mem2 = bendArray(nType, nBend)%bendMembr(2)
         mem3 = bendArray(nType, nBend)%bendMembr(3)

         v1%x = MolArray(nType)%mol(nMol)%x(mem1) - MolArray(nType)%mol(nMol)%x(mem2)
         v1%y = MolArray(nType)%mol(nMol)%y(mem1) - MolArray(nType)%mol(nMol)%y(mem2)
         v1%z = MolArray(nType)%mol(nMol)%z(mem1) - MolArray(nType)%mol(nMol)%z(mem2)
         r1 = dsqrt(v1%x*v1%x + v1%y*v1%y + v1%z*v1%z)

         v2%x = MolArray(nType)%mol(nMol)%x(mem3) - MolArray(nType)%mol(nMol)%x(mem2)
         v2%y = MolArray(nType)%mol(nMol)%y(mem3) - MolArray(nType)%mol(nMol)%y(mem2)
         v2%z = MolArray(nType)%mol(nMol)%z(mem3) - MolArray(nType)%mol(nMol)%z(mem2)
         r2 = dsqrt(v2%x*v2%x + v2%y*v2%y + v2%z*v2%z)

!     Randomly generate a bond angle and
         call GenerateBendAngle(angle_new, k_bend, ang_eq, Prob)
         if (grnd() .lt. 0.5d0) then
            call Generate_UnitCone(v2, r1, angle_new, v3)
            disp(1)%molType = nType
            disp(1)%molIndx = nMol
            disp(1)%atmIndx = mem1

            disp(1)%x_old => MolArray(nType)%mol(nMol)%x(mem1)
            disp(1)%y_old => MolArray(nType)%mol(nMol)%y(mem1)
            disp(1)%z_old => MolArray(nType)%mol(nMol)%z(mem1)

            disp(1)%x_new = v3%x + MolArray(nType)%mol(nMol)%x(mem2)
            disp(1)%y_new = v3%y + MolArray(nType)%mol(nMol)%y(mem2)
            disp(1)%z_new = v3%z + MolArray(nType)%mol(nMol)%z(mem2)
         else
            call Generate_UnitCone(v1, r2, angle_new, v3)
            disp(1)%molType = nType
            disp(1)%molIndx = nMol
            disp(1)%atmIndx = mem3

            disp(1)%x_old => MolArray(nType)%mol(nMol)%x(mem3)
            disp(1)%y_old => MolArray(nType)%mol(nMol)%y(mem3)
            disp(1)%z_old => MolArray(nType)%mol(nMol)%z(mem3)

            disp(1)%x_new = v3%x + MolArray(nType)%mol(nMol)%x(mem2)
            disp(1)%y_new = v3%y + MolArray(nType)%mol(nMol)%y(mem2)
            disp(1)%z_new = v3%z + MolArray(nType)%mol(nMol)%z(mem2)
         end if

!     Calculate the Energy Associated with this move and ensure it conforms to
!     the cluster criteria.
         E_Inter = 0d0
         E_Intra = 0d0
         call Shift_EnergyCalc(E_Inter, E_Intra, disp, PairList, dETable, useIntra, rejMove)
         if (rejMove) return

!     Calculate Acceptance and determine if the move is accepted or not
         if (E_Inter .le. 0d0) then
            disp(1)%x_old = disp(1)%x_new
            disp(1)%y_old = disp(1)%y_new
            disp(1)%z_old = disp(1)%z_new
            E_T = E_T + E_Inter + E_Intra
            ETable = ETable + dETable
            acc_x = acc_x + 1d0
            if (distCriteria) then
               call NeighborUpdate_Distance(PairList, nIndx)
            else
               call NeighborUpdate(PairList, nIndx)
            end if
         elseif (exp(-beta*E_Inter) .gt. grnd()) then
            disp(1)%x_old = disp(1)%x_new
            disp(1)%y_old = disp(1)%y_new
            disp(1)%z_old = disp(1)%z_new
            E_T = E_T + E_Inter + E_Intra
            ETable = ETable + dETable
            acc_x = acc_x + 1d0
            if (distCriteria) then
               call NeighborUpdate_Distance(PairList, nIndx)
            else
               call NeighborUpdate(PairList, nIndx)
            end if
         end if

      end subroutine
!==============================================================================
