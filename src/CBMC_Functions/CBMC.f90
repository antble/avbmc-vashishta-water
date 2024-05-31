module CBMC_Module
contains
!=============================================================
   subroutine CBMC(E_T, acc_x, atmp_x)
      use SimParameters
      use CBMC_Variables
      use CBMC_Utility
      use Coords
!      use E_Interface
      use EnergyPointers, only: Shift_ECalc, Update_SubEnergies, Update_SubEnergies_Vashishta
      use EnergyCriteria
      use EnergyTables
      use Forcefield
      use DistanceCriteria
      use IndexingFunctions
      use PairStorage, only: UpdateDistArray, useDistStore
      use Pressure_LJ_Electro, only: Shift_PressCalc_Inter
      use UmbrellaSamplingNew, only: useUmbrella, GetUmbrellaBias_Disp
      implicit none

      real(dp), intent(inout) :: E_T, acc_x, atmp_x

      logical, parameter :: useIntra(1:4) = [.true., .true., .true., .true.]

      logical :: rejMove
      logical :: regrown(1:maxAtoms)
      logical :: regrowDirection

      integer :: i, j, iPath, nType, nMol, nIndx, nMove, nPath, nAtom
      integer :: iAtom, nAtomLoc, nDisp
      integer :: nGrow, GrowFrom(1:maxAtoms), GrowPrev(1:maxAtoms), GrowNum(1:maxAtoms), TorNum(1:maxAtoms)
      integer :: TorList(1:maxAtoms, 1:maxBranches), GrowList(1:maxAtoms, 1:maxBranches)

      real(dp) :: grnd, ranNum, sumInt
      real(dp) :: dx, dy, dz
      real(dp) :: E_Diff, biasDiff
      real(dp) :: PairList(1:maxMol)
      real(dp) :: dETable(1:maxMol)
      real(dp) :: rosenProb_New, rosenProb_Old, rosenRatio
      real(dp) :: E_Inter, E_Intra

      type(displacement) :: disp(1:maxAtoms)

      rejMove = .false.

!     Randomly Select a Molecule from the cluster and obtain its molecule type and index
      ranNum = grnd()
      sumInt = probTypeCBMC(1)
      nType = 1
      do while (sumInt .lt. ranNum)
         nType = nType + 1
         sumInt = probTypeCBMC(nType)
      end do
      if (NPART(nType) .eq. 0d0) then
         return
      end if
      atmp_x = atmp_x + 1d0
      nMol = floor(NPART(nType)*grnd() + 1d0)
      nIndx = molArray(nType)%mol(nMol)%indx

      regrown = .true.
      select case (regrowType(nType))
      case (1)
         regrown = .false.
         regrown(1) = .true.
         call Simple_Partial_ConfigGen(nType, nIndx, nMol, rosenProb_New, rosenProb_Old, rejMove)
         if (rejMove) then
            return
         end if
      case (2)
         !Begin by selecting an atom in the chain to serve as the regrowth point
         nAtom = floor(grnd()*pathArray(nType)%pathMax(1) + 1d0)
         !Choose the direction to regrow the chain
         if (grnd() .gt. 0.5d0) then
            regrowDirection = .true.
         else
            regrowDirection = .false.
         end if
         !If regrowDirection = .true. then regrowth is performed in the positive chain direction, otherwise
         !it is done in the negative direction.  In this step the
         do i = 1, pathArray(nType)%pathMax(1)
            if (pathArray(nType)%path(1, i) .eq. nAtom) then
               nAtomLoc = i
            end if
         end do

         if (regrowDirection) then
            do iPath = nAtomLoc, pathArray(nType)%pathMax(1)
               iAtom = pathArray(nType)%path(1, iPath)
               regrown(iAtom) = .false.
            end do
         else
            do iPath = 1, nAtomLoc
               iAtom = pathArray(nType)%path(1, iPath)
               regrown(iAtom) = .false.
            end do
         end if
         !In the event that all atoms have been selected for regrowth, randomly pick one chain's terminal
         !atoms to remain grown.
         if (all(regrown .eqv. .false.)) then
            regrown(nAtom) = .true.
         end if

         call StraightChain_Partial_ConfigGen(nType, nMol, regrown(1:maxatoms), regrowDirection, nAtom, rosenProb_New, rejMove)
         if (rejMove) then
            return
         end if

         call StraightChain_Partial_ConfigGen_Reverse(nType, nMol, regrown(1:maxatoms), regrowDirection, nAtom, rosenProb_Old)

      case (3)
         nPath = floor(grnd()*pathArray(nType)%nPaths + 1d0)
         nAtomLoc = floor(grnd()*pathArray(nType)%pathMax(nPath) + 1d0)
         nAtom = pathArray(nType)%path(nPath, nAtomLoc)

         call Schedule_BranchedMol_Growth(nType, nPath, nAtomLoc, nAtom, nGrow, regrown, GrowFrom, GrowPrev, &
                                          GrowNum, GrowList, TorNum, TorList)

         call BranchedMol_Partial_ConfigGen(nType, nMol, regrown(1:maxatoms), nGrow, &
                                            GrowFrom, GrowPrev, GrowNum, GrowList, TorNum, TorList, rosenProb_New, rejMove)
         if (rejMove) then
            return
         end if

         call BranchedMol_Partial_ConfigGen_Reverse(nType, nMol, regrown(1:maxatoms), nGrow, &
                                                    GrowFrom, GrowPrev, GrowNum, GrowList, TorNum, TorList, rosenProb_Old)

         if (all(regrown .eqv. .false.)) then
            regrown(nAtom) = .true.
         end if

      case default
         write (*, *) "ERROR! CBMC is not implimented for this regrowth type! :", regrowType(nType)
         stop
      end select
!"
      nDisp = 0
      do iAtom = 1, nAtoms(nType)
         if (.not. regrown(iAtom)) then
            nDisp = nDisp + 1
            disp(nDisp)%molType = int(nType, atomIntType)
            disp(nDisp)%molIndx = int(nMol, atomIntType)
            disp(nDisp)%atmIndx = int(iAtom, atomIntType)

            disp(nDisp)%x_old => MolArray(nType)%mol(nMol)%x(iAtom)
            disp(nDisp)%y_old => MolArray(nType)%mol(nMol)%y(iAtom)
            disp(nDisp)%z_old => MolArray(nType)%mol(nMol)%z(iAtom)

            disp(nDisp)%x_new = newMol%x(iAtom)
            disp(nDisp)%y_new = newMol%y(iAtom)
            disp(nDisp)%z_new = newMol%z(iAtom)

!        else
!          nDisp = nDisp + 1
!          disp(iAtom)%molType = int(nType,2)
!          disp(iAtom)%molIndx = int(nMol,2)
!          disp(iAtom)%atmIndx = int(iAtom,2)

!          disp(iAtom)%x_old => MolArray(nType)%mol(nMol)%x(iAtom)
!          disp(iAtom)%y_old => MolArray(nType)%mol(nMol)%y(iAtom)
!          disp(iAtom)%z_old => MolArray(nType)%mol(nMol)%z(iAtom)

!          disp(iAtom)%x_new = disp(iAtom)%x_old
!          disp(iAtom)%y_new = disp(iAtom)%y_old
!          disp(iAtom)%z_new = disp(iAtom)%z_old
         end if
      end do

      E_Inter = 0d0
      E_Intra = 0d0
!      call Shift_EnergyCalc(E_Inter, E_Intra, disp(1:nDisp), PairList, dETable, useIntra, rejMove)
      call Shift_ECalc(E_Inter, E_Intra, disp(1:nDisp), PairList, dETable, useIntra, rejMove)
      if (rejMove) then
         return
      end if

      rosenRatio = rosenProb_Old/rosenProb_New

      biasDiff = 0E0
      if (useUmbrella) then
         call GetUmbrellaBias_Disp(disp(1:nDisp), biasDiff, rejMove)
         if (rejMove) then
            return
         end if
      end if

!      if(rosenRatio*exp(-beta*E_Inter + biasDiff) .gt. grnd()) then
      if (log(rosenRatio) - beta*E_Inter + biasDiff .gt. log(grnd())) then
         if (calcPressure) then
            call Shift_PressCalc_Inter(P_Diff, disp(1:nDisp))
            pressure = pressure + P_Diff
!          write(*,*) pressure, P_Diff
         end if
         do iAtom = 1, nDisp
            disp(iAtom)%x_old = disp(iAtom)%x_new
            disp(iAtom)%y_old = disp(iAtom)%y_new
            disp(iAtom)%z_old = disp(iAtom)%z_new
         end do
         E_T = E_T + E_Inter + E_Intra
         ETable = ETable + dETable
         acc_x = acc_x + 1d0
         if (useDistStore) then
            call UpdateDistArray
         end if
         if (distCriteria) then
!          call NeighborUpdate_Distance(PairList, nIndx)
            call NeighborUpdate_Distance(PairList, nIndx)
         else
            call NeighborUpdate(PairList, nIndx)
         end if

         !call Update_SubEnergies
         select case (trim(adjustl(ForceFieldName)))
         case('Vashishta')
            call Update_SubEnergies_Vashishta
            print *, "CBMC: Vashishta sub energies update ", (trim(adjustl(ForceFieldName)))
         case default
            call Update_SubEnergies
            print *, "CBMC: Default sub energies update ", (trim(adjustl(ForceFieldName)))
         end select
      end if

   end subroutine
!============================================================
end module

