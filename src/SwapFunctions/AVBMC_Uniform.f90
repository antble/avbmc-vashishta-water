! module AVBMC_Module
!    use VarPrecision
!    real(dp), allocatable :: swapProb(:)
! contains
!===================================================================================
      subroutine AVBMC(E_T, acc_x, atmp_x)
         use SimParameters
         use Constants
         implicit none
         real(dp), intent(inout) :: E_T, atmp_x, acc_x
         real(dp) :: grnd
         atmp_x = atmp_x + 1d0
         if (grnd() .lt. 0.5d0) then
            call AVBMC_Uniform_In(E_T, acc_x)
         else
            call AVBMC_Uniform_Out(E_T, acc_x)
         end if

      end
!===================================================================================
      subroutine AVBMC_Uniform_In(E_T, acc_x)
         use SimParameters
         use Constants
         use ForceField
         ! use E_Interface
         use Coords
         use UmbrellaFunctions
         use ForceField
         use IndexingFunctions
         use EnergyCriteria
         use DistanceCriteria
         use InterEnergy_LJ_Electro
         use EnergyTables
         use AcceptRates
         use CBMC_Variables
         implicit none

         real(dp), intent(inout) :: E_T
         real(dp), intent(inout) :: acc_x
         logical rejMove
         integer :: NDiff(1:nMolTypes)
         integer :: i, nMove, nTargType, nTargMol, nTargIndx, nTarget
         integer :: nType, nIndx, bIndx, nNei
         real(dp) :: grnd
         real(dp) :: dx, dy, dz, r
         real(dp) :: genProbRatio
         real(dp) :: E_Inter, E_Intra, bias_Diff
         real(dp) :: biasOld, biasNew
         real(dp) :: PairList(1:maxMol)
         real(dp) :: dETable(1:maxMol)
         real(dp) :: x1, y1, z1

         if (NTotal .eq. maxMol) return

!      Choose the type of molecule to be inserted
         if (nMolTypes .eq. 1) then
            nType = 1
         else
            nType = floor(nMolTypes*grnd() + 1d0)
         end if
         atmpSwapIn(nType) = atmpSwapIn(nType) + 1d0

         if (NPART(nType) .eq. NMAX(nTYPE)) then
            return
         end if

         NDiff = 0
         NDiff(nType) = 1
         nMove = floor(NTotal*grnd() + 1d0)
         call Get_MolIndex(nMove, NPart, nTargType, nTargMol)
         nTarget = MolArray(nTargType)%mol(nTargMol)%indx

!      Generate the configuration for the newly inserted molecule
         select case (regrowType(nType))
         case (0)
            call Ridgid_ConfigGen(nType)
         case (1)
            call Simple_ConfigGen(nType)
         case default
            write (6, *) "Error! EBias can not regrow a molecule of regrow type:", nType
            stop
         end select

!      Uniformly generate a random distance
         r = Dist_Critr*grnd()**(1d0/3d0)
         if (r .lt. global_r_min) return

         call Generate_UnitSphere(dx, dy, dz)
         dx = r*dx
         dy = r*dy
         dz = r*dz

!      Position the first atom of the new molecule at the randomly generated location.
         x1 = molArray(nTargType)%mol(nTargMol)%x(1) + dx - newMol%x(1)
         y1 = molArray(nTargType)%mol(nTargMol)%y(1) + dy - newMol%y(1)
         z1 = molArray(nTargType)%mol(nTargMol)%z(1) + dz - newMol%z(1)
         do i = 1, nAtoms(nType)
            newMol%x(i) = newMol%x(i) + x1
            newMol%y(i) = newMol%y(i) + y1
            newMol%z(i) = newMol%z(i) + z1
         end do

!      Perform a check to see if the cluster criteria is statisfied or not.
         if (.not. distCriteria) then
            rejMove = .false.
            call QuickNei_ECalc_Inter(nTargType, nTargMol, rejMove)
            if (rejMove) then
               clusterCritRej = clusterCritRej + 1d0
               return
            end if
         end if

!      Calculate the Energy Difference Associated with the move
         E_Inter = 0d0
         E_Intra = 0d0
         call SwapIn_EnergyCalc(E_Inter, E_Intra, PairList, dETable, rejMove)
         if (rejMove) then
            return
         end if

!      if(.not. distCriteria) then
!        rejMove = .false.
!        call SwapIn_EnergyCriteria(nType,PairList,rejMove)
!        if(rejMove) then
!          clusterCritRej = clusterCritRej + 1d0
!          return
!        endif
!      endif

!     Determine the reverse probability of this move.
         call Uniform_GetNeighborCount(nTarget, nNei)

!     Calculate the umbrella sampling bias.
         bIndx = getBiasIndex(NPart, NMAX)
         biasOld = NBias(bIndx)
         bIndx = getNewBiasIndex(NPart, NMAX, NDiff)
         biasNew = NBias(bIndx)
         bias_diff = biasNew - biasOld

!     Calculate acceptance probability and determine if the move is accepted or not
!      genProbRatio = (ProbTarg_Out * ProbSel_Out * avbmc_vol * dble(nMolTypes) * gas_dens(nType)) / ProbTarg_In
         genProbRatio = (dble(NTotal)*avbmc_vol*dble(nMolTypes)*gas_dens(nType))/(dble(NTotal + 1)*dble(nNei + 1))

         if (genProbRatio*exp(-beta*E_Inter + bias_diff) .gt. grnd()) then
            acptSwapIn(nType) = acptSwapIn(nType) + 1d0
            do i = 1, nAtoms(nType)
               molArray(nType)%mol(NPART(nType) + 1)%x(i) = newMol%x(i)
               molArray(nType)%mol(NPART(nType) + 1)%y(i) = newMol%y(i)
               molArray(nType)%mol(NPART(nType) + 1)%z(i) = newMol%z(i)
            end do
            E_T = E_T + E_Inter + E_Intra
            acc_x = acc_x + 1d0
            isActive(molArray(nType)%mol(NPART(nType) + 1)%indx) = .true.
            nIndx = molArray(nType)%mol(NPART(nType) + 1)%indx
            if (distCriteria) then
               call NeighborUpdate_Distance(PairList, nIndx)
            else
               call NeighborUpdate(PairList, nIndx)
            end if
            NTotal = NTotal + 1
            ETable = ETable + dETable
            NPART(nType) = NPART(nType) + 1
!         call Create_NeiETable
            call Update_SubEnergies
         end if
      end subroutine
!===================================================================================
      subroutine AVBMC_Uniform_Out(E_T, acc_x)
         use SimParameters
         use Constants
         ! use E_Interface
         use Coords
         use UmbrellaFunctions
         use ForceField
         use IndexingFunctions
         use EnergyCriteria
         use DistanceCriteria
         use EnergyTables
         use AcceptRates
         use UmbrellaFunctions
         implicit none

         real(dp), intent(inout) :: E_T
         real(dp), intent(inout) :: acc_x

         logical :: rejMove
         integer :: nTarget, nIndx, bIndx, iType, nNei
         integer :: nSel, nType, nMol, nTargMol, nTargType
         integer :: NDiff(1:nMolTypes)
         real(dp) :: grnd
         real(dp) :: genProbRatio
         real(dp) :: bias_diff
         real(dp) :: biasOld, biasNew
         real(dp) :: E_Inter, E_Intra
         real(dp) :: dETable(1:maxMol)
         real(dp) :: rx, ry, rz, dist

         if (NTotal .eq. 1) return

!     Pick a Random Target Particle to Delete
         nTarget = floor(NTotal*grnd() + 1d0)
         call Get_MolIndex(nTarget, NPART, nTargType, nTargMol)
         nTarget = MolArray(nTargType)%mol(nTargMol)%indx

         call Uniform_ChooseNeighbor(nTarget, nSel, nNei)
         nType = typeList(nSel)
         nMol = subIndxList(nSel)

         atmpSwapOut(nType) = atmpSwapOut(nType) + 1d0
         if (NPART(nType) .eq. NMIN(nType)) then
            return
         end if

         NDiff = 0
         NDiff(nType) = -1
         bIndx = getBiasIndex(NPart, NMAX)
         biasOld = NBias(bIndx)
         bIndx = getNewBiasIndex(NPart, NMAX, NDiff)
         biasNew = NBias(bIndx)
         bias_diff = biasNew - biasOld

!     Check to see that the appropriate atoms are within the insertion distance
!     in order to ensure the move is reversible. If not reject the move since
!     the reverse probility is equal to 0.
         if (.not. distCriteria) then
            rx = molArray(nTargType)%mol(nTargMol)%x(1) - molArray(nType)%mol(nMol)%x(1)
            ry = molArray(nTargType)%mol(nTargMol)%y(1) - molArray(nType)%mol(nMol)%y(1)
            rz = molArray(nTargType)%mol(nTargMol)%z(1) - molArray(nType)%mol(nMol)%z(1)
            dist = rx**2 + ry**2 + rz**2
            if (dist .gt. Dist_Critr_sq) then
               return
            end if
         end if

!     Check to see if the deletion of the particle will break the cluster
         rejMove = .false.
         call SwapOut_EnergyCriteria(nSel, rejMove)
         if (rejMove) then
            clusterCritRej = clusterCritRej + 1d0
            return
         end if

!      Calculate the Energy Difference Associated with the move.
         E_Inter = 0d0
         E_Intra = 0d0
         call SwapOut_EnergyCalc(E_Inter, E_Intra, nType, nMol, dETable)
!      genProbRatio = ProbTargIn / (ProbTargOut * ProbSel * dble(nMolTypes) * avbmc_vol * gas_dens(nType))
         genProbRatio = (dble(nNei)*dble(NTotal))/(dble(NTotal - 1)*dble(nMolTypes)*avbmc_vol*gas_dens(nType))

!      Calculate Acceptance and determine if the move is accepted or not
         if (genProbRatio*exp(-beta*E_Inter + bias_diff) .gt. grnd()) then
            acptSwapOut(nType) = acptSwapOut(nType) + 1d0
            molArray(nType)%mol(nMol)%x(1:nAtoms(nType)) = molArray(nType)%mol(NPART(nType))%x(1:nAtoms(nType))
            molArray(nType)%mol(nMol)%y(1:nAtoms(nType)) = molArray(nType)%mol(NPART(nType))%y(1:nAtoms(nType))
            molArray(nType)%mol(nMol)%z(1:nAtoms(nType)) = molArray(nType)%mol(NPART(nType))%z(1:nAtoms(nType))
            E_T = E_T + E_Inter + E_Intra
            nIndx = molArray(nType)%mol(nMol)%indx
            call NeighborUpdate_Delete(nIndx)
            isActive(molArray(nType)%mol(NPART(nType))%indx) = .false.
            acc_x = acc_x + 1d0
            ETable = ETable - dETable
            ETable(nIndx) = ETable(molArray(nType)%mol(NPART(nType))%indx)
            ETable(molArray(nType)%mol(NPART(nType))%indx) = 0d0
            NPART(nType) = NPART(nType) - 1
            NTotal = NTotal - 1
!         call Create_NeiETable
            call Update_SubEnergies
!         call DEBUG_Output_NeighborList
         end if
      end subroutine
!=================================================================================
      pure subroutine Uniform_GetNeighborCount(nTarget, nNei)
         use SimParameters
         use Coords

         implicit none
         integer, intent(in) :: nTarget
         integer, intent(out) :: nNei
         integer :: j

         nNei = 0
         do j = 1, maxMol
            if (j .ne. nTarget) then
               if (NeighborList(nTarget, j)) then
                  nNei = nNei + 1
               end if
            end if
         end do

      end subroutine
!=================================================================================
      subroutine Uniform_ChooseNeighbor(nTarget, nSel, nNei)
         use SimParameters
         use Coords
         implicit none
         integer, intent(in) :: nTarget
         integer, intent(out) :: nSel, nNei
         integer :: i, j
         integer :: ListCur(1:60)
         real(dp) :: grnd

         nNei = 0
         ListCur = 0
         do j = 1, maxMol
            if (NeighborList(nTarget, j)) then
               if (j .ne. nTarget) then
                  nNei = nNei + 1
                  ListCur(nNei) = j
               end if
            end if
         end do

         j = floor(nNei*grnd() + 1d0)
         nSel = ListCur(j)

         if (nSel .eq. 0) then
            write (35, *) "ERROR"
            write (35, *) "NPART:", NPART
            do i = 1, maxMol
               write (35, *) (NeighborList(i, j), j=1, maxMol)
            end do
         end if

      end subroutine
!=================================================================================
! end module