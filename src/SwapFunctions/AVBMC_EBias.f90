!****** Energy Biased Aggregation-Volume-Bias Monte Carlo (AVBMC) algorithm *******
!   This file contains the nessisary functions to impliment the energy biased swap
!   move for cluster simulations.
!===================================================================================
      subroutine AVBMC(E_T, acc_x, atmp_x)
         use SimParameters
         use Constants
         implicit none
         interface
            subroutine AVBMC_EBias_In(E_T, acc_x)
               implicit none
               real(dp), intent(inout) :: E_T, acc_x
            end subroutine
         end interface
         interface
            subroutine AVBMC_EBias_Out(E_T, acc_x)
               implicit none
               real(dp), intent(inout) :: E_T, acc_x
            end subroutine
         end interface

         real(dp), intent(inout) :: E_T, atmp_x, acc_x
         real(dp) :: grnd

         atmp_x = atmp_x + 1d0
         if (grnd() .lt. 0.5d0) then
            call AVBMC_EBias_In(E_T, acc_x)
         else
            call AVBMC_EBias_Out(E_T, acc_x)
         end if

      end subroutine
!===================================================================================
      subroutine AVBMC_EBias_In(E_T, acc_x)
         use SimParameters
         use Constants
         use ForceField
         use E_Interface
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

         interface
            subroutine EBias_Insert_ChooseTarget(nInsType, nTarget, nTargType, nMol, ProbSel)
               implicit none
               integer, intent(in) :: nInsType
               integer, intent(out) :: nTarget, nTargType, nMol
               real(dp), intent(out) :: ProbSel
            end subroutine
         end interface

         interface
            subroutine EBias_Insert_ReverseProbTarget(nTarget, nType, newNeiETable, ProbRev)
               implicit none
               integer, intent(in) :: nTarget, nType
               real(dp), intent(in) :: newNeiETable(:)
               real(dp), intent(out) :: ProbRev
            end subroutine
         end interface

         interface
            subroutine EBias_Insert_ReverseProbSel(nTarget, nType, dE, ProbRev)
               implicit none
               integer, intent(in) :: nTarget, nType
               real(dp), intent(in) :: dE(:)
               real(dp), intent(out) :: ProbRev
            end subroutine
         end interface

         interface
            subroutine Insert_NewNeiETable_Distance(nType, PairList, dE, newNeiTable)
               implicit none
               integer, intent(in) :: nType
               real(dp), intent(in) :: PairList(:)
               real(dp), intent(inout) :: dE(:), newNeiTable(:)
            end subroutine
         end interface

         interface
            subroutine Insert_NewNeiETable(nType, PairList, dE, newNeiTable)
               implicit none
               integer, intent(in) :: nType
               real(dp), intent(in) :: PairList(:)
               real(dp), intent(inout) :: dE(:), newNeiTable(:)
            end subroutine
         end interface

         real(dp), intent(inout) :: E_T
         real(dp), intent(inout) :: acc_x
         logical rejMove
         integer :: NDiff(1:nMolTypes)
         integer :: i, nTargType, nTargMol, nTargIndx, nTarget
         integer :: nType, nIndx, bIndx
         integer :: atmType1, atmType2
         real(dp) :: grnd
         real(dp) :: dx, dy, dz, r
         real(dp) :: genProbRatio
         real(dp) :: E_Inter, E_Intra, bias_Diff
         real(dp) :: biasOld, biasNew
         real(dp) :: PairList(1:maxMol)
         real(dp) :: dETable(1:maxMol)
         real(dp) :: newNeiETable(1:maxMol)
         real(dp) :: x1, y1, z1
         real(dp) :: ProbTarg_In, ProbTarg_Out, ProbSel_Out
         real(dp) :: rmin_ij

         if (NTotal .eq. maxMol) return

!      Choose the type of molecule to be inserted
         if (nMolTypes .eq. 1) then
            nType = 1
         else
            nType = floor(nMolTypes*grnd() + 1d0)
         end if
         atmpSwapIn(nType) = atmpSwapIn(nType) + 1d0
         atmpInSize(NTotal) = atmpInSize(NTotal) + 1d0
         if (NPART(nType) .eq. NMAX(nTYPE)) then
            return
         end if

         NDiff = 0
         NDiff(nType) = 1
         call EBias_Insert_ChooseTarget(nType, nTarget, nTargType, nTargMol, ProbTarg_In)
         nTargIndx = MolArray(nTargType)%mol(nTargMol)%indx

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
         atmType1 = atomArray(nType, 1)
         atmType2 = atomArray(nTargType, 1)
         rmin_ij = dsqrt(r_min_tab(atmType1, atmType2))
         if (r .lt. rmin_ij) return

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
         if (distCriteria) then
            call Insert_NewNeiETable_Distance(nType, PairList, dETable, newNeiETable)
         else
            call Insert_NewNeiETable(nType, PairList, dETable, newNeiETable)
         end if
         call EBias_Insert_ReverseProbTarget(nTarget, nType, newNeiETable, ProbTarg_Out)
         call EBias_Insert_ReverseProbSel(nTarget, nType, dETable, ProbSel_Out)

!     Calculate the umbrella sampling bias.
         bIndx = getBiasIndex(NPart, NMAX)
         biasOld = NBias(bIndx)
         bIndx = getNewBiasIndex(NPart, NMAX, NDiff)
         biasNew = NBias(bIndx)
         bias_diff = biasNew - biasOld

!     Calculate acceptance probability and determine if the move is accepted or not
         genProbRatio = (ProbTarg_Out*ProbSel_Out*avbmc_vol*dble(nMolTypes)*gas_dens(nType))/ProbTarg_In
!      genProbRatio = (ProbTarg_Out * ProbSel_Out * avbmc_vol * gas_dens(nType)) / ProbTarg_In
         if (genProbRatio*exp(-beta*E_Inter + bias_diff) .gt. grnd()) then
            acptSwapIn(nType) = acptSwapIn(nType) + 1d0
            acptInSize(NTotal) = acptInSize(NTotal) + 1d0
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
      subroutine AVBMC_EBias_Out(E_T, acc_x)
         use SimParameters
         use Constants
         use E_Interface
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

         interface
            subroutine EBias_Remove_ChooseNeighbor(nTarget, biasArray, nSel, ProbTarget)
               integer, intent(in) :: nTarget
               real(dp), intent(in) :: biasArray(:)
               integer, intent(out) ::  nSel
               real(dp), intent(out) :: ProbTarget
            end subroutine
         end interface

         interface
            subroutine EBias_Remove_ReverseProbTarget(nTarget, nSel, nType, dE, ProbSel)
               implicit none
               integer, intent(in) :: nSel, nType
               integer, intent(in) :: nTarget
               real(dp), intent(in) :: dE(:)
               real(dp), intent(out) :: ProbSel
            end subroutine
         end interface

         real(dp), intent(inout) :: E_T
         real(dp), intent(inout) :: acc_x

         logical :: rejMove
         integer :: nTarget, nIndx, bIndx, iType
         integer :: nSel, nType, nMol, nTargMol, nTargType
         integer :: NDiff(1:nMolTypes)
         real(dp) :: grnd
         real(dp) :: genProbRatio
         real(dp) :: bias_diff
         real(dp) :: biasOld, biasNew
         real(dp) :: E_Inter, E_Intra
         real(dp) :: dETable(1:maxMol)
         real(dp) :: ProbTargOut, ProbSel, ProbTargIn
         real(dp) :: rx, ry, rz, dist
         real(dp) :: biasArray(1:nMolTypes)

         if (NTotal .eq. 1) return

!     Pick a Random Target Particle to Delete
!      nTarget = floor(NTotal*grnd() + 1d0)
!      call Get_MolIndex(nTarget, NPART, nTargType, nTargMol)
!      nTargIndx = MolArray(nTargType)%mol(nTargMol)%indx
         call Create_NeiETable
         call EBias_Remove_ChooseTarget(nTarget, nTargType, nTargMol, ProbTargOut)
         biasArray = 0d0
         bIndx = getBiasIndex(NPart, NMAX)
         biasOld = NBias(bIndx)
         do iType = 1, nMolTypes
            if (NPART(iType) .eq. 0) cycle
            NDiff = 0
            NDiff(iType) = -1
            if (NPART(iType) - 1 .lt. NMIN(iType)) cycle
            bIndx = getNewBiasIndex(NPart, NMAX, NDiff)
            biasNew = NBias(bIndx)
            bias_diff = biasNew - biasOld
            biasArray(iType) = bias_diff
         end do
         call EBias_Remove_ChooseNeighbor(nTarget, biasArray, nSel, ProbSel)
         nType = typeList(nSel)
         nMol = subIndxList(nSel)

         atmpSwapOut(nType) = atmpSwapOut(nType) + 1d0
         if (NPART(nType) .eq. NMIN(nType)) then
            return
         end if

         bias_diff = biasArray(nType)

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
         call EBias_Remove_ReverseProbTarget(nTarget, nSel, nType, dETable, ProbTargIn)
         genProbRatio = ProbTargIn/(ProbTargOut*ProbSel*dble(nMolTypes)*avbmc_vol*gas_dens(nType))

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
      subroutine EBias_Insert_ChooseTarget(nInsType, nTarget, nTargType, nMol, ProbSel)
         use SimParameters
         use Coords
         use EnergyTables
         implicit none
         integer, intent(in) :: nInsType
         integer, intent(out) :: nTarget, nTargType, nMol
         real(dp), intent(out) :: ProbSel

         integer :: i, iType
         integer :: cnt(1:nMolTypes)
         real(dp) :: avgE(1:nMolTypes)
         real(dp) :: ProbTable(1:maxMol)
         real(dp) :: grnd, norm
         real(dp) :: ranNum, sumInt

         ProbTable = 0d0
         avgE = 0d0
         cnt = 0
         do i = 1, maxMol
            if (isActive(i)) then
               iType = typeList(i)
               avgE(iType) = avgE(iType) + ETable(i)
               cnt(iType) = cnt(iType) + 1
            end if
         end do

         do iType = 1, nMolTypes
            if (cnt(iType) .ne. 0) then
               avgE(iType) = avgE(iType)/dble(cnt(iType))
            end if
         end do
         do i = 1, maxMol
            if (isActive(i)) then
               iType = typeList(i)
               ProbTable(i) = exp(biasAlpha(nInsType, iType)*(ETable(i) - avgE(iType)))
            end if
         end do

         norm = sum(ProbTable)
         ranNum = norm*grnd()

         sumInt = 0d0
         nTarget = 0
         do while (sumInt .lt. ranNum)
            nTarget = nTarget + 1
            sumInt = sumInt + ProbTable(nTarget)
         end do
         nTargType = typeList(nTarget)
         ProbSel = ProbTable(nTarget)/norm

         nMol = 0
         do i = 1, nTargType - 1
            nMol = nMol + NMAX(i)
         end do
         nMol = nTarget - nMol

      end subroutine
!=================================================================================
      subroutine EBias_Insert_ReverseProbTarget(nTarget, nType, newNeiETable, ProbRev)
         use SimParameters
         use Coords
         use EnergyTables
         implicit none
         integer, intent(in) :: nTarget, nType
         real(dp), intent(in) :: newNeiETable(:)
         real(dp), intent(out) :: ProbRev

         integer :: i, nIndx
         integer :: cnt(1:nMolTypes)
         real(dp) :: ProbTable(1:maxMol)
         real(dp) :: norm, EMax

         nIndx = molArray(nType)%mol(NPART(nType) + 1)%indx
         ProbTable = 0d0
         EMax = newNeiETable(nTarget)
!      EMax = 0d0
         cnt = 0
         do i = 1, maxMol
            if (isActive(i) .or. (i .eq. nIndx)) then
               ProbTable(i) = exp(beta*(newNeiETable(i) - Emax))
            end if
         end do
!      ProbTable(nIndx) = exp(beta*(newNeiETable(nIndx)-Emax))

         norm = sum(ProbTable)
         ProbRev = ProbTable(nTarget)/norm

      end subroutine
!=================================================================================
      subroutine EBias_Insert_ReverseProbSel(nTarget, nType, dE, ProbRev)
         use SimParameters
         use Coords
         use EnergyTables
         use UmbrellaFunctions
         implicit none
         integer, intent(in) :: nTarget, nType
         real(dp), intent(in) :: dE(:)
         real(dp), intent(out) :: ProbRev

         integer :: i, iType, nIndx, bIndx
         integer :: NDiff(1:nMolTypes)
         real(dp) :: biasArray(1:nMolTypes)
         real(dp) :: norm, bias_diff, biasOld, biasNew

!      biasArray = 0d0
!      NDiff = 0
!      NDiff(nType) = +1
!      bIndx = getNewBiasIndex(NPart,NMAX, NDiff)
!      biasOld = NBias(bIndx)
!      do iType = 1, nMolTypes
!        NDiff = 0
!        NDiff(nType) = +1
!        NDiff(iType) = NDiff(iType) - 1
!        if(NPART(iType) - 1 .lt. NMIN(iType)) cycle
!        bIndx = getNewBiasIndex(NPart,NMAX, NDiff)
!        biasNew = NBias(bIndx)
!        bias_diff = biasNew - biasOld
!        biasArray(iType) = bias_diff
!      enddo

         nIndx = molArray(nType)%mol(NPART(nType) + 1)%indx
         norm = 0d0
         do i = 1, maxMol
            if (.not. isActive(i)) then
               cycle
            end if
            if (NeighborList(i, nTarget)) then
               if (i .ne. nTarget) then
!            iType = typeList(i)
!            norm = norm + exp(beta*(ETable(i)+dE(i)-dE(nIndx)) + biasArray(iType))
                  norm = norm + exp(beta*(ETable(i) + dE(i) - dE(nIndx)))
               end if
            end if
         end do
!      norm = norm + exp(biasArray(nType))
         norm = norm + 1d0
!      ProbRev = exp(biasArray(nType))/norm
         ProbRev = 1d0/norm

      end subroutine
!=================================================================================
      subroutine EBias_Remove_ChooseTarget(nTarget, nType, nMol, ProbTarget)
         use SimParameters
         use Coords
         use EnergyTables
         implicit none
         integer, intent(out) :: nTarget, nType, nMol
         real(dp), intent(out) :: ProbTarget

         integer :: i
         real(dp) :: ProbTable(1:maxMol)
         real(dp) :: scaleE
         real(dp) :: grnd, norm
         real(dp) :: ranNum, sumInt

         ProbTable = 0d0
         scaleE = 0d0
         do i = 1, maxMol
            if (isActive(i)) then
               scaleE = NeiETable(i)
               exit
            end if
         end do
!      write(6,*) scaleE
         do i = 1, maxMol
            if (isActive(i)) then
               ProbTable(i) = exp(beta*(NeiETable(i) - scaleE))
            end if
         end do

         norm = sum(ProbTable)
         ranNum = norm*grnd()

         sumInt = 0d0
         nTarget = 0
         do while (sumInt .lt. ranNum)
            nTarget = nTarget + 1
            sumInt = sumInt + ProbTable(nTarget)
         end do

         nType = typeList(nTarget)
         ProbTarget = ProbTable(nTarget)/norm

         nMol = 0
         do i = 1, nType - 1
            nMol = nMol + NMAX(i)
         end do
         nMol = nTarget - nMol

      end subroutine
!=================================================================================
      subroutine EBias_Remove_ChooseNeighbor(nTarget, biasArray, nSel, ProbTarget)
         use SimParameters
         use Coords
         use EnergyTables
         implicit none
         integer, intent(in) :: nTarget
         real(dp), intent(in) :: biasArray(:)
         integer, intent(out) ::  nSel
         real(dp), intent(out) :: ProbTarget

         integer :: i, iType
         real(dp) :: ProbTable(1:maxMol)
         real(dp) :: grnd, norm
         real(dp) :: ranNum, sumInt
         real(dp) :: dummy

         ProbTable = 0d0
         dummy = biasArray(1)
         do i = 1, maxMol
            if (.not. isActive(i)) then
               cycle
            end if
            if (NeighborList(i, nTarget)) then
               if (i .ne. nTarget) then
!            iType = typeList(i)
!            ProbTable(i) = exp(beta * ETable(i) + biasArray(iType))
                  ProbTable(i) = exp(beta*ETable(i))
               end if
            end if
         end do

         norm = sum(ProbTable)
         ranNum = norm*grnd()

         sumInt = 0d0
         nSel = 0
         do while (sumInt .lt. ranNum)
            nSel = nSel + 1
            sumInt = sumInt + ProbTable(nSel)
         end do

         ProbTarget = ProbTable(nSel)/norm

      end subroutine
!=================================================================================
      subroutine EBias_Remove_ReverseProbTarget(nTarget, nSel, nType, dE, ProbSel)
         use SimParameters
         use Coords
         use EnergyTables
         implicit none
         integer, intent(in) :: nSel, nType
         integer, intent(in) :: nTarget
         real(dp), intent(in) :: dE(:)
         real(dp), intent(out) :: ProbSel

         integer :: i, iType
         integer :: cnt(1:nMolTypes)

         real(dp) :: avgE(1:nMolTypes)
         real(dp) :: ProbTable(1:maxMol)
         real(dp) :: grnd, norm

         if (NTotal .eq. 2) then
            ProbSel = 1d0
            return
         end if

         ProbTable = 0d0
         avgE = 0d0
         cnt = 0
         do i = 1, maxMol
            if (isActive(i)) then
               if (i .ne. nSel) then
                  iType = typeList(i)
                  avgE(iType) = avgE(iType) + ETable(i) - dE(i)
                  cnt(iType) = cnt(iType) + 1

               end if
            end if
         end do

         do iType = 1, nMolTypes
            if (cnt(iType) .ne. 0) then
               avgE(iType) = avgE(iType)/dble(cnt(iType))
            end if
         end do

         do i = 1, maxMol
            if (isActive(i)) then
               if (i .ne. nSel) then
                  iType = typeList(i)
                  ProbTable(i) = exp(biasAlpha(nType, iType)*((ETable(i) - dE(i)) - avgE(iType)))
               end if
            end if
         end do

         norm = sum(ProbTable)
         ProbSel = ProbTable(nTarget)/norm

      end subroutine
!=================================================================================
