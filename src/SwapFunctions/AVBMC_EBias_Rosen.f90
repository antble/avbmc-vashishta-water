!****** Energy Biased Aggregation-Volume-Bias Monte Carlo (AVBMC) algorithm *******
!   This file contains the nessisary functions to impliment the energy biased swap
!   move for cluster simulations.
!===================================================================================
module AVBMC_Module
   use VarPrecision
   real(dp), allocatable :: swapProb(:)
   logical :: detailbalance_check = .false. 
contains
!===================================================================================
   subroutine AVBMC(E_T, acc_x, atmp_x)
      use SimParameters
      use Constants
      implicit none
      real(dp), intent(inout) :: E_T, atmp_x, acc_x
      real(dp) :: rndnum, grnd
      logical :: accept

      prevMoveAccepted = .false.

      if (detailbalance_check) then 
         call AVBMC_CheckReversibility(E_T, maxMol, accept)
         return 
      end if 

      if (grnd() .lt. 0.5d0) then
         call AVBMC_EBias_Rosen_In(E_T, maxMol, acc_x, atmp_x)
      else
         call AVBMC_EBias_Rosen_Out(E_T, maxMol, acc_x, atmp_x)
      end if
      
   end subroutine
!===================================================================================
   subroutine AVBMC_CheckReversibility(E_T, arrayMax, accept)
      use AVBMC_CBMC
      use AcceptRates
      use CBMC_Utility
      use CBMC_Variables
      use Constants
      use Coords
      use SimParameters
      use ForceField
      use EnergyPointers, only: SwapIn_ECalc, SwapOut_ECalc, Update_SubEnergies, Update_SubEnergies_Vashishta
      use ForceField
      use IndexingFunctions
      use EnergyCriteria
      use DistanceCriteria
      use EnergyTables
      use NeighborTable
      use SwapBoundary
      implicit none

      integer, intent(in) :: arrayMax
      real(dp), intent(inout) :: E_T
      logical, intent(out) :: accept 
      logical :: rejMove, rejMove2
      logical :: isIncluded(1:arrayMax)
      integer :: i, nTargType, nTargMol, nTargIndx, nTarget
      integer :: nType, nIndx,  nSel, nMol, targid, forcetargid, forceid
      integer :: NDiff(1:nMolTypes)
      real(dp) :: grnd
      real(dp) :: genProbRatio, rosenRatio, rosenI, rosenD
      real(dp) :: E_Inter1, E_Inter2, E_Intra, biasDiff
      real(dp) :: PairList(1:arrayMax)
      real(dp) :: dETable(1:arrayMax)
      real(dp) :: newNeiETable(1:arrayMax)
      real(dp) :: ProbTarg_In, ProbTarg_Out, ProbSel_Out  ! insertion
      real(dp) :: ProbTargOut, ProbSel, ProbTargIn     ! deletion
      real(dp) :: Boltzterm, P12, P21, A12, A21, detailBalance, E_1, E_2
      real(dp) :: acc_x, atmp_x

      print *,"====================================================="
      print *, "Checking reversibility...INSERTION"
      E_1 = E_T
      call AVBMC_EBias_Rosen_In(E_T, arrayMax, acc_x, atmp_x, targid, A12, E_Inter1, accept, rosenI)
      E_2 = E_T
      print *, "INSERTION Accept =", accept
      print *, "TargetID=",targid
      if (accept) then 
         print *, "E_1 + E_Inter1 = E_2,  log(A12), rosenI", E_1, E_Inter1, E_2,  log(A12), rosenI
      else 
         print *, "E_1,  E_2, E_Inter1,  log(A12), rosenI", E_1, E_2, E_Inter1, log(A12), rosenI
      end if 

      if (.not. accept) then 
         print *, "returning ..."
         print *,"====================================================="
         return 
      end if
      print *, "Current number of particles", NPART(1)
      print *, "Checking reversibility...DELETION"
      call AVBMC_EBias_Rosen_Out(E_T, arrayMax, acc_x, atmp_x, targid, forceid, A21, E_Inter2, accept, rosenD)
      if (accept) then 
         print *, "E_2 + E_Inter2 = E_1,  log(A21), rosenD", E_2, E_Inter2, E_1, log(A21), rosenD
      else 
         print *, "E_T, E_Inter2, E_1, log(A21), rosenD", E_T, E_Inter2, E_1,  log(A21), rosenD
      end if 
      
      print *, "DELETION Accept ?", accept
      print *, "Target id=",targid
      print *, "forceID=",forceid
      if (.not. accept) then 
         print *, "returning ..."
         print *,"====================================================="
         return 
      end if
      print *,"====================================================="
      if( ((A12 == 0E0_dp) .or. (A21 == 0E0_dp)) .and. (A21 /= A12) ) then
         print *, "Zero Probability observed in move that should be reverisbile"
         print *, "This implies a calculation error in the detailed balance."
         print *, "A12, A21:", A12, A21
         error stop
      end if

      P12 = log(A12) - beta*E_Inter1
      print *, "P12=", P12
      if(P12 > 0E0_dp) P12 = 0E0_dp
      P21 = log(A21) - beta*E_Inter2
      print *, "P21=",P21
      if(P21 > 0E0_dp) P21 = 0E0_dp

      print *, "log(rosenI)-log(rosenD) |", log(rosenI)-log(rosenD), log(rosenI), log(rosenD)
      detailBalance = log(A21) + P12 - P21 - beta*E_Inter2 + log(rosenI) - log(rosenD)
      print *, "E_inter", E_Inter2, E_Inter1
      print *, "A12, A21, A12*A21, log(A21)", A12, A21, A12*A21 , log(A21)
      print *, "ln(P12), ln(P21):", P12, P21
      print *, "Detailed balance",  detailBalance, log(A21), P12, -P21, -beta*E_Inter1

      if (abs(detailBalance) > 1e-6_dp) then 
         stop "Violation of Detailed Balance Detected!" 
      end if 
      ! reversibility_count = reversibility_count + 1d0
   end subroutine 

   subroutine AVBMC_EBias_Rosen_In(E_T, arrayMax, acc_x, atmp_x, targid, A12, dE12, accept, rosenI)
      use AcceptRates
      use AVBMC_RejectionVar
      use AVBMC_CBMC
      use CBMC_Utility
      use CBMC_Variables
      use Constants
      use Coords
      use SimParameters
      use ForceField
      use EnergyPointers, only: SwapIn_ECalc, SwapOut_ECalc, Update_SubEnergies, Quick_Nei_ECalc, Update_SubEnergies_Vashishta
      use UmbrellaFunctions
      use ForceField
      use IndexingFunctions
      use EnergyCriteria
      use DistanceCriteria
      use EnergyTables
      use NeighborTable
      use SwapBoundary
      use Pressure_LJ_Electro, only: NewMol_PressCalc_Inter
      use PairStorage, only: UpdateDistArray, PrintDistArray
      use UmbrellaSamplingNew, only: GetUmbrellaBias_SwapIn
      
      implicit none


      integer, intent(in) :: arrayMax
      real(dp), intent(inout) :: E_T
      real(dp), intent(inout) :: acc_x, atmp_x
      integer, intent(out), optional :: targid 
      real(dp), intent(out), optional :: A12, dE12, rosenI
      logical, intent(out), optional :: accept
      logical :: rejMove

      logical :: isIncluded(1:arrayMax)
      integer :: NDiff(1:nMolTypes)
      integer :: i, nTargType, nTargMol, nTargIndx, nTarget
      integer :: nType, nIndx
      real(dp) :: grnd
      real(dp) :: genProbRatio, rosenRatio
      real(dp) :: E_Inter, E_Intra, biasDiff
      real(dp) :: PairList(1:arrayMax)
      real(dp) :: dETable(1:arrayMax)
      real(dp) :: newNeiETable(1:arrayMax)
      real(dp) :: ProbTarg_In, ProbTarg_Out, ProbSel_Out
      real(dp) :: Boltzterm, E_insertion, E_deletion
      real(dp) :: sumInt, ranNum
      character(len=100) :: fl_name

      integer :: nSel, nMol

      if(present(accept)) then 
         accept=.false. 
      end if

      if (NTotal .eq. maxMol) then
         boundaryRej = boundaryRej + 1d0
         totalRej = totalRej + 1d0
         return
      end if

!      Choose the type of molecule to be inserted
      if (nMolTypes .eq. 1) then
         nType = 1
      else
!        nType = floor(nMolTypes*grnd() + 1d0)
         ranNum = grnd()
         sumInt = swapProb(1)
         nType = 1
         do while (sumInt .lt. ranNum)
            nType = nType + 1
            sumInt = sumInt + swapProb(nType)
         end do
      end if

      NDiff = 0
      NDiff(nType) = 1
      rejMove = boundaryFunction(NPART, NDiff)
      if (rejMove) then
         boundaryRej = boundaryRej + 1d0
         totalRej = totalRej + 1d0
         return
      end if

 
      atmp_x = atmp_x + 1d0
      atmpSwapIn(nType) = atmpSwapIn(nType) + 1d0
      atmpInSize(NTotal) = atmpInSize(NTotal) + 1d0

      call EBias_Insert_ChooseTarget(nType, nTarget, nTargType, nTargMol, ProbTarg_In)
      nTargIndx = MolArray(nTargType)%mol(nTargMol)%indx
      if (present(targid)) then 
         targid = nTarget ! output
      end if
!      Generate the configuration for the newly inserted molecule
      nIndx = MolArray(nType)%mol(NPART(nType) + 1)%indx
      call Rosen_CreateSubset(nTarget, isIncluded)

      select case (regrowType(nType))
      case (0)
         call Ridgid_RosenConfigGen(nType, nIndx, nTarget, nTargType, isIncluded, rosenRatio, rejMove)
      case (1)
         call Simple_RosenConfigGen(nType, nIndx, nTarget, nTargType, isIncluded, rosenRatio, rejMove)
      case (2)
         call StraightChain_RosenConfigGen(nType, nIndx, nTarget, nTargType, isIncluded, rosenRatio, rejMove)
      case (3)
         call BranchedMol_RosenConfigGen(nType, nIndx, nTarget, nTargType, isIncluded, rosenRatio, rejMove)
      case (4)
         call Simple_RosenConfigGen_Vashishta(nType, nIndx, nTarget, nTargType, isIncluded, rosenRatio, rejMove)
      case (5)
         call Ridgid_RosenConfigGen_Vashishta(nType, nIndx, nTarget, nTargType, isIncluded, rosenRatio, rejMove)
      case default
         write (*, *) "Error! EBias can not regrow a molecule of regrow type:", regrowType(nType)
         stop
      end select
      if (rejMove) then
         totalRej = totalRej + 1d0
         ovrlapRej = ovrlapRej + 1d0
         return
      end if

      ! call DEBUG_Output_NewConfig

!      Perform a check to see if the cluster criteria is statisfied or not.
      if (.not. distCriteria) then
         rejMove = .false.
         call Quick_Nei_ECalc(nTargType, nTargMol, rejMove)
         if (rejMove) then
            totalRej = totalRej + 1d0
            critriaRej = critriaRej + 1d0
            clusterCritRej = clusterCritRej + 1d0
            return
         end if
      end if

!     Calculate the umbrella sampling bias.
      NPART_new = NPART + NDiff
      NTotal_New = NTotal + 1
      call GetUmbrellaBias_SwapIn(biasDiff, rejMove)
      if (rejMove) then
         if (present(targid)) then 
            accept = .false.
         end if
         boundaryRej = boundaryRej + 1d0
         totalRej = totalRej + 1d0
         return
      end if
!      Calculate the Energy Difference Associated with the move
      E_Inter = 0d0
      E_Intra = 0d0
      call SwapIn_ECalc(E_Inter, E_Intra, PairList, dETable, rejMove)
      ! print *, "SwapIn_ECalc E_Intra", E_Intra

      if (present(targid)) then 
      print *, "SwapIn_ECalc E_Inter, E_Intra, biasDiff:=", E_Inter, E_Intra, biasDiff
         dE12 = E_Inter 
         rosenI = rosenRatio
      end if
      if (rejMove) then
         totalRej = totalRej + 1d0
         ovrlapRej = ovrlapRej + 1d0
         return
      end if

!     Determine the reverse probability of this move.
      if (distCriteria) then
         call Insert_NewNeiETable_Distance_V2(nType, PairList, dETable, newNeiETable)
      else
         call Insert_NewNeiETable(nType, PairList, dETable, newNeiETable)
      end if
      call EBias_Insert_ReverseProbTarget(nTarget, nType, newNeiETable, ProbTarg_Out)
      call EBias_Insert_ReverseProbSel(nTarget, nType, dETable, ProbSel_Out)



!     Calculate acceptance probability and determine if the move is accepted or not
      genProbRatio = (ProbTarg_Out*ProbSel_Out*avbmc_vol*gas_dens(nType))/(ProbTarg_In*rosenRatio)
      if (present(targid)) then 
         print *, "ProbTarg_Out=",ProbTarg_Out
         print *, "ProbSel_Out=",ProbSel_Out
         print *, "avbmc_vol=",avbmc_vol
         print *, "gas_dens(nType)=",gas_dens(nType)
         print *, "ProbTarg_In=", ProbTarg_In
         print *, "rosenRatio",rosenRatio
         print *, "A12", genProbRatio
         A12 = genProbRatio
      end if

      select case (trim(adjustl(ForceFieldName)))
      case('Vashishta')
         Boltzterm = exp(-beta*(E_Inter - E_Intra) + biasDiff)
      case default
         Boltzterm = exp(-beta*(E_Inter) + biasDiff)
      end select

      if (genProbRatio*Boltzterm .gt. grnd()) then
         if (present(targid)) then 
            accept = .true.
         end if

         if (calcPressure) then
            call NewMol_PressCalc_Inter(P_Diff)
            pressure = pressure + P_Diff
         end if
!         call PrintDistArray
         acptSwapIn(nType) = acptSwapIn(nType) + 1d0
         acptInSize(NTotal) = acptInSize(NTotal) + 1d0
         do i = 1, nAtoms(nType)
            molArray(nType)%mol(NPART(nType) + 1)%x(i) = newMol%x(i)
            molArray(nType)%mol(NPART(nType) + 1)%y(i) = newMol%y(i)
            molArray(nType)%mol(NPART(nType) + 1)%z(i) = newMol%z(i)
         end do
         
         select case (trim(adjustl(ForceFieldName)))
         case('Vashishta')
            E_T = E_T + E_Inter
         case default
            E_T = E_T + E_Inter + E_Intra
         end select
         
         acc_x = acc_x + 1d0
         isActive(nIndx) = .true.
         if (distCriteria) then
            call NeighborUpdate_SwapIn_Distance(PairList, nType)
         else
            call NeighborUpdate(PairList, nIndx)
         end if
         call UpdateDistArray
         NTotal = NTotal + 1
         ETable = ETable + dETable
         NPART(nType) = NPART(nType) + 1

         select case (trim(adjustl(ForceFieldName)))
         case('Vashishta')
            call Update_SubEnergies_Vashishta
         case default
            call Update_SubEnergies
         end select
         !call Update_SubEnergies
         prevMoveAccepted = .true.
      else
         totalRej = totalRej + 1d0
         dbalRej = dbalRej + 1d0
      end if
   end subroutine
!===================================================================================
   subroutine AVBMC_EBias_Rosen_Out(E_T, arrayMax, acc_x, atmp_x, forcetargid, forceid, A21, dE21, accept, rosenD)
      use AVBMC_CBMC
      use AVBMC_RejectionVar
      use SimParameters
      use Constants
!      use E_Interface
      use EnergyPointers, only: SwapOut_ECalc, Update_SubEnergies, Update_SubEnergies_Vashishta
      use Coords
      use UmbrellaFunctions
      use ForceField
      use IndexingFunctions
      use EnergyCriteria
      use DistanceCriteria
      use EnergyTables
      use AcceptRates
      use UmbrellaFunctions
      use CBMC_Variables
      use NeighborTable
      use SwapBoundary
      use PairStorage, only: UpdateDistArray_SwapOut
      use UmbrellaSamplingNew, only: GetUmbrellaBias_SwapOut
      use Pressure_LJ_Electro, only: Mol_PressCalc_Inter
      implicit none

      integer, intent(in) :: arrayMax
      real(dp), intent(inout) :: E_T
      real(dp), intent(inout) :: acc_x, atmp_x
      integer, intent(inout), optional :: forcetargid, forceid 
      real(dp), intent(out), optional :: A21, dE21, rosenD
      logical, intent(inout), optional :: accept
      logical :: rejMove
      integer :: nTarget, nIndx
      integer :: nSel, nType, nMol, nTargMol, nTargType
      integer :: NDiff(1:nMolTypes)
      real(dp) :: grnd, Boltzterm 
      real(dp) :: genProbRatio
      real(dp) :: biasDiff
      real(dp) :: E_Inter, E_Intra
      real(dp) :: dETable(1:arrayMax)
      real(dp) :: ProbTargOut, ProbSel, ProbTargIn
      real(dp) :: rx, ry, rz, dist, rosenRatio
      real(dp) :: ranNum, sumInt 
      real(dp) :: dist1, dist2, rndNum

      
      if (present(forcetargid)) then 
            accept = .false.
         end if
      if (NTotal .eq. 1) then
         boundaryRej_out = boundaryRej_out + 1d0
         totalRej_out = totalRej_out + 1d0
         return
      end if


!     Pick a Random Particle Type
      if (nMolTypes .eq. 1) then
         nType = 1
      else
!        nType = floor(nMolTypes*grnd() + 1d0)
         ranNum = grnd()
         sumInt = swapProb(1)
         nType = 1
         do while (sumInt .lt. ranNum)
            nType = nType + 1
            sumInt = sumInt + swapProb(nType)
         end do
      end if

      NDiff = 0
      NDiff(nType) = -1
      rejMove = boundaryFunction(NPART, NDiff)
      if (rejMove) then
         boundaryRej_out = boundaryRej_out + 1d0
         totalRej_out = totalRej_out + 1d0

         print *, "rejMove", rejMove
         stop
         return
      end if
      
      call Create_NeiETable(nType)
      if(present(forcetargid)) then
         call EBias_Remove_ChooseTarget(nTarget, nTargType, nTargMol, ProbTargOut, forcetargid)
      else 
         call EBias_Remove_ChooseTarget(nTarget, nTargType, nTargMol, ProbTargOut)
      end if
      
      !! call EBias_Remove_ChooseTarget(nTarget, nTargType, nTargMol, ProbTargOut)
      !! call EBias_Remove_ChooseNeighbor(nTarget, nType, nSel, ProbSel)
      if(present(forcetargid)) then
         call EBias_Remove_ChooseNeighbor(nTarget, nType, nSel, ProbSel, NPART(1))
         nSel = NPART(1)
         forceid = nSel
      else 
         call EBias_Remove_ChooseNeighbor(nTarget, nType, nSel, ProbSel)
      end if
!      nType = typeList(nSel)
      nMol = subIndxList(nSel)
      atmp_x = atmp_x + 1d0
      atmpSwapOut(nType) = atmpSwapOut(nType) + 1d0

!     Check to see that the appropriate atoms are within the insertion distance
!     in order to ensure the move is reversible. If not reject the move since
!     the reverse probility is equal to 0.
      if (.not. distCriteria) then
         rx = molArray(nTargType)%mol(nTargMol)%x(1) - molArray(nType)%mol(nMol)%x(1)
         ry = molArray(nTargType)%mol(nTargMol)%y(1) - molArray(nType)%mol(nMol)%y(1)
         rz = molArray(nTargType)%mol(nTargMol)%z(1) - molArray(nType)%mol(nMol)%z(1)
         dist = rx*rx + ry*ry + rz*rz
         if (dist .gt. Dist_Critr_sq) then
            critriaRej_out = critriaRej_out + 1d0
            totalRej_out = totalRej_out + 1d0
            return
         end if
      end if
      
      !============================================================================
      ! The purpose of this process is to prevent a water with dissociated hydrogen
      ! from being deleted
      !============================================================================
      ! check if the molecule is OH using an rmax=1.4 bonding distance definition.
      if (trim(adjustl(ForceFieldName)) == 'Vashishta') then 
         rx = molArray(nTargType)%mol(nMol)%x(1) - molArray(nType)%mol(nMol)%x(2)
         ry = molArray(nTargType)%mol(nMol)%y(1) - molArray(nType)%mol(nMol)%y(2)
         rz = molArray(nTargType)%mol(nMol)%z(1) - molArray(nType)%mol(nMol)%z(2)
         dist1 = sqrt(rx*rx + ry*ry + rz*rz)

         rx = molArray(nTargType)%mol(nMol)%x(1) - molArray(nType)%mol(nMol)%x(3)
         ry = molArray(nTargType)%mol(nMol)%y(1) - molArray(nType)%mol(nMol)%y(3)
         rz = molArray(nTargType)%mol(nMol)%z(1) - molArray(nType)%mol(nMol)%z(3)
         dist2 = sqrt(rx*rx + ry*ry + rz*rz)
         !! Early Tests.
         ! if (dist2 > 1.4 .or. dist1 > 1.4) then
         ! if (dist2 > 1.3 .or. dist1 > 1.3) then
         ! if (dist2 > 1.1 .or. dist1 > 1.1) then
         if (dist2 > 2.0 .or. dist1 > 2.0) then
            rejMove = .true. 
            OHrej_out = OHrej_out + 1d0
            totalRej_out = totalRej_out + 1d0
            return
         end if 
      end if
      
!     Calculate the umbrella sampling bias.
!      NDiff = 0
!      NDiff(nType) = -1
      NPART_new = NPART + NDiff
      NTotal_New = NTotal - 1
      call GetUmbrellaBias_SwapOut(nType, nMol, biasDiff, rejMove)
      if (rejMove) then
         boundaryRej_out = boundaryRej_out + 1d0
         totalRej_out = totalRej_out + 1d0
         return
      end if
      

!     Check to see if the deletion of the particle will break the cluster
      rejMove = .false.
      call SwapOut_EnergyCriteria(nSel, rejMove)
      
      if (rejMove) then
         critriaRej_out = critriaRej_out + 1d0
         totalRej_out = totalRej_out + 1d0
         return
      end if

      select case (regrowType(nType))
      case (0)
         call Ridgid_RosenConfigGen_Reverse(nType, nMol, nTarget, nTargType, rosenRatio)
      case (1)
         call Simple_RosenConfigGen_Reverse(nType, nMol, nTarget, nTargType, rosenRatio)
      case (2)
         call StraightChain_RosenConfigGen_Reverse(nType, nMol, nTarget, nTargType, rosenRatio)
      case (3)
         call BranchedMol_RosenConfigGen_Reverse(nType, nMol, nTarget, nTargType, rosenRatio)
      case (4)
         call Simple_RosenConfigGen_Reverse_Vashishta(nType, nMol, nTarget, nTargType, rosenRatio)
      case (5)
         call Ridgid_RosenConfigGen_Reverse_Vashishta(nType, nMol, nTarget, nTargType, rosenRatio)
      case default
         write (*, *) "Error! EBias can not regrow a molecule of regrow type:", nType
         stop
      end select
!      Calculate the Energy Difference Associated with the move.
      E_Inter = 0d0
      E_Intra = 0d0
!      call SwapOut_EnergyCalc(E_Inter, E_Intra, nType, nMol, dETable)
      call SwapOut_ECalc(E_Inter, E_Intra, nType, nMol, dETable)
      ! print *, "SwapOut_ECalc E_Intra", E_Intra
      if (present(forcetargid)) then 
         print *, "SwapOut_ECalc E_Inter, E_Intra, biasDiff:=", E_Inter, E_Intra, biasDiff
         dE21 = E_Inter
         rosenD = rosenRatio
      end if
      call EBias_Remove_ReverseProbTarget(nTarget, nSel, nType, dETable, ProbTargIn)

!      genProbRatio = (ProbTargIn * rosenRatio) / (ProbTargOut * ProbSel * dble(nMolTypes) * avbmc_vol * gas_dens(nType))
      genProbRatio = (ProbTargIn*rosenRatio)/(ProbTargOut*ProbSel*avbmc_vol*gas_dens(nType))
      select case (trim(adjustl(ForceFieldName)))
      case('Vashishta')
         Boltzterm = exp(-beta*(E_Inter + E_Intra) + biasDiff)
      case default
         Boltzterm = exp(-beta*(E_Inter) + biasDiff)
      end select
      if (present(forcetargid)) then 
         print *, "ProbTargIn*rosenRatio",ProbTargIn
         print *, "rosenRatio",rosenRatio
         print *, "ProbTargOut",ProbTargOut
         print *, "ProbSel",ProbSel
         print *, "avbmc_vol",avbmc_vol
         print *, "gas_dens(nType)",gas_dens(nType)
         print *, "A21=genProbRatio", genProbRatio
         print *, "Boltzmann term:", Boltzterm
         A21 = genProbRatio
      end if
!      Calculate Acceptance and determine if the move is accepted or not
      ! if (genProbRatio*exp(-beta*(E_Inter + E_Intra) + biasDiff) .gt. grnd()) then
      if (genProbRatio*Boltzterm .gt. grnd()) then
         if (present(forcetargid)) then 
            accept = .true.
         end if

         if (calcPressure) then
            call Mol_PressCalc_Inter(nType, nMol, P_Diff)
            pressure = pressure - P_Diff
         end if
         acptSwapOut(nType) = acptSwapOut(nType) + 1d0
         molArray(nType)%mol(nMol)%x(1:nAtoms(nType)) = molArray(nType)%mol(NPART(nType))%x(1:nAtoms(nType))
         molArray(nType)%mol(nMol)%y(1:nAtoms(nType)) = molArray(nType)%mol(NPART(nType))%y(1:nAtoms(nType))
         molArray(nType)%mol(nMol)%z(1:nAtoms(nType)) = molArray(nType)%mol(NPART(nType))%z(1:nAtoms(nType))
         
         select case (trim(adjustl(ForceFieldName)))
         case('Vashishta')
            E_T = E_T + E_Inter
         case default
            E_T = E_T + E_Inter + E_Intra
         end select
         
         nIndx = molArray(nType)%mol(nMol)%indx
         call NeighborUpdate_Delete(nIndx, molArray(nType)%mol(NPART(nType))%indx)

         call UpdateDistArray_SwapOut(nType, nMol)
         isActive(molArray(nType)%mol(NPART(nType))%indx) = .false.
         ETable = ETable - dETable
         ETable(nIndx) = ETable(molArray(nType)%mol(NPART(nType))%indx)
         ETable(molArray(nType)%mol(NPART(nType))%indx) = 0d0
         NPART(nType) = NPART(nType) - 1
         NTotal = NTotal - 1
         acc_x = acc_x + 1d0

         select case (trim(adjustl(ForceFieldName)))
         case('Vashishta')
            call Update_SubEnergies_Vashishta
         case default
            call Update_SubEnergies
         end select
         prevMoveAccepted = .true.
!         call DEBUG_Output_NeighborList
      else
         dbalRej_out = dbalRej_out + 1d0
         totalRej_out = totalRej_out + 1d0
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
      real(dp) :: ProbTable(1:maxMol)
      real(dp) :: norm, EMax

      nIndx = molArray(nType)%mol(NPART(nType) + 1)%indx
      ProbTable = 0E0_dp
      EMax = -huge(dp)
      do i = 1, maxMol
         if (neiCount(i) .gt. 0) then
            if (newNeiETable(i) .gt. EMax) then
               EMax = newNeiETable(i)
            end if
         end if
      end do
      do i = 1, maxMol
         if (neiCount(i) .gt. 0) then
            if (beta*(newNeiETable(i) - Emax) .gt. -30d0) then
               ProbTable(i) = exp(beta*(newNeiETable(i) - Emax))
            end if
         end if
      end do

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

      integer :: iMol, iIndx, nIndx
      real(dp) :: norm

      nIndx = molArray(nType)%mol(NPART(nType) + 1)%indx
      norm = 0d0
      do iMol = 1, NPART(nType)
         iIndx = molArray(nType)%mol(iMol)%indx
         if (NeighborList(iIndx, nTarget)) then
            norm = norm + exp(beta*(ETable(iIndx) + dE(iIndx) - dE(nIndx)))
         end if
      end do
      norm = norm + 1E0
      ProbRev = 1E0/norm
   end subroutine
!=================================================================================
   subroutine EBias_Remove_ChooseTarget(nTarget, nTargType, nTargMol, ProbTarget, forcetargid)
      use Coords
      use EnergyTables
      use SimParameters
      implicit none
      integer, intent(out) :: nTarget, nTargType, nTargMol
      real(dp), intent(out) :: ProbTarget
      integer, intent(in), optional :: forcetargid
      integer :: i
      real(dp) :: ProbTable(1:maxMol)
      real(dp) :: grnd, norm
      real(dp) :: ranNum, sumInt
      
      ProbTable = 0d0

      do i = 1, maxMol
!        if(isActive(i)) then
         if (neiCount(i) .gt. 0) then
            ProbTable(i) = exp(beta*NeiETable(i))
         end if
      end do
      norm = sum(ProbTable)
      ranNum = norm*grnd()

      sumInt = ProbTable(1)
      nTarget = 1
      do while (sumInt .lt. ranNum)
         nTarget = nTarget + 1
         sumInt = sumInt + ProbTable(nTarget)
      end do

      if(present(forcetargid)) then 
         nTarget = forcetargid 
      end if 

      nTargType = typeList(nTarget)
      ProbTarget = ProbTable(nTarget)/norm

      nTargMol = 0
      do i = 1, nTargType - 1
         nTargMol = nTargMol + NMAX(i)
      end do
      nTargMol = nTarget - nTargMol
      ! print *, "EBias_Remove_ChooseTarget", nTarget, forcetargid
   end subroutine
!=================================================================================
   subroutine EBias_Remove_ChooseNeighbor(nTarget, nType, nSel, ProbTarget, forceID)
      use SimParameters
      use Coords
      use EnergyTables
      implicit none
      integer, intent(in) :: nTarget, nType
      integer, intent(out) ::  nSel
      real(dp), intent(out) :: ProbTarget
      integer, intent(in), optional :: forceID

      integer :: iIndx, iMol
      real(dp) :: ProbTable(1:maxMol)
      real(dp) :: grnd, norm
      real(dp) :: ranNum, sumInt

      ProbTable = 0d0
!      iIndx = molArray(nType)%mol(1)%indx
      do iMol = 1, NPART(nType)
         iIndx = molArray(nType)%mol(iMol)%indx
         if (NeighborList(iIndx, nTarget)) then
!          iType = typeList(i)
            ProbTable(iIndx) = exp(beta*ETable(iIndx))
!          ProbTable(i) = exp(beta * ETable(i))
         end if
!        iIndx = iIndx + 1
      end do

      norm = sum(ProbTable)
      ranNum = norm*grnd()
      sumInt = ProbTable(1)
      nSel = 1
      do while (sumInt .lt. ranNum)
         nSel = nSel + 1
         sumInt = sumInt + ProbTable(nSel)
      end do
      ProbTarget = ProbTable(nSel)/norm
      if (present(forceID)) then 
         ProbTarget = ProbTable(forceID)/norm
      end if 
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
      real(dp) :: norm

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
end module
