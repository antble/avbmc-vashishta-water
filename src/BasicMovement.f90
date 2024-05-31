!===========================================================================================
module SimpleMCMoves_Module
contains
!===========================================================================================
   subroutine SingleAtom_Translation(E_T, acc_x, atmp_x)
      use CBMC_Variables
      use Coords
!      use E_Interface
      use DistanceCriteria
      use EnergyPointers, only: Detailed_ECalc, Shift_ECalc, Update_SubEnergies, Update_SubEnergies_Vashishta
      use EnergyCriteria
      use EnergyTables
      use Forcefield
      use IndexingFunctions
      use PairStorage, only: UpdateDistArray, useDistStore
      use Pressure_LJ_Electro, only: Shift_PressCalc_Inter
      use SimParameters
      use AcceptRates
      use IndexingFunctions
      use UmbrellaSamplingNew, only: useUmbrella, GetUmbrellaBias_Disp
      implicit none

      real(dp), intent(inout) :: E_T, acc_x, atmp_x

      logical, parameter :: useIntra(1:4) = [.true., .true., .true., .true.]
      logical :: rejMove, errRtn, updated
      integer :: nType, nMol, nIndx, nMove, nAtom, iAtom, jAtom
      real(dp) :: biasDiff, biasEnergy
      real(dp) :: grnd
      real(dp) :: dx, dy, dz
      real(dp) :: E_Diff, E_Inter, E_Intra
      type(displacement) :: disp(1:maxAtoms)
      real(dp) :: PairList(1:maxMol)
      real(dp) :: dETable(1:maxMol)

      updated = .false.
      prevMoveAccepted = .false.
      rejMove = .false.
      atmp_x = atmp_x + 1E0
      
!     Randomly Select a Particle from the cluster and obtain its molecule type and index
      nMove = floor(NTotal*grnd() + 1E0)
      call Get_MolIndex(nMove, NPart, nType, nMol)
      atmpTrans(nType) = atmpTrans(nType) + 1E0
      nIndx = MolArray(nType)%mol(nMol)%indx
      nAtom = floor(nAtoms(nType)*grnd() + 1E0)
      
!     Generate a random translational displacement
      dx = max_dist_single(nType)*(2E0*grnd() - 1E0)
      dy = max_dist_single(nType)*(2E0*grnd() - 1E0)
      dz = max_dist_single(nType)*(2E0*grnd() - 1E0)
      ! print *, dx, dy, dz
   
      ! Construct the Displacement Vectors for each atom in the molecule that was chosen.
      do iAtom = 1, nAtoms(nType)
         disp(iAtom)%molType = int(nType, atomIntType)
         disp(iAtom)%molIndx = int(nMol, atomIntType)
         disp(iAtom)%atmIndx = int(iAtom, atomIntType)

         disp(iAtom)%x_old => MolArray(nType)%mol(nMol)%x(iAtom)
         disp(iAtom)%y_old => MolArray(nType)%mol(nMol)%y(iAtom)
         disp(iAtom)%z_old => MolArray(nType)%mol(nMol)%z(iAtom)

         ! update the selected atom
         if (iAtom == nAtom) then 
            disp(iAtom)%x_new = disp(iAtom)%x_old + dx
            disp(iAtom)%y_new = disp(iAtom)%y_old + dy
            disp(iAtom)%z_new = disp(iAtom)%z_old + dz
         else 
            disp(iAtom)%x_new = disp(iAtom)%x_old
            disp(iAtom)%y_new = disp(iAtom)%y_old 
            disp(iAtom)%z_new = disp(iAtom)%z_old
         end if 
      end do
      
!     Calculate the Energy Associated with this move and ensure it conforms to
!     the cluster criteria.
      E_Inter = 0E0
      E_Intra = 0E0
      
      call Shift_ECalc(E_Inter, E_Intra, disp(1:nAtoms(nType)), PairList, dETable, useIntra, rejMove)
      if (rejMove) return

      E_Diff = E_Inter 
!     Calculate Acceptance and determine if the move is accepted or not     
      if(E_Diff .le. 0E0) then
         if(calcPressure) then
            call Shift_PressCalc_Inter(P_Diff, disp(1:1))
            pressure = pressure + P_Diff
         endif

         ! update the atom position
         do iAtom = 1, nAtoms(nType)
            disp(iAtom)%x_old = disp(iAtom)%x_new
            disp(iAtom)%y_old = disp(iAtom)%y_new
            disp(iAtom)%z_old = disp(iAtom)%z_new
         end do

         acc_x = acc_x + 1E0
         if(useDistStore) then
            call UpdateDistArray
         endif
         if( distCriteria ) then
            call NeighborUpdate_Distance(PairList, nIndx)
            
         else
            call NeighborUpdate(PairList, nIndx)
         endif  

         
         select case (trim(adjustl(ForceFieldName)))
         case('Vashishta')
            call Update_MoleculeArray(updated)
            if (updated) then
               ! reset the energy tabulation
               call Detailed_ECalc(E_T, errRtn)
               ! remove the previously added E_Inter due to Shift_ECalc
               E_Inter_T = E_Inter_T -  E_Inter
            end if 
            if (.not. updated) then 
               E_T = E_T + E_Diff 
               ETable = ETable + dETable
            end if
            call Update_SubEnergies_Vashishta
         case default
            E_T = E_T + E_Diff
            ETable = ETable + dETable
            call Update_SubEnergies
         end select      
         prevMoveAccepted = .true.
      elseif(-beta*E_Diff .gt. log(grnd())) then
         if(calcPressure) then
            call Shift_PressCalc_Inter(P_Diff, disp(1:1))
            pressure = pressure + P_Diff
         endif
         
         do iAtom = 1, nAtoms(nType)
            disp(iAtom)%x_old = disp(iAtom)%x_new
            disp(iAtom)%y_old = disp(iAtom)%y_new
            disp(iAtom)%z_old = disp(iAtom)%z_new
         end do
         
         acc_x = acc_x + 1E0
         if(useDistStore) then
            call UpdateDistArray
         endif

         if( distCriteria ) then
            call NeighborUpdate_Distance(PairList, nIndx)
            
         else
            call NeighborUpdate(PairList, nIndx)
         endif

         
         select case (trim(adjustl(ForceFieldName)))
         case('Vashishta')
            call Update_MoleculeArray(updated)
            ! reset the energy tabulation
            if (updated) then
               call Detailed_ECalc(E_T, errRtn)
               ! remove the previously added E_Inter due to Shift_ECalc
               E_Inter_T = E_Inter_T -  E_Inter
            end if
            if (.not. updated) then
               E_T = E_T + E_Diff
               ETable = ETable + dETable
            end if 
            call Update_SubEnergies_Vashishta
         case default
            E_T = E_T + E_Diff
            ETable = ETable + dETable
            call Update_SubEnergies
         end select 
         prevMoveAccepted = .true.
      end if
      ! print *, "Boltzterm=",exp(-beta*E_Diff), E_Inter, beta
      ! print *, "E_Inter/outputEConv ", E_Inter/outputEConv 

   end subroutine
!===========================================================================================
   subroutine Translation(E_T, acc_x, atmp_x)
      use AcceptRates
      use SimParameters
      use IndexingFunctions
      use Coords
      use Forcefield
!      use E_Interface
      use EnergyPointers, only: Shift_ECalc, Update_SubEnergies, Update_SubEnergies_Vashishta
      use EnergyCriteria
      use DistanceCriteria
      use EnergyTables
      use PairStorage, only: UpdateDistArray
      use UmbrellaSamplingNew, only: useUmbrella, GetUmbrellaBias_Disp
      implicit none

      real(dp), intent(inout) :: E_T, acc_x, atmp_x
      logical, parameter :: useIntra(1:4) = [.false., .false., .false., .false.]
      character(len=100) :: format_string, fl_name, out1
      logical :: rejMove
      integer :: iAtom, nType, nMol, iMol, nIndx, nMove, selected_atom
      real(dp) :: biasDiff, biasEnergy
      real(dp) :: grnd
      real(dp) :: dx, dy, dz
      real(dp) :: E_Inter, E_Intra
      type(displacement) :: disp(1:maxAtoms)
      real(dp) :: PairList(1:maxMol)
      real(dp) :: dETable(1:maxMol)

      format_string = "(A,I5,A)"
      prevMoveAccepted = .false.

      if (NTotal .eq. 1) return
      
      rejMove = .false.
!     Randomly Select a Particle from the cluster and obtain its molecule type and index
      nMove = floor(NTotal*grnd() + 1E0)
      call Get_MolIndex(nMove, NPart, nType, nMol)

      atmp_x = atmp_x + 1E0
      atmpTrans(nType) = atmpTrans(nType) + 1E0
      nIndx = MolArray(nType)%mol(nMol)%indx
!     Generate a random translational displacement
      dx = max_dist(nType)*(2E0*grnd() - 1E0)
      dy = max_dist(nType)*(2E0*grnd() - 1E0)
      dz = max_dist(nType)*(2E0*grnd() - 1E0)

!     Construct the Displacement Vectors for each atom in the molecule that was chosen.
      do iAtom = 1, nAtoms(nType)
         disp(iAtom)%molType = int(nType, atomIntType)
         disp(iAtom)%molIndx = int(nMol, atomIntType)
         disp(iAtom)%atmIndx = int(iAtom, atomIntType)

         disp(iAtom)%x_old => MolArray(nType)%mol(nMol)%x(iAtom)
         disp(iAtom)%y_old => MolArray(nType)%mol(nMol)%y(iAtom)
         disp(iAtom)%z_old => MolArray(nType)%mol(nMol)%z(iAtom)

         ! update then selected atom
         disp(iAtom)%x_new = disp(iAtom)%x_old + dx
         disp(iAtom)%y_new = disp(iAtom)%y_old + dy
         disp(iAtom)%z_new = disp(iAtom)%z_old + dz 
      end do
      
      
!     Calculate the Energy Associated with this move and ensure it conforms to
!     the cluster criteria.
      E_Inter = 0E0
      E_Intra = 0E0
      !call Shift_EnergyCalc(E_Inter, E_Intra, disp(1:nAtoms(nType)), PairList, dETable, useIntra, rejMove)
      call Shift_ECalc(E_Inter, E_Intra, disp(1:nAtoms(nType)), PairList, dETable, useIntra, rejMove)
      if (rejMove) then
         return
      end if

      biasDiff = 0E0
!      write(*,*) useUmbrella
      if (useUmbrella) then
         call GetUmbrellaBias_Disp(disp(1:nAtoms(nType)), biasDiff, rejMove)
         if (rejMove) then
            return
         end if
      end if
      biasEnergy = beta*E_Inter - biasDiff
!     Calculate Acceptance and determine if the move is accepted or not
      if (biasEnergy .le. 0E0) then
         acptTrans(nType) = acptTrans(nType) + 1E0
         call Update_Shift(disp(1:nAtoms(nType)), nType, nIndx, E_T, E_Inter, acc_x, atmp_x, PairList, dETable)
!      elseif(exp(-biasEnergy) .gt. grnd()) then
      elseif (-biasEnergy .gt. log(grnd())) then
         acptTrans(nType) = acptTrans(nType) + 1E0
         call Update_Shift(disp(1:nAtoms(nType)), nType, nIndx, E_T, E_Inter, acc_x, atmp_x, PairList, dETable)
      end if
   end subroutine
!===========================================================================================
   subroutine Rotation(E_T, acc_x, atmp_x)
      use SimParameters
      implicit none
      real(dp), intent(inout) :: atmp_x, acc_x, E_T
      real(dp) :: ran_num, grnd

      if (NTotal .eq. 1) then
!       acc_x=acc_x+1E0
         return
      end if

      prevMoveAccepted = .false.

      ran_num = grnd()
      if (ran_num .lt. 1E0/3E0) then
         call Rot_xy(E_T, acc_x, atmp_x)
      elseif (ran_num .lt. 2E0/3E0) then
         call Rot_xz(E_T, acc_x, atmp_x)
      else
         call Rot_yz(E_T, acc_x, atmp_x)
      end if

   end subroutine

!=======================================================
   subroutine Rot_xy(E_T, acc_x, atmp_x)
      use AcceptRates
      use SimParameters
      use IndexingFunctions
      use Coords
      use Forcefield
!      use E_Interface
      use EnergyPointers, only: Shift_ECalc, Update_SubEnergies, Update_SubEnergies_Vashishta
      use EnergyCriteria
      use DistanceCriteria
      use EnergyTables
!      use PairStorage, only: UpdateDistArray
      use UmbrellaSamplingNew, only: useUmbrella, GetUmbrellaBias_Disp
      implicit none

      real(dp), intent(inout) :: E_T, acc_x, atmp_x
      logical, parameter :: useIntra(1:4) = [.false., .false., .false., .false.]

      logical :: rejMove
      integer :: iAtom, nMove
      integer :: atmType, nMol, nType, nIndx
      real(dp) :: E_Inter, E_Intra
      real(dp) :: grnd
      real(dp) :: biasDiff, biasEnergy
      real(dp) :: c_term, s_term
      real(dp) :: x_scale, y_scale
      real(dp) :: xcm, ycm, angle
      type(displacement) :: disp(1:maxAtoms)
      real(dp) :: PairList(1:maxMol)
      real(dp) :: dETable(1:maxMol)

!     Randomly Select a Particle from the cluster and obtain its molecule type and index
      nMove = floor(NTotal*grnd() + 1E0)
      call Get_MolIndex(nMove, NPart, nType, nMol)
      if (nAtoms(nType) .eq. 1) then
         return
      end if
      atmp_x = atmp_x + 1E0
      atmpRot(nType) = atmpRot(nType) + 1E0
      nIndx = MolArray(nType)%mol(nMol)%indx

!     Set the Displacement Array
      do iAtom = 1, nAtoms(nType)
         disp(iAtom)%molType = int(nType, atomIntType)
         disp(iAtom)%molIndx = int(nMol, atomIntType)
         disp(iAtom)%atmIndx = int(iAtom, atomIntType)

         disp(iAtom)%x_old => MolArray(nType)%mol(nMol)%x(iAtom)
         disp(iAtom)%y_old => MolArray(nType)%mol(nMol)%y(iAtom)
         disp(iAtom)%z_old => MolArray(nType)%mol(nMol)%z(iAtom)
      end do

!     Uniformly choose a random rotational displacement ranging from -max_rot to +max_rot
      angle = max_rot(nType)*(2E0*grnd() - 1E0)
      c_term = cos(angle)
      s_term = sin(angle)

!     Determine the center of mass which will act as the pivot point for the rotational motion.
      xcm = 0E0
      ycm = 0E0
      do iAtom = 1, nAtoms(nType)
         atmType = atomArray(nType, iAtom)
         xcm = xcm + atomData(atmType)%mass*disp(iAtom)%x_old
         ycm = ycm + atomData(atmType)%mass*disp(iAtom)%y_old
      end do

      xcm = xcm/totalMass(nType)
      ycm = ycm/totalMass(nType)

!     Generate a random translational displacement 
      do iAtom = 1, nAtoms(nType)
         disp(iAtom)%z_new = disp(iAtom)%z_old
         x_scale = disp(iAtom)%x_old - xcm
         y_scale = disp(iAtom)%y_old - ycm
         disp(iAtom)%x_new = c_term*x_scale - s_term*y_scale + xcm
         disp(iAtom)%y_new = s_term*x_scale + c_term*y_scale + ycm
      end do

!     Calculate the Energy Difference Associated with the move
      E_Inter = 0E0
      E_Intra = 0E0
      rejMove = .false.
!      call Shift_EnergyCalc(E_Inter, E_Intra, disp(1:nAtoms(nType)), PairList, dETable, useIntra, rejMove)
      call Shift_ECalc(E_Inter, E_Intra, disp(1:nAtoms(nType)), PairList, dETable, useIntra, rejMove)
      if (rejMove) return

      biasDiff = 0E0
      if (useUmbrella) then
         call GetUmbrellaBias_Disp(disp(1:nAtoms(nType)), biasDiff, rejMove)
         if (rejMove) then
            return
         end if
      end if
      biasEnergy = beta*E_Inter - biasDiff

!      Calculate Acceptance and determine if the move is accepted or not
      if (biasEnergy .le. 0E0) then
         acptRot(nType) = acptRot(nType) + 1E0
         call Update_Shift(disp(1:nAtoms(nType)), nType, nIndx, E_T, E_Inter, acc_x, atmp_x, PairList, dETable)
!      elseif(exp(-biasEnergy) .gt. grnd()) then
      elseif (-biasEnergy .gt. log(grnd())) then
         acptRot(nType) = acptRot(nType) + 1E0
         call Update_Shift(disp(1:nAtoms(nType)), nType, nIndx, E_T, E_Inter, acc_x, atmp_x, PairList, dETable)
      end if
   end subroutine
!=======================================================
   subroutine Rot_xz(E_T, acc_x, atmp_x)
      use AcceptRates
      use SimParameters
      use IndexingFunctions
      use Coords
      use Forcefield
!      use E_Interface
      use EnergyPointers, only: Shift_ECalc, Update_SubEnergies, Update_SubEnergies_Vashishta
      use EnergyCriteria
      use DistanceCriteria
      use EnergyTables
!      use PairStorage, only: UpdateDistArray
      use UmbrellaSamplingNew, only: useUmbrella, GetUmbrellaBias_Disp
      implicit none

      real(dp), intent(inout) :: E_T, acc_x, atmp_x

      logical, parameter :: useIntra(1:4) = [.false., .false., .false., .false.]

      logical :: rejMove
      integer :: iAtom, nMove
      integer :: atmType, nMol, nType, nIndx
      real(dp) :: angle
      real(dp) :: E_Inter, E_Intra
      real(dp) :: grnd
      real(dp) :: biasDiff, biasEnergy
      real(dp) :: c_term, s_term
      real(dp) :: x_scale, z_scale
      real(dp) :: xcm, zcm
      type(displacement) :: disp(1:maxAtoms)
      real(dp) :: PairList(1:maxMol)
      real(dp) :: dETable(1:maxMol)

      prevMoveAccepted = .false.

!     Randomly Select a Particle from the cluster and obtain its molecule type and index
      nMove = floor(NTotal*grnd() + 1E0)
      call Get_MolIndex(nMove, NPart, nType, nMol)
      if (nAtoms(nType) .eq. 1) then
         return
      end if
      atmp_x = atmp_x + 1E0
      atmpRot(nType) = atmpRot(nType) + 1E0
      nIndx = MolArray(nType)%mol(nMol)%indx

!     Set the Displacement Array
      do iAtom = 1, nAtoms(nType)
         disp(iAtom)%molType = int(nType, atomIntType)
         disp(iAtom)%molIndx = int(nMol, atomIntType)
         disp(iAtom)%atmIndx = int(iAtom, atomIntType)

         disp(iAtom)%x_old => MolArray(nType)%mol(nMol)%x(iAtom)
         disp(iAtom)%y_old => MolArray(nType)%mol(nMol)%y(iAtom)
         disp(iAtom)%z_old => MolArray(nType)%mol(nMol)%z(iAtom)
      end do

!     Uniformly choose a random rotational displacement ranging from -max_rot to +max_rot
      angle = max_rot(nType)*(2E0*grnd() - 1E0)
      c_term = cos(angle)
      s_term = sin(angle)

!     Determine the center of mass which will act as the pivot point for the rotational motion.
      xcm = 0E0
      zcm = 0E0
      do iAtom = 1, nAtoms(nType)
         atmType = atomArray(nType, iAtom)
         xcm = xcm + atomData(atmType)%mass*disp(iAtom)%x_old
         zcm = zcm + atomData(atmType)%mass*disp(iAtom)%z_old
      end do

      xcm = xcm/totalMass(nType)
      zcm = zcm/totalMass(nType)

!     Generate a random translational displacement
      do iAtom = 1, nAtoms(nType)
         disp(iAtom)%y_new = disp(iAtom)%y_old
         x_scale = disp(iAtom)%x_old - xcm
         z_scale = disp(iAtom)%z_old - zcm
         disp(iAtom)%x_new = c_term*x_scale - s_term*z_scale + xcm
         disp(iAtom)%z_new = s_term*x_scale + c_term*z_scale + zcm
      end do

!     Calculate the Energy Difference Associated with the move
      E_Inter = 0E0
      E_Intra = 0E0
      rejMove = .false.
!      call Shift_EnergyCalc(E_Inter, E_Intra, disp(1:nAtoms(nType)), PairList, dETable, useIntra, rejMove)
      call Shift_ECalc(E_Inter, E_Intra, disp(1:nAtoms(nType)), PairList, dETable, useIntra, rejMove)
      if (rejMove) return

      biasDiff = 0E0
      if (useUmbrella) then
         call GetUmbrellaBias_Disp(disp(1:nAtoms(nType)), biasDiff, rejMove)
         if (rejMove) then
            return
         end if
      end if
      biasEnergy = beta*E_Inter - biasDiff

!      Calculate Acceptance and determine if the move is accepted or not
      if (biasEnergy .le. 0E0) then
         acptRot(nType) = acptRot(nType) + 1E0
         call Update_Shift(disp(1:nAtoms(nType)), nType, nIndx, E_T, E_Inter, acc_x, atmp_x, PairList, dETable)
!      elseif(exp(-biasEnergy) .gt. grnd()) then
      elseif (-biasEnergy .gt. log(grnd())) then
         acptRot(nType) = acptRot(nType) + 1E0
         call Update_Shift(disp(1:nAtoms(nType)), nType, nIndx, E_T, E_Inter, acc_x, atmp_x, PairList, dETable)
      end if

   end subroutine
!=======================================================
   subroutine Rot_yz(E_T, acc_x, atmp_x)
      use AcceptRates
      use SimParameters
      use IndexingFunctions
      use Coords
      use Forcefield
!      use E_Interface
      use EnergyPointers, only: Shift_ECalc, Update_SubEnergies, Update_SubEnergies_Vashishta
      use EnergyCriteria
      use DistanceCriteria
      use EnergyTables
!      use PairStorage, only: UpdateDistArray
      use UmbrellaSamplingNew, only: useUmbrella, GetUmbrellaBias_Disp
      implicit none

      real(dp), intent(inout) :: E_T, acc_x, atmp_x
      logical, parameter :: useIntra(1:4) = [.false., .false., .false., .false.]

      logical :: rejMove
      integer :: iAtom, nMove
      integer :: atmType, nMol, nType, nIndx
      real(dp) :: angle
      real(dp) :: E_Inter, E_Intra
      real(dp) :: grnd
      real(dp) :: biasDiff, biasEnergy
      real(dp) :: c_term, s_term
      real(dp) :: y_scale, z_scale
      real(dp) :: ycm, zcm
      type(displacement) :: disp(1:maxAtoms)
      real(dp) :: PairList(1:maxMol)
      real(dp) :: dETable(1:maxMol)

!     Randomly Select a Particle from the cluster and obtain its molecule type and index
      nMove = floor(NTotal*grnd() + 1E0)
      call Get_MolIndex(nMove, NPart, nType, nMol)
      if (nAtoms(nType) .eq. 1) then
         return
      end if
      atmp_x = atmp_x + 1E0
      atmpRot(nType) = atmpRot(nType) + 1E0
      nIndx = MolArray(nType)%mol(nMol)%indx

!     Set the Displacement Array
      do iAtom = 1, nAtoms(nType)
         disp(iAtom)%molType = int(nType, atomIntType)
         disp(iAtom)%molIndx = int(nMol, atomIntType)
         disp(iAtom)%atmIndx = int(iAtom, atomIntType)

         disp(iAtom)%x_old => MolArray(nType)%mol(nMol)%x(iAtom)
         disp(iAtom)%y_old => MolArray(nType)%mol(nMol)%y(iAtom)
         disp(iAtom)%z_old => MolArray(nType)%mol(nMol)%z(iAtom)
      end do

!     Uniformly choose a random rotational displacement ranging from -max_rot to +max_rot
      angle = max_rot(nType)*(2E0*grnd() - 1E0)
      c_term = cos(angle)
      s_term = sin(angle)

!     Determine the center of mass which will act as the pivot point for the rotational motion.
      ycm = 0E0
      zcm = 0E0
      do iAtom = 1, nAtoms(nType)
         atmType = atomArray(nType, iAtom)
         ycm = ycm + atomData(atmType)%mass*disp(iAtom)%y_old
         zcm = zcm + atomData(atmType)%mass*disp(iAtom)%z_old
      end do

      ycm = ycm/totalMass(nType)
      zcm = zcm/totalMass(nType)

!     Generate a random translational displacement
      do iAtom = 1, nAtoms(nType)
         disp(iAtom)%x_new = disp(iAtom)%x_old
         y_scale = disp(iAtom)%y_old - ycm
         z_scale = disp(iAtom)%z_old - zcm
         disp(iAtom)%y_new = c_term*y_scale - s_term*z_scale + ycm
         disp(iAtom)%z_new = s_term*y_scale + c_term*z_scale + zcm
      end do

!     Calculate the Energy Difference Associated with the move
      E_Inter = 0E0
      E_Intra = 0E0
      rejMove = .false.
!      call Shift_EnergyCalc(E_Inter, E_Intra, disp(1:nAtoms(nType)), PairList, dETable, useIntra, rejMove)
      call Shift_ECalc(E_Inter, E_Intra, disp(1:nAtoms(nType)), PairList, dETable, useIntra, rejMove)
      if (rejMove) return

      biasDiff = 0E0
      if (useUmbrella) then
         call GetUmbrellaBias_Disp(disp(1:nAtoms(nType)), biasDiff, rejMove)
         if (rejMove) then
            return
         end if
      end if
      biasEnergy = beta*E_Inter - biasDiff
!      Calculate Acceptance and determine if the move is accepted or not
      if (biasEnergy .le. 0E0) then
         acptRot(nType) = acptRot(nType) + 1E0
         call Update_Shift(disp(1:nAtoms(nType)), nType, nIndx, E_T, E_Inter, acc_x, atmp_x, PairList, dETable)
!      elseif(exp(-biasEnergy) .gt. grnd()) then
      elseif (-biasEnergy .gt. log(grnd())) then
         acptRot(nType) = acptRot(nType) + 1E0
         call Update_Shift(disp(1:nAtoms(nType)), nType, nIndx, E_T, E_Inter, acc_x, atmp_x, PairList, dETable)
      end if
   end subroutine

!=======================================================
   subroutine Update_Shift(disp, nType, nIndx, E_T, E_Inter, acc_x, atmp_x, PairList, dETable)
      use AcceptRates
      use SimParameters
      use Coords
      use Forcefield, only: nAtoms, ForceFieldName
      use EnergyPointers, only: Shift_ECalc, Update_SubEnergies, Update_SubEnergies_Vashishta
      use EnergyCriteria
      use DistanceCriteria
      use EnergyTables
      use PairStorage, only: UpdateDistArray, useDistStore
      use UmbrellaSamplingNew, only: useUmbrella, GetUmbrellaBias_Disp
      use Pressure_LJ_Electro, only: Shift_PressCalc_Inter
      implicit none

      real(dp), intent(inout) :: E_T, acc_x, atmp_x
      integer, intent(in) :: nIndx, nType
      type(displacement), intent(inout) :: disp(:)
      real(dp), intent(in) :: PairList(1:maxMol)
      real(dp), intent(in) :: dETable(1:maxMol)
      real(dp), intent(in) :: E_Inter

      integer :: iAtom

      if (calcPressure) then
         call Shift_PressCalc_Inter(P_Diff, disp(1:nAtoms(nType)))
         pressure = pressure + P_Diff
      end if

      do iAtom = 1, nAtoms(nType)
         disp(iAtom)%x_old = disp(iAtom)%x_new
         disp(iAtom)%y_old = disp(iAtom)%y_new
         disp(iAtom)%z_old = disp(iAtom)%z_new
      end do
      E_T = E_T + E_Inter
      ETable = ETable + dETable
      acc_x = acc_x + 1E0
      if (useDistStore) then
         call UpdateDistArray
      end if
      if (distCriteria) then
!        call NeighborUpdate_Distance(PairList,nIndx)
         call NeighborUpdate_Distance(PairList, nIndx)
      else
         call NeighborUpdate(PairList, nIndx)
      end if
      
      select case (trim(adjustl(ForceFieldName)))
      case('Vashishta')
         call Update_SubEnergies_Vashishta
         !print *, "Vashishta sub energies update in Update_Shift -- ", (trim(adjustl(ForceFieldName)))
      case default
         call Update_SubEnergies
         !print *, "Default sub energies update in Update_Shift -- ", (trim(adjustl(ForceFieldName)))
      end select
      ! call Update_SubEnergies
      prevMoveAccepted = .true.

   end subroutine

!=======================================================
!     Experimental Temperature Change Move.  Not guarenteed to give accurate results.
   subroutine TemperatureMove(E_T, acc_x, atmp_x)
      use AcceptRates
      use SimParameters
      use UmbrellaSamplingNew, only: useUmbrella, GetUmbrellaBias_Temperature
      implicit none

      real(dp), intent(inout) :: E_T, acc_x, atmp_x
      real(dp), parameter :: power = (6d0/2d0)

      logical :: rejMove
      integer :: iAtom, nMove
      integer :: atmType, nMol, nType, nIndx
      real(dp) :: grnd, betaNew
      real(dp) :: biasNew, biasOld, biasDiff, biasEnergy

      if (NTotal .ne. 1) then
         return
      end if

!     Randomly Select a Particle from the cluster and obtain its molecule type and index
      atmp_x = atmp_x + 1E0_dp
      TempNew = temperature + 15E0_dp*(2E0_dp*grnd() - 1E0_dp)
      betaNew = 1E0_dp/TempNew
      biasOld = 0E0_dp
      biasNew = 0E0_dp
      if (useUmbrella) then
         call GetUmbrellaBias_Temperature(biasDiff, rejMove)
         if (rejMove) then
            return
         end if
      end if
      biasEnergy = (betaNew - beta)*E_T - biasDiff

!      Calculate Acceptance and determine if the move is accepted or not
!      if(biasEnergy .le. 0E0_dp) then
!        acc_x = acc_x + 1E0_dp
!        temperature = TempNew
!        beta = betaNew
      if ((TempNew/temperature)**power*exp(-biasEnergy) .gt. grnd()) then
         acc_x = acc_x + 1E0_dp
         temperature = TempNew
         beta = betaNew
      end if
!      write(*,*) temperature

   end subroutine
!===========================================================================================
end module
