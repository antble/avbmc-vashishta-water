!=================================================================================
module DistanceCriteria

   logical, allocatable :: ClusterMember(:), flipped(:)
contains
!=================================================================================
!     Extensive Cluster Criteria Check.  Used at the start and end of the simulation.
!     This ensures that all particles in the cluster are properly connected to each other.
!     This function also calculates the initial Neighborlist that is used throughout the simulation.
   subroutine Detailed_DistanceCriteria(PairList, rejMove)
      use Coords
      use IndexingFunctions
      use ParallelVar
      use SimParameters
      implicit none
      logical, intent(out) :: rejMove
      real(dp), intent(inout) :: PairList(:, :)

!      logical :: ClusterMember(1:maxMol)
      integer :: h, cnt
      integer :: iType, jType, iMol, jMol, iIndx, jIndx
      integer :: globIndx1, globIndx2

      if (.not. allocated(ClusterMember)) then
         allocate (ClusterMember(1:maxMol))
      end if
      if (.not. allocated(flipped)) then
         allocate (flipped(1:maxMol))
      end if

      rejMove = .false.
      NeighborList = .false.
      if (NTotal .eq. 1) then
         return
      end if

      do iType = 1, nMolTypes
         do iMol = 1, NPART(iType)
            globIndx1 = MolArray(iType)%mol(iMol)%globalIndx(1)
            iIndx = MolArray(iType)%mol(iMol)%indx
            do jType = 1, nMolTypes
               do jMol = 1, NPART(jType)
                  globIndx2 = MolArray(jType)%mol(jMol)%globalIndx(1)
                  jIndx = MolArray(jType)%mol(jMol)%indx
                  if (jIndx .eq. iIndx) then
                     cycle
                  end if
                  if (PairList(iIndx, jIndx) .lt. Dist_Critr_sq) then
                     NeighborList(iIndx, jIndx) = .true.
                     NeighborList(jIndx, iIndx) = .true.
                  else
                     NeighborList(iIndx, jIndx) = .false.
                     NeighborList(jIndx, iIndx) = .false.
                  end if
               end do
            end do
         end do
      end do

      if (all(NeighborList .eqv. .false.)) then
         write (nout, *) "------- Cluster Criteria Not Met!!!! --------"
         rejMove = .true.
         return
      end if

      do iType = 1, nMolTypes
         if (NPART(iType) .gt. 0) then
            iIndx = MolArray(iType)%mol(1)%indx
            ClusterMember(iIndx) = .true.
            exit
         end if
      end do

      do h = 1, maxMol
         do iIndx = 1, maxMol
            do jIndx = 1, maxMol
               if (NeighborList(iIndx, jIndx)) then
                  if (ClusterMember(iIndx)) then
                     ClusterMember(jIndx) = .true.
!                cnt = cnt + 1
                  end if
                  if (ClusterMember(jIndx)) then
                     ClusterMember(iIndx) = .true.
!                cnt = cnt + 1
                  end if
               end if
            end do
         end do
      end do

      cnt = 0
      do iType = 1, nMolTypes
         if (NPART(iType) .lt. NMAX(iType)) then
            do iMol = NPART(iType) + 1, NMAX(iType)
               iIndx = MolArray(iType)%mol(iMol)%indx
               ClusterMember(iIndx) = .true.
               cnt = cnt + 1
            end do
         end if
      end do

      if (any(ClusterMember .eqv. .false.)) then
         rejMove = .true.
         write (nout, *) "------- Cluster Criteria Not Met!!!! --------"
      end if

   end subroutine
!=================================================================================
!    This subroutine update the molecular array of the system after a translational move
!    This function also calculates the initial Neighborlist that is used throughout the simulation.
   ! subroutine Update_MoleculeArray(PairList, nIndx)
   subroutine Update_MoleculeArray(updated) !dETable
      use Coords
      use ForceField
      use ForceFieldPara_Vashishta
      use IndexingFunctions
      use SimParameters
      use EnergyTables
      use DistanceCriteria_PairStore, only: NeighborUpdate_Distance_PairStore
      use PairStorage, only: useDistStore
      implicit none

      logical, intent(out) :: updated
      ! integer, intent(in) :: nIndx
      ! real(dp), intent(inout) :: dETable(:)
      ! real(dp), intent(in) :: PairList(:)
      integer :: nType, nMol, jMol, jType, jIndx, shortest_indx, shortest_jAtom
      integer :: num, num_update
      real(dp) :: O_x1, O_y1, O_z1, rx1, ry1, rz1, rx2, ry2, rz2, r1_ij, r2_ij, r_ij
      integer :: iType, iMol, iIndx, iAtom, jAtom, ctr 
      real(dp) ::  rx, ry, rz, tmp_x1, tmp_y1, tmp_z1, tmp_x2, tmp_y2, tmp_z2, shortest_dist
      real(dp) :: r_ij1, r_ij2, O_x2, O_y2, O_z2

      ctr = 0
      iType = 1
      jType = 1
      
      do iMol = 1, NPART(jType)
         iIndx = MolArray(jType)%mol(iMol)%indx
         ! oxygen
         O_x1 = MolArray(iType)%mol(iIndx)%x(1)
         O_y1 = MolArray(iType)%mol(iIndx)%y(1)
         O_z1 = MolArray(iType)%mol(iIndx)%z(1)
         ! first hydrogen 
         tmp_x1 = MolArray(iType)%mol(iIndx)%x(2)
         tmp_y1 = MolArray(iType)%mol(iIndx)%y(2)
         tmp_z1 = MolArray(iType)%mol(iIndx)%z(2)
         rx1 = tmp_x1 - O_x1
         ry1 = tmp_y1 - O_y1
         rz1 = tmp_z1 - O_z1
         r1_ij = sqrt(rx1**2 + ry1**2 + rz1**2)  
         ! Second hydrogen
         tmp_x2 = MolArray(iType)%mol(iIndx)%x(3)
         tmp_y2 = MolArray(iType)%mol(iIndx)%y(3)
         tmp_z2 = MolArray(iType)%mol(iIndx)%z(3)
         rx2 = tmp_x2 - O_x1
         ry2 = tmp_y2 - O_y1
         rz2 = tmp_z2 - O_z1
         r2_ij = sqrt(rx2**2 + ry2**2 + rz2**2) 


         ! num_update is the number of hydrogen to be updated 
         num_update = 2
         ! if ((r1_ij > 1.4) .and. (r2_ij > 1.4)) then
         !    num_update = 2 
         ! else if ((r1_ij < 1.4 .and. r2_ij < 1.4)) then 
         !    cycle 
         ! end if 

         ! check the neighboring molecule's hydrogen
         do num = 1, num_update 
            ! first hydrogen 
            tmp_x1 = MolArray(iType)%mol(iIndx)%x(2)
            tmp_y1 = MolArray(iType)%mol(iIndx)%y(2)
            tmp_z1 = MolArray(iType)%mol(iIndx)%z(2)
            rx1 = tmp_x1 - O_x1
            ry1 = tmp_y1 - O_y1
            rz1 = tmp_z1 - O_z1
            r1_ij = sqrt(rx1**2 + ry1**2 + rz1**2)  
            ! Second hydrogen
            tmp_x2 = MolArray(iType)%mol(iIndx)%x(3)
            tmp_y2 = MolArray(iType)%mol(iIndx)%y(3)
            tmp_z2 = MolArray(iType)%mol(iIndx)%z(3)
            rx2 = tmp_x2 - O_x1
            ry2 = tmp_y2 - O_y1
            rz2 = tmp_z2 - O_z1
            r2_ij = sqrt(rx2**2 + ry2**2 + rz2**2) 

            shortest_indx = iIndx
            shortest_dist = 99
            shortest_jAtom = 2
            ! LOOP OVER ALL THE HYDROGEN THAT ARE NEIGHBORS OF THE TARGET OXYGEN
            do jMol = 1, NPART(jType)
               if (jMol == iIndx) then 
                  cycle 
               end if 
               O_x2 = MolArray(iType)%mol(jMol)%x(1)
               O_y2 = MolArray(iType)%mol(jMol)%y(1)
               O_z2 = MolArray(iType)%mol(jMol)%z(1)

               jIndx = MolArray(jType)%mol(jMol)%indx
               ! if (NeighborList(iIndx, jIndx)) then
                  do jAtom = 2, nAtoms(jType)
                     ! distance of the hydrogen to its target oxygen
                     rx = MolArray(iType)%mol(jMol)%x(jAtom) - O_x1
                     ry = MolArray(iType)%mol(jMol)%y(jAtom) - O_y1
                     rz = MolArray(iType)%mol(jMol)%z(jAtom) - O_z1 
                     r_ij1 = sqrt(rx**2 + ry**2 + rz**2) 
                     ! distance of the hydrogen to its mother oxygen
                     rx = MolArray(iType)%mol(jMol)%x(jAtom) - O_x2
                     ry = MolArray(iType)%mol(jMol)%y(jAtom) - O_y2
                     rz = MolArray(iType)%mol(jMol)%z(jAtom) - O_z2 
                     r_ij2 = sqrt(rx**2 + ry**2 + rz**2) 

                     ! swap position only if the hydrogen is closer to the target oxygen
                     ! if (r_ij1 > 1.4) then 
                     !    cycle
                     ! end if 
                     if (r_ij1 < r_ij2) then! .and. (r_ij1 <= 1.4 .or. r_ij2 <= 1.4)) then
                        if (r_ij1 <= shortest_dist) then
                           shortest_indx = jMol
                           shortest_dist = r_ij1
                           shortest_jAtom = jAtom
                        end if
                     else 
                        cycle
                     end if 
                  end do
            end do

            ! if ((shortest_dist <= 1.4)) then ! .and. (r1_ij > 1.4 .or. r2_ij > 1.4)) then
            if (shortest_dist .ne. 99) then
               ! replace the first hydrogen 
               jAtom = shortest_jAtom
               if ((shortest_dist .le. r1_ij) .and. (shortest_dist .ge. r2_ij)) then
                  ctr = ctr + 1
                  ! replace the first hydrogen with the 
                  MolArray(iType)%mol(iIndx)%x(2) = MolArray(iType)%mol(shortest_indx)%x(jAtom)
                  MolArray(iType)%mol(iIndx)%y(2) = MolArray(iType)%mol(shortest_indx)%y(jAtom)
                  MolArray(iType)%mol(iIndx)%z(2) = MolArray(iType)%mol(shortest_indx)%z(jAtom)
                  MolArray(iType)%mol(shortest_indx)%x(jAtom) = tmp_x1
                  MolArray(iType)%mol(shortest_indx)%y(jAtom) = tmp_y1 
                  MolArray(iType)%mol(shortest_indx)%z(jAtom) = tmp_z1 
               else if ((shortest_dist .ge. r1_ij) .and. (shortest_dist .le. r2_ij) ) then
                  ctr = ctr + 1
                  ! replace the first hydrogen with the
                  MolArray(iType)%mol(iIndx)%x(3) = MolArray(iType)%mol(shortest_indx)%x(jAtom)
                  MolArray(iType)%mol(iIndx)%y(3) = MolArray(iType)%mol(shortest_indx)%y(jAtom)
                  MolArray(iType)%mol(iIndx)%z(3) = MolArray(iType)%mol(shortest_indx)%z(jAtom)
                  MolArray(iType)%mol(shortest_indx)%x(jAtom) = tmp_x2
                  MolArray(iType)%mol(shortest_indx)%y(jAtom) = tmp_y2 
                  MolArray(iType)%mol(shortest_indx)%z(jAtom) = tmp_z2  
               else if ((shortest_dist .le. r1_ij) .and. (shortest_dist  .le. r2_ij)) then 
                  ctr = ctr + 1
                  ! compare the distance between the original hydrogen and the new hydrogen
                  rx = MolArray(iType)%mol(shortest_indx )%x(jAtom) - tmp_x2
                  ry = MolArray(iType)%mol(shortest_indx )%y(jAtom) - tmp_y2 
                  rz = MolArray(iType)%mol(shortest_indx )%z(jAtom) - tmp_z2
                  r2_ij = sqrt(rx**2 + ry**2 + rz**2) 

                  rx = MolArray(iType)%mol(shortest_indx)%x(jAtom) - tmp_x1
                  ry = MolArray(iType)%mol(shortest_indx)%y(jAtom) - tmp_y1
                  rz = MolArray(iType)%mol(shortest_indx)%z(jAtom) - tmp_z1
                  r1_ij = sqrt(rx**2 + ry**2 + rz**2) 
                  if (r1_ij < r2_ij) then
                     MolArray(iType)%mol(iIndx)%x(2) = MolArray(iType)%mol(shortest_indx)%x(jAtom)
                     MolArray(iType)%mol(iIndx)%y(2) = MolArray(iType)%mol(shortest_indx)%y(jAtom)
                     MolArray(iType)%mol(iIndx)%z(2) = MolArray(iType)%mol(shortest_indx)%z(jAtom)

                     MolArray(iType)%mol(shortest_indx)%x(jAtom) = tmp_x1
                     MolArray(iType)%mol(shortest_indx)%y(jAtom) = tmp_y1 
                     MolArray(iType)%mol(shortest_indx)%z(jAtom) = tmp_z1 
                  else
                     MolArray(iType)%mol(iIndx)%x(3) = MolArray(iType)%mol(shortest_indx)%x(jAtom)
                     MolArray(iType)%mol(iIndx)%y(3) = MolArray(iType)%mol(shortest_indx)%y(jAtom)
                     MolArray(iType)%mol(iIndx)%z(3) = MolArray(iType)%mol(shortest_indx)%z(jAtom)
                     MolArray(iType)%mol(shortest_indx)%x(jAtom) = tmp_x2
                     MolArray(iType)%mol(shortest_indx)%y(jAtom) = tmp_y2 
                     MolArray(iType)%mol(shortest_indx)%z(jAtom) = tmp_z2 
                  end if 
               end if 
            end if
         end do  
      end do 
      updated = .false.
      if (ctr > 0) then 
         updated = .true.
      end if
   end subroutine

!===================================================================================
!     This function determines if a given translational move will destroy a cluster.
   subroutine Shift_DistanceCriteria(PairList, nIndx, rejMove)
      use Coords
      use IndexingFunctions
      use SimParameters
      implicit none

      logical, intent(out) :: rejMove
      real(dp), intent(in) :: PairList(:)
      integer, intent(in) :: nIndx

      logical :: neiFlipped, memberAdded
      logical :: ClusterMember(1:maxMol)
      logical :: flipped(1:maxMol)
      integer :: iIndx, jIndx, h
      integer :: nType, nMol, i, jType, jMol, globIndx1, globIndx2
      integer :: curNeigh(1:60), neiMax

      rejMove = .false.
      if (NTotal .eq. 1) return
      ClusterMember = .false.
      flipped = .false.

      nType = typeList(nIndx)
      nMol = subIndxList(nIndx)
      globIndx1 = molArray(nType)%mol(nMol)%globalIndx(1)

!     This section dermines which molecules are neighbored with the new trial position.  In the event
!     that the molecule's new location has no neghibors all further calcualtions are skipped and the move is
!     rejected.

      memberAdded = .false.
      do jIndx = 1, maxMol
         if (.not. isActive(jIndx)) then
            cycle
         end if
         if (nIndx .eq. jIndx) then
            cycle
         end if
         if (PairList(jIndx) .lt. Dist_Critr_sq) then
            ClusterMember(jIndx) = .true.
            memberAdded = .true.
         end if
      end do

      if (.not. memberAdded) then
         rejMove = .true.
         return
      end if

!      This part of the code tabulates all the neighbors located around the particle's old position.
      neiMax = 0
      curNeigh = 0
      do iIndx = 1, maxMol
         if (NeighborList(iIndx, nIndx)) then
            if (iIndx .ne. nIndx) then
               if (isActive(iIndx)) then
                  neiMax = neiMax + 1
                  curNeigh(neiMax) = iIndx
               end if
            end if
         end if
      end do

!     This section performs a quick check to see if the molecules that were neighbored with the old position
!     are part of the new cluster.  If all the old neighbors are indeed part of the cluster then no furth
!     calculations are needed.
      neiFlipped = .true.
      do iIndx = 1, neiMax
         if (.not. clusterMember(curNeigh(iIndx))) then
            neiFlipped = .false.
            exit
         end if
      end do

      if (neiFlipped) then
         rejMove = .false.
         return
      end if

      do h = 1, NTotal
!        cnt = 0
         memberAdded = .false.
         do iIndx = 1, maxMol
            if (ClusterMember(iIndx) .neqv. flipped(iIndx)) then
               do jIndx = 1, maxMol
                  if (NeighborList(iIndx, jIndx)) then
                     if (jIndx .ne. nIndx) then
                        ClusterMember(jIndx) = .true.
                        memberAdded = .true.
                     end if
                  end if
               end do
               flipped(iIndx) = .true.
            end if
         end do

         neiFlipped = .true.
         do i = 1, neiMax
            if (.not. clusterMember(curNeigh(i))) then
               neiFlipped = .false.
               exit
            end if
         end do
         if (neiFlipped) then
            exit
         else
            if (.not. memberAdded) then
               exit
            end if
         end if
      end do

      if (.not. neiFlipped) then
         rejMove = .true.
      end if

   end subroutine

!=================================================================================
!     This function determines if removing a particle from the cluster will result in the destruction of the cluster criteria.
   subroutine SwapOut_DistanceCriteria(nSwap, rejMove)
      use SimParameters
      use Coords
      use IndexingFunctions
      implicit none

      logical, intent(out) :: rejMove
      integer, intent(inout) :: nSwap
      integer :: iIndx, jIndx, h, cnt

      rejMove = .false.
      ClusterMember = .false.
      flipped = .false.

      do iIndx = 1, maxMol
         if (.not. isActive(iIndx)) then
            ClusterMember(iIndx) = .true.
            flipped(iIndx) = .true.
         end if
      end do

!      In order to initialize the cluster criteria search, a single particle must be chosen as the starting point.
      do iIndx = 1, maxMol
         if (isActive(iIndx)) then
            if (nSwap .ne. iIndx) then
               ClusterMember(iIndx) = .true.
               exit
            end if
         end if
      end do

      cnt = 0
      do iIndx = 1, maxMol
         if (isActive(iIndx) .eqv. .false.) then
            ClusterMember(iIndx) = .true.
            flipped(iIndx) = .true.
            cnt = cnt + 1
         end if
      end do

      do h = 1, maxMol
         do iIndx = 1, maxMol
            if (iIndx .ne. nSwap) then
               if (ClusterMember(iIndx) .neqv. flipped(iIndx)) then
                  do jIndx = 1, maxMol
                     if (NeighborList(iIndx, jIndx)) then
                        ClusterMember(jIndx) = .true.
                        cnt = cnt + 1
                     end if
                  end do
                  flipped(iIndx) = .true.
               end if
            end if
         end do
         if (cnt .eq. maxMol - 1) then
            exit
         end if
      end do

      ClusterMember(nSwap) = .true.

      if (any(ClusterMember .eqv. .false.)) then
         rejMove = .true.
      end if

   end subroutine
!=================================================================================
!     This function updates the neighborlist if a move is accepted.
   subroutine NeighborUpdate_Distance(PairList, nIndx)
      use Coords
      use IndexingFunctions
      use SimParameters
      use DistanceCriteria_PairStore, only: NeighborUpdate_Distance_PairStore
      use PairStorage, only: useDistStore
      implicit none
      integer, intent(in) :: nIndx
      real(dp), intent(in) :: PairList(:)
      integer :: nType, nMol, jMol, jType, jIndx

      if (useDistStore) then
         call NeighborUpdate_Distance_PairStore(nIndx)
         return
      end if

      nType = typeList(nIndx)
      nMol = subIndxList(nIndx)
      do jIndx = 1, maxMol
         if (.not. isActive(jIndx)) then
            cycle
         end if
         if (jIndx .ne. nIndx) then
            if (PairList(jIndx) .lt. Dist_Critr_sq) then
               NeighborList(nIndx, jIndx) = .true.
               NeighborList(jIndx, nIndx) = .true.
            else
               NeighborList(nIndx, jIndx) = .false.
               NeighborList(jIndx, nIndx) = .false.
            end if
         end if
      end do
!      NeighborList(nIndx,nIndx) = .false.

   end subroutine
!=================================================================================
!     This function updates the neighborlist if a move is accepted.
   subroutine NeighborUpdate_SwapIn_Distance(PairList, nType)
      use Coords
      use IndexingFunctions
      use SimParameters
      use DistanceCriteria_PairStore, only: NeighborUpdate_SwapIn_Distance_PairStore
      use PairStorage, only: useDistStore
      implicit none
      integer, intent(in) :: nType
      real(dp), intent(in) :: PairList(:)
      integer ::  nIndx, nMol, jMol, jType, jIndx

      if (useDistStore) then
         call NeighborUpdate_SwapIn_Distance_PairStore(nType)
         return
      end if

      nMol = NPART(nType) + 1
      nIndx = MolArray(nType)%mol(nMol)%indx
!      write(35,*) nType, nMol, nIndx , globIndx1
      do jType = 1, nMolTypes
         do jMol = 1, NPART(jType)
            jIndx = MolArray(jType)%mol(jMol)%indx
            if (jIndx .ne. nIndx) then
               if (PairList(jIndx) .lt. Dist_Critr_sq) then
                  NeighborList(nIndx, jIndx) = .true.
                  NeighborList(jIndx, nIndx) = .true.
               else
                  NeighborList(nIndx, jIndx) = .false.
                  NeighborList(jIndx, nIndx) = .false.
               end if
            end if
         end do
      end do

   end subroutine
!=================================================================================
end module

