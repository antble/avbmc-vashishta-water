!=================================================================================
module DistanceCriteria

   logical, allocatable :: ClusterMember(:), flipped(:)
contains
!=================================================================================
!     Extensive Cluster Criteria Check.  Used at the start and end of the simulation.
!     This ensures that all particles in the cluster are properly connected to each other.
!     This function also calculates the initial Neighborlist that is used throughout the simulation.
   subroutine Detailed_DistanceCriteria(PairList, rejMove)
      use SimParameters
      use Coords
      use IndexingFunctions
      use ParallelVar
      implicit none
      logical, intent(out) :: rejMove
      real(dp), intent(inout) :: PairList(:, :)

!      logical :: ClusterMember(1:maxMol)
      integer :: i, j, h, cnt
      integer :: iType, jType, iMol, jMol, iIndx, jIndx

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
         do jType = iType, nMolTypes
            do iMol = 1, NPART(iType)
               iIndx = MolArray(iType)%mol(iMol)%indx
               do jMol = 1, NPART(jType)
                  jIndx = MolArray(jType)%mol(jMol)%indx
                  if (PairList(iIndx, jIndx) .le. Dist_Critr_sq) then
                     NeighborList(iIndx, jIndx) = .true.
                     NeighborList(jIndx, iIndx) = .true.
                  end if
               end do
            end do
         end do
      end do

      cnt = 0
      do i = 1, maxMol
         if (isActive(i) .eqv. .false.) then
            ClusterMember(i) = .true.
            cnt = cnt + 1
         end if
      end do

      ClusterMember(1) = .true.
      do h = 1, maxMol
         do i = 1, maxMol
            do j = 1, maxMol
               if (NeighborList(i, j)) then
                  if (ClusterMember(i)) then
                     ClusterMember(j) = .true.
!                cnt = cnt + 1
                  end if
                  if (ClusterMember(j)) then
                     ClusterMember(i) = .true.
!                cnt = cnt + 1
                  end if
               end if
            end do
         end do
      end do

      if (any(ClusterMember .eqv. .false.)) then
         rejMove = .true.
         write (nout, *) "------- Cluster Criteria Not Met! --------"
      end if

   end subroutine
!=================================================================================
!     This function determines if a given translational move will destroy a cluster.
   subroutine Shift_DistanceCriteria(PairList, nIndx, rejMove)
      use SimParameters
      use Coords
      use IndexingFunctions
      implicit none

      logical, intent(out) :: rejMove
      real(dp), intent(in) :: PairList(:)
      integer, intent(in) :: nIndx

!      logical :: ClusterMember(1:maxMol)
!      logical :: flipped(1:maxMol)
      integer :: i, j, h, cnt
      integer :: nType, jType

      rejMove = .false.
      if (NTotal .eq. 1) return
      ClusterMember = .false.
      flipped = .false.

      nType = Get_MolType(nIndx, NMAX)

!     This block performs a quick check to see if any neighbors were lost
!     in the process of making this move.  If all neighbors were retained
!     then the cluster remains in tact and the detailed calculations are no
!     longer needed.
      do j = 1, maxMol + 1
         if (j .eq. maxMol + 1) then
            return
         end if
         if (j .ne. nIndx) then
            if (NeighborList(nIndx, j)) then
!            jType = Get_MolType(j,NMAX)
               jType = typeList(j)
               if (PairList(j) .gt. Dist_Critr_sq) then
                  exit
               end if
            end if
         end if
      end do

!     This section dermines which molecules are neighbored with the new trial position
      cnt = 0
      do j = 1, maxMol
         if (j .ne. nIndx) then
!          jType = Get_MolType(j,NMAX)
            jType = typeList(j)
            if (PairList(j) .le. Dist_Critr_sq) then
               ClusterMember(j) = .true.
               cnt = cnt + 1
            end if
         end if
      end do

!     This section checked to see if there were any neighbors for the molecule's new position.
!     If cnt is equal to 0 the new position has no neighbors which implies the cluster is broken.
      if (cnt .eq. 0) then
         rejMove = .true.
         return
      end if

!
      do i = 1, maxMol
         if (isActive(i) .eqv. .false.) then
            ClusterMember(i) = .true.
            flipped(i) = .true.
            cnt = cnt + 1
         end if
      end do

      ClusterMember(nIndx) = .true.
      flipped(nIndx) = .true.
      cnt = cnt + 1

      do h = 1, NTotal
         do i = 1, maxMol
            if (isActive(i)) then
               if (ClusterMember(i) .neqv. flipped(i)) then
                  do j = 1, maxMol
                     if (.not. isActive(j)) then
                        cycle
                     end if
                     if (NeighborList(i, j)) then
                        ClusterMember(j) = .true.
                        cnt = cnt + 1
                     end if
                  end do
                  flipped(i) = .true.
               end if
            end if
         end do
         if (cnt .eq. maxMol) exit
      end do

      if (any(ClusterMember .eqv. .false.)) then
         rejMove = .true.
      end if

   end subroutine

!=================================================================================
!     This function determines if a given translational move will destroy a cluster.
   subroutine SwapOut_DistanceCriteria(nSwap, rejMove)
      use SimParameters
      use Coords
      use IndexingFunctions
      implicit none

      logical, intent(out) :: rejMove
      integer, intent(inout) :: nSwap

!      logical :: ClusterMember(1:maxMol)
!      logical :: flipped(1:maxMol)
      integer i, j, h, cnt

      rejMove = .false.
      ClusterMember = .false.
      flipped = .false.

      do i = 1, maxMol
         if (isActive(i) .eqv. .true.) then
            if (nSwap .ne. i) then
               ClusterMember(i) = .true.
               exit
            end if
         end if
      end do

      cnt = 0
      do i = 1, maxMol
         if (isActive(i) .eqv. .false.) then
            ClusterMember(i) = .true.
            flipped(i) = .true.
            cnt = cnt + 1
         end if
      end do

      do h = 1, maxMol
         do i = 1, maxMol
            if (i .ne. nSwap) then
               if (ClusterMember(i) .neqv. flipped(i)) then
                  do j = 1, maxMol
                     if (NeighborList(i, j)) then
                        ClusterMember(j) = .true.
                        cnt = cnt + 1
                     end if
                  end do
                  flipped(i) = .true.
               end if
            end if
         end do
         if (cnt .eq. maxMol - 1) exit
      end do

      ClusterMember(nSwap) = .true.

      if (any(ClusterMember .eqv. .false.)) then
         rejMove = .true.
      end if

   end subroutine
!=================================================================================
!     This function updates the neighborlist if a move is accepted.
   subroutine NeighborUpdate_Distance(PairList, nIndx)
      use SimParameters
      use IndexingFunctions
      use Coords
      implicit none
      integer, intent(in) :: nIndx
      real(dp), intent(in) :: PairList(:)
      integer iType, j, jType

!      NeighborList(nIndx,:)=.false.
!      NeighborList(:,nIndx)=.false.
      iType = Get_MolType(nIndx, NMAX)
      do j = 1, maxMol
         if (isActive(j)) then
            if (j .ne. nIndx) then
               jType = Get_MolType(j, NMAX)
               if (PairList(j) .le. Dist_Critr_sq) then
                  NeighborList(nIndx, j) = .true.
                  NeighborList(j, nIndx) = .true.
               else
                  NeighborList(nIndx, j) = .false.
                  NeighborList(j, nIndx) = .false.
               end if
            end if
         end if
      end do

   end subroutine
!=================================================================================
end module

