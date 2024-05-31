!=================================================================================
module EnergyCriteria
contains
!=================================================================================
!     Extensive Cluster Criteria Check.  Used at the start and end of the simulation.
!     This ensures that all particles in the cluster are properly connected to each other.
!     This function also calculates the initial Neighborlist that is used throughout the simulation.
   subroutine Detailed_EnergyCriteria(PairList, rejMove)
      use SimParameters
      use Coords
      use IndexingFunctions
      use ParallelVar
      implicit none
      logical, intent(out) :: rejMove
      real(dp), intent(inout) :: PairList(1:maxMol, 1:maxMol)

      logical :: ClusterMember(1:maxMol)
      integer :: i, j, h, cnt
      integer :: iType, jType, iMol, jMol, iIndx, jIndx

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
                  if (PairList(iIndx, jIndx) .le. Eng_Critr(iType, jType)) then
                     NeighborList(iIndx, jIndx) = .true.
                     NeighborList(jIndx, iIndx) = .true.
                  end if
               end do
            end do
         end do
      end do

      if (all(NeighborList .eqv. .false.)) then
         rejMove = .true.
         write (nout, *) "------- Cluster Criteria Not Met! --------"
      end if

      cnt = 0
      do i = 1, maxMol
         if (isActive(i) .eqv. .false.) then
            ClusterMember(i) = .true.
            cnt = cnt + 1
         end if
      end do

      do i = 1, maxMol
         if (isActive(i) .eqv. .true.) then
            ClusterMember(i) = .true.
            exit
         end if
      end do

!      ClusterMember(1)=.true.
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

      do i = 1, maxMol
         NeighborList(i, i) = .false.
      end do

      if (any(ClusterMember .eqv. .false.)) then
         rejMove = .true.
         write (nout, *) "------- Cluster Criteria Not Met! --------"
      end if

   end subroutine
!=================================================================================
!     This function determines if a given translational move will destroy a cluster.
   subroutine Shift_EnergyCriteria(PairList, nIndx, rejMove)
      use SimParameters
      use Coords
      use IndexingFunctions
      implicit none

      logical, intent(out) :: rejMove
      real(dp), intent(in) :: PairList(1:maxMol)
      integer, intent(in) :: nIndx

      logical :: neiFlipped, memberAdded
      logical :: ClusterMember(1:maxMol)
      logical :: flipped(1:maxMol)
      integer :: i, j, h
      integer :: nType, jType
      integer :: curNeigh(1:60), neiMax

      rejMove = .false.
      if (NTotal .eq. 1) return
      ClusterMember = .false.
      flipped = .false.

      nType = Get_MolType(nIndx, NMAX)

!     This section dermines which molecules are neighbored with the new trial position.  In the event
!     that the molecule's new location has no neghibors all further calcualtions are skipped and the move is
!     rejected.

      memberAdded = .false.
      do j = 1, maxMol
         jType = typeList(j)
         if (PairList(j) .le. Eng_Critr(jType, nType)) then
            ClusterMember(j) = .true.
            memberAdded = .true.
         end if
      end do

      if (.not. memberAdded) then
         rejMove = .true.
         return
      end if

!      This part of the code tabulates all the neighbors
      neiMax = 0
      curNeigh = 0
      do i = 1, maxMol
         if (NeighborList(i, nIndx)) then
            if (i .ne. nIndx) then
               if (isActive(i)) then
                  neiMax = neiMax + 1
                  curNeigh(neiMax) = i
               end if
            end if
         end if
      end do

!     This section performs a quick check to see if the molecules that were neighbored with the old position
!     are part of the new cluster.  If all the old neighbors are indeed part of the cluster then no furth
!     calculations are needed.
      neiFlipped = .true.
      do i = 1, neiMax
         if (.not. clusterMember(curNeigh(i))) then
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
         do i = 1, maxMol
            if (ClusterMember(i) .neqv. flipped(i)) then
               do j = 1, maxMol
                  if (NeighborList(i, j)) then
                     if (j .ne. nIndx) then
                        ClusterMember(j) = .true.
                        memberAdded = .true.
                     end if
                  end if
               end do
               flipped(i) = .true.
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
!     This function determines if a given translational move will destroy a cluster.
   pure subroutine SwapIn_EnergyCriteria(nType, PairList, rejMove)
      use SimParameters
      use Coords
      use IndexingFunctions
      implicit none

      logical, intent(out) :: rejMove
      real(dp), intent(in) :: PairList(1:maxMol)
      integer, intent(in) :: nType
      integer :: j, jType

      do j = 1, maxMol
         if (isActive(j)) then
            jType = Get_MolType(j, NMAX)
            if (PairList(j) .le. Eng_Critr(nType, jType)) then
               rejMove = .false.
               return
            end if
         end if
      end do

      rejMove = .true.

   end subroutine
!=================================================================================
!     This function determines if a given translational move will destroy a cluster.
   subroutine SwapOut_EnergyCriteria(nSwap, rejMove)
      use SimParameters
      use Coords
      use IndexingFunctions
      implicit none

      logical, intent(out) :: rejMove
      integer, intent(inout) :: nSwap

      logical :: neiFlipped, memberAdded
      logical :: clusterMember(1:maxMol)
      logical :: flipped(1:maxMol)
      integer :: i, j, h
      integer :: curNeigh(1:60), neiMax

      rejMove = .false.
!      if(NTotal-1 .eq. 1) return

      ClusterMember = .false.
      flipped = .false.
      neiFlipped = .false.

!      write(2,*) "----------------------------"
!      write(2,*) NPART(1), nSwap
      do i = 1, maxMol
         if (.not. isActive(i)) then
            cycle
         end if
         do j = 1, maxMol
            if (.not. isActive(j)) then
               cycle
            end if
!          write(2,*) i, j, NeighborList(i, nSwap), NeighborList(nSwap, i)
         end do
      end do

      do i = 1, maxMol
         if (.not. isActive(i)) then
            cycle
         end if
         if (NeighborList(i, nSwap)) then
            if (i .ne. nSwap) then
               clusterMember(i) = .true.
               exit
            end if
         end if
      end do

      neiMax = 0
      curNeigh = 0
      do i = 1, maxMol
         if (.not. isActive(i)) then
            cycle
         end if
         if (NeighborList(i, nSwap)) then
            if (i .ne. nSwap) then
               neiMax = neiMax + 1
               curNeigh(neiMax) = i
            end if
         end if
      end do

      do h = 1, NTotal
         memberAdded = .false.
         do i = 1, maxMol
            if (ClusterMember(i) .neqv. flipped(i)) then
               if (i .ne. nSwap) then
                  do j = 1, maxMol
                     if (NeighborList(i, j)) then
                        if (j .ne. nSwap) then
                           clusterMember(j) = .true.

                           memberAdded = .true.
                        end if
                     end if
                  end do
                  flipped(i) = .true.
               end if
            end if
         end do
!        do i = 1, maxMol
!          if(isActive(i)) then
!            write(2,*) i, clusterMember(i), flipped(i)
!          endif
!        enddo
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
!     This function updates the neighborlist if a move is accepted.
   subroutine NeighborUpdate(PairList, nIndx)
      use SimParameters
      use IndexingFunctions
      use Coords
      implicit none
      integer iType, j, jType, nIndx
      real(dp) :: PairList(1:maxMol)

!      do j=1,maxMol
!        if(.not. isActive(j)) cycle
!        if(j .ne. nIndx) then
!            NeighborList(nIndx,j)=.true.
!            NeighborList(j,nIndx)=.true.
!          else
!            NeighborList(nIndx,j)=.false.
!            NeighborList(j,nIndx)=.false.
!          endif
!        endif
!      enddo

      iType = Get_MolType(nIndx, NMAX)
      do j = 1, maxMol
         if (.not. isActive(j)) then
            cycle
         end if
         if (j .ne. nIndx) then
            jType = Get_MolType(j, NMAX)
            if (PairList(j) .le. Eng_Critr(iType, jType)) then
               NeighborList(nIndx, j) = .true.
               NeighborList(j, nIndx) = .true.
            else
               NeighborList(nIndx, j) = .false.
               NeighborList(j, nIndx) = .false.
            end if
         end if
      end do
      NeighborList(nIndx, nIndx) = .false.

   end subroutine
!=================================================================================
!     This function updates the neighborlist if a move is accepted.
   subroutine NeighborUpdate_Delete(nIndx, nSwapIndx)
      use SimParameters
      use IndexingFunctions
      use Coords
      implicit none
      integer, intent(in) :: nIndx, nSwapIndx
      integer :: nType
      integer :: i

!      write(35,*) nIndx, nSwapIndx
      nType = typeList(nIndx)
!      nSwapIndx = molArray(nType)%mol(NPART(nType))%indx

!      NeighborList(nIndx,:) = NeighborList(nSwapIndx,:)
!      NeighborList(:,nIndx) = NeighborList(:,nSwapIndx)

      if (nIndx .eq. nSwapIndx) then
         NeighborList(nSwapIndx, :) = .false.
         NeighborList(:, nSwapIndx) = .false.
         return
      end if

      do i = 1, maxMol
!       if(.not. isActive(i)) then
!          NeighborList(i,nIndx) = .false.
!          NeighborList(nIndx,i) = .false.
!          cycle
!        endif
         if (NeighborList(i, nSwapIndx)) then
            if (i .ne. nIndx) then
               NeighborList(i, nIndx) = .true.
               NeighborList(nIndx, i) = .true.
            end if
         else
            NeighborList(i, nIndx) = .false.
            NeighborList(nIndx, i) = .false.
         end if
      end do

      NeighborList(nIndx, nIndx) = .false.

      NeighborList(nSwapIndx, :) = .false.
      NeighborList(:, nSwapIndx) = .false.

   end subroutine
!=================================================================================
   subroutine MultipleSwap_EnergyCriteria(nType2, nIndx1, PairList, isIncluded, rejMove)
      use SimParameters
      use Coords
      use IndexingFunctions
      implicit none

      logical, intent(out) :: rejMove
      real(dp), intent(in) :: PairList(1:maxMol)
      logical, intent(in) :: isIncluded(:)
      integer, intent(in) :: nType2, nIndx1

      logical :: neiFlipped, memberAdded
      logical :: ClusterMember(1:maxMol)
      logical :: flipped(1:maxMol)
      integer :: i, j, h
      integer :: jType
      integer :: curNeigh(1:60), neiMax

      rejMove = .false.
      if (NTotal .eq. 1) return
      ClusterMember = .false.
      flipped = .false.

!     This section dermines which molecules are neighbored with the new trial position.  In the event
!     that the molecule's new location has no neghibors all further calcualtions are skipped and the move is
!     rejected.

      memberAdded = .false.
      do j = 1, maxMol
         if (isIncluded(j)) then
            jType = typeList(j)
            if (PairList(j) .le. Eng_Critr(jType, nType2)) then
               ClusterMember(j) = .true.
               memberAdded = .true.
            end if
         end if
      end do

      if (.not. memberAdded) then
         rejMove = .true.
         return
      end if

!      This part of the code tabulates all the neighbors
      neiMax = 0
      curNeigh = 0
      do i = 1, maxMol
         if (NeighborList(i, nIndx1)) then
            if (i .ne. nIndx1) then
               if (isIncluded(i)) then
                  neiMax = neiMax + 1
                  curNeigh(neiMax) = i
               end if
            end if
         end if
      end do

!     This section performs a quick check to see if the molecules that were neighbored with the old position
!     are part of the new cluster.  If all the old neighbors are indeed part of the cluster then no furth
!     calculations are needed.
      neiFlipped = .true.

      do i = 1, neiMax
         if (.not. clusterMember(curNeigh(i))) then
            neiFlipped = .false.
            exit
         end if
      end do

      if (neiFlipped) then
         rejMove = .false.
         return
      end if

      do h = 1, NTotal
         memberAdded = .false.
         do i = 1, maxMol
            if (isIncluded(i)) then
               if (ClusterMember(i) .neqv. flipped(i)) then
                  do j = 1, maxMol
                     if (isIncluded(i)) then
                        if (NeighborList(i, j)) then
                           ClusterMember(j) = .true.
                           memberAdded = .true.
                        end if
                     end if
                  end do
                  flipped(i) = .true.
               end if
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
end module

