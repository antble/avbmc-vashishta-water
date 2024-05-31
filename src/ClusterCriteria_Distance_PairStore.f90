!=================================================================================
module DistanceCriteria_PairStore

   logical, allocatable :: ClusterMember(:), flipped(:)
contains
!=================================================================================
!     Extensive Cluster Criteria Check.  Used at the start and end of the simulation.
!     This ensures that all particles in the cluster are properly connected to each other.
!     This function also calculates the initial Neighborlist that is used throughout the simulation.
   subroutine Detailed_DistanceCriteria(rejMove)
      use Coords
      use IndexingFunctions
      use ParallelVar
      use PairStorage
      use SimParameters
      implicit none
      logical, intent(out) :: rejMove
!      real(dp), intent(inout) :: PairList(:, :)

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

!      do iType=1,nMolTypes
!      do jType=iType,nMolTypes
!        do iMol = 1 ,NPART(iType)
!        iIndx = MolArray(iType)%mol(iMol)%indx
!        do jMol = 1 ,NPART(jType)
!          jIndx = MolArray(jType)%mol(jMol)%indx
!          if(PairList(iIndx,jIndx) .le. Dist_Critr_sq ) then
!            NeighborList(iIndx,jIndx)=.true.
!            NeighborList(jIndx,iIndx)=.true.
!          endif
!        enddo
!        enddo
!      enddo
!      enddo

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
                  if (rPair(globIndx1, globIndx2)%p%r_sq .lt. Dist_Critr_sq) then
!                write(*,*) globIndx1, globIndx2, rPair(globIndx1, globIndx2) % p % r_sq, rPair(globIndx1, globIndx2) % p % r
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
!     This function determines if a given translational move will destroy a cluster.
   subroutine Shift_DistanceCriteria(nIndx, rejMove)
      use Coords
      use IndexingFunctions
      use PairStorage, only: rPairNew, newDist
      use SimParameters
      implicit none

      logical, intent(out) :: rejMove
!      real(dp), intent(in) :: PairList(:)
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
         jType = typeList(jIndx)
         jMol = subIndxList(jIndx)
         globIndx2 = molArray(jType)%mol(jMol)%globalIndx(1)
         if (rPairNew(globIndx1, globIndx2)%p%r_sq .lt. Dist_Critr_sq) then
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
   subroutine NeighborUpdate_Distance_PairStore(nIndx)
      use Coords
      use IndexingFunctions
      use PairStorage
      use SimParameters
      implicit none
      integer, intent(in) :: nIndx
!      real(dp), intent(in) :: PairList(:)
      integer :: nType, nMol, jMol, jType, jIndx, globIndx1, globIndx2

      nType = typeList(nIndx)
      nMol = subIndxList(nIndx)
      globIndx1 = MolArray(nType)%mol(nMol)%globalIndx(1)
!      write(35,*) nIndx
      do jIndx = 1, maxMol
         if (.not. isActive(jIndx)) then
            cycle
         end if
         if (jIndx .ne. nIndx) then
            jType = typeList(jIndx)
            jMol = subIndxList(jIndx)
            globIndx2 = MolArray(jType)%mol(jMol)%globalIndx(1)
!          write(35,*) globIndx1, globIndx2, rPairNew(globIndx1,globIndx2) % p % r_sq, Dist_Critr_sq
!          if( rPair(globIndx1,globIndx2) % p % r_sq  .le.  Dist_Critr_sq ) then
            if (rPairNew(globIndx1, globIndx2)%p%r_sq .lt. Dist_Critr_sq) then
!            write(35,*) globIndx1, globIndx2, rPair(globIndx1,globIndx2) % p % r_sq
               NeighborList(nIndx, jIndx) = .true.
               NeighborList(jIndx, nIndx) = .true.
            else
               NeighborList(nIndx, jIndx) = .false.
               NeighborList(jIndx, nIndx) = .false.
            end if
         end if
      end do
!      NeighborList(nIndx,nIndx) = .false.

!      do jType = 1, nMolTypes
!        do jMol = 1, NPART(jType)
!          jIndx = MolArray(jType)%mol(jMol)%indx
!          if(jIndx .ne. nIndx) then
!            globIndx2 = MolArray(jType)%mol(jMol)%globalIndx(1)
!            if( rPair(globIndx1,globIndx2) % p % r_sq  .le.  Dist_Critr_sq ) then
!              write(2,*) globIndx1, globIndx2, rPair(globIndx1,globIndx2) % p % r_sq
!              NeighborList(nIndx, jIndx) = .true.
!              NeighborList(jIndx, nIndx) = .true.
!            else
!              NeighborList(nIndx, jIndx) = .false.
!              NeighborList(jIndx, nIndx) = .false.
!            endif
!          endif
!        enddo
!      enddo
!      NeighborList(nIndx,nIndx) = .false.
!      write(2,*)

   end subroutine
!=================================================================================
!     This function updates the neighborlist if a move is accepted.
   subroutine NeighborUpdate_SwapIn_Distance_PairStore(nType)
      use Coords
      use IndexingFunctions
      use PairStorage
      use SimParameters
      implicit none
      integer, intent(in) :: nType
!      real(dp), intent(in) :: PairList(:)
      integer ::  nIndx, nMol, jMol, jType, jIndx, globIndx1, globIndx2

      nMol = NPART(nType) + 1
      nIndx = MolArray(nType)%mol(nMol)%indx
      globIndx1 = MolArray(nType)%mol(nMol)%globalIndx(1)
!      write(35,*) nType, nMol, nIndx , globIndx1
      do jType = 1, nMolTypes
         do jMol = 1, NPART(jType)
            jIndx = MolArray(jType)%mol(jMol)%indx
            if (jIndx .ne. nIndx) then
               globIndx2 = MolArray(jType)%mol(jMol)%globalIndx(1)
!            write(35,*) globIndx1, globIndx2, rPairNew(globIndx1,globIndx2) % p % r_sq
               if (rPairNew(globIndx1, globIndx2)%p%r_sq .lt. Dist_Critr_sq) then

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

