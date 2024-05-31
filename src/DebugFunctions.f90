!=======================================================
      subroutine DEBUG_Output_NeighborList
         use SimParameters
         use Coords
         implicit none
         integer i, j
         write (35, *) "-------------------------------------------"
         write (35, *) "NeighborList:"
         write (35, *) "NPART=", NPART
         do i = 1, maxMol
            write (35, *) (NeighborList(i, j), j=1, maxMol)
         end do

      end subroutine
!=======================================================
      subroutine DEBUG_NeighborQualityCheck(errRtn)
         use SimParameters
         use Coords
         implicit none
         logical, intent(out) :: errRtn
         logical :: neiFound
         integer :: i, j
         integer :: iType, iMol

         errRtn = .false.
         do i = 1, maxMol
            if (isActive(i)) then
               do j = 1, maxMol
                  if (isActive(j)) then
                     if (NeighborList(i, j) .neqv. NeighborList(j, i)) then
                        write (35, *) "ERROR! Matrix symmetry broken!"
                        write (35, *) i, j, NeighborList(i, j), NeighborList(j, i)
                        errRtn = .true.
                     end if
                  end if
               end do
            end if
         end do

         do i = 1, maxMol
            if (.not. isActive(i)) then
               do j = 1, maxMol
                  if (NeighborList(i, j)) then
                     write (35, *) "ERROR! Inactive particle does not have its neighborlist set to false"
                     write (35, *) i, j, NeighborList(i, j), NeighborList(j, i)
                     errRtn = .true.
                  end if
               end do
            end if
         end do

         if (NTotal .gt. 1) then
            do iType = 1, nMolTypes
               do iMol = 1, NPART(iType)
                  neiFound = .false.
                  i = molArray(iType)%mol(iMol)%indx
                  do j = 1, maxMol
                     if (.not. isActive(j)) then
                        cycle
                     end if
                     if (NeighborList(i, j)) then
                        neiFound = .true.
                        exit
                     end if
                  end do
                  if (.not. neiFound) then
                     write (35, *) "ERROR! Particle has no neighbors!"
                     write (35, *) i, NeighborList(i, :)
                     errRtn = .true.
                  end if
               end do
            end do
         end if

         do iType = 1, nMolTypes
            do iMol = 1, NPART(iType)
               i = molArray(iType)%mol(iMol)%indx
               if (.not. isActive(i)) then
                  write (35, *) "ERROR! Particle is not active when it should be!"
                  write (35, *) i, iType, iMol, isActive(i)
                  errRtn = .true.
               end if
            end do
         end do

         do iType = 1, nMolTypes
            do iMol = NPART(iType) + 1, NMAX(iType)
               i = molArray(iType)%mol(iMol)%indx
               if (isActive(i)) then
                  write (35, *) "ERROR! Particle is active when it shouldn't be!"
                  write (35, *) i, iType, iMol, isActive(i)
                  errRtn = .true.
               end if
            end do
         end do

      end subroutine
!=======================================================
      subroutine DEBUG_Output_NewConfig
         use VarPrecision
         use SimParameters
         use Coords
         use Forcefield
         implicit none
         integer :: iType, iMol, iAtom, atmType
         write (35, *) "-------------------------------------------"
         write (35, *) "Trial Configuration:"
         write (35, *) "NPART=", NPART
         write (35, *)
         do iType = 1, nMolTypes
            do iMol = 1, NPart(iType)
               do iAtom = 1, nAtoms(iType)
                  atmType = atomArray(iType, iAtom)
                  write (35, *) atomData(atmType)%Symb, &
                     MolArray(iType)%mol(iMol)%x(iAtom), &
                     MolArray(iType)%mol(iMol)%y(iAtom), &
                     MolArray(iType)%mol(iMol)%z(iAtom)
               end do
            end do
         end do

         iType = newMol%molType
         do iAtom = 1, nAtoms(iType)
            atmType = atomArray(iType, iAtom)
            write (35, *) atomData(atmType)%Symb, newMol%x(iAtom), newMol%y(iAtom), newMol%z(iAtom)
         end do

      end subroutine
!==================================================================
