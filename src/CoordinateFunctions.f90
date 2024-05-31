!========================================================================
module CoordinateFunctions
   character(len=50) :: configFile
   integer, parameter :: maxLineLen = 500
contains
!========================================================================
   subroutine AllocateCoordinateArrays
      use ForceField
      use SimParameters
      use Coords
      use EnergyTables
      use PairStorage
      implicit none
      integer :: i, j, k, cnt, AllocationStatus
      integer :: iType

      if (nMolTypes .eq. 0) then
         write (*, *) "ERROR! The program attempted to allocate the coordinate arrays without &
&                    prior to number of molecule types being defined!"
         stop
      end if

      allocate (MolArray(1:nMolTypes), stat=AllocationStatus)
      do i = 1, nMolTypes
         allocate (MolArray(i)%mol(1:NMAX(i)), stat=AllocationStatus)
         do j = 1, NMAX(i)
            allocate (MolArray(i)%mol(j)%x(1:nAtoms(i)), stat=AllocationStatus)
            allocate (MolArray(i)%mol(j)%y(1:nAtoms(i)), stat=AllocationStatus)
            allocate (MolArray(i)%mol(j)%z(1:nAtoms(i)), stat=AllocationStatus)
            allocate (MolArray(i)%mol(j)%globalIndx(1:nAtoms(i)), stat=AllocationStatus)
         end do
      end do

      nTotalAtoms = 0
      do iType = 1, nMolTypes
         nTotalAtoms = nTotalAtoms + NMAX(iType)*nAtoms(iType)
      end do
      write (35, *) "Total Atoms in the System:", nTotalAtoms

      allocate (atomIndicies(1:nTotalAtoms), stat=AllocationStatus)
      cnt = 0
      do i = 1, nMolTypes
         do j = 1, NMAX(i)
            do k = 1, nAtoms(i)
               cnt = cnt + 1
               MolArray(i)%mol(j)%globalIndx(k) = cnt
               atomIndicies(cnt)%nType = i
               atomIndicies(cnt)%nMol = j
               atomIndicies(cnt)%nAtom = k
               atomIndicies(cnt)%atmType = atomArray(i, k)
            end do
         end do
      end do

      do i = 1, size(atomIndicies)
         write (35, *) "Atom indices",i, atomIndicies(i)%nType, atomIndicies(i)%nMol, atomIndicies(i)%nAtom, atomIndicies(i)%atmType
      end do

      allocate (gasConfig(1:nMolTypes), stat=AllocationStatus)
      do i = 1, nMolTypes
         allocate (gasConfig(i)%x(1:nAtoms(i)), stat=AllocationStatus)
         allocate (gasConfig(i)%y(1:nAtoms(i)), stat=AllocationStatus)
         allocate (gasConfig(i)%z(1:nAtoms(i)), stat=AllocationStatus)
      end do

      maxMol = sum(NMAX)
      allocate (typeList(1:maxMol), stat=AllocationStatus)
      allocate (subIndxList(1:maxMol), stat=AllocationStatus)
      allocate (ETable(1:maxMol), stat=AllocationStatus)
      allocate (NeiETable(1:maxMol), stat=AllocationStatus)
      allocate (neiCount(1:maxMol), stat=AllocationStatus)

      cnt = 0
      do i = 1, nMolTypes
         do j = 1, NMAX(i)
            cnt = cnt + 1
            typeList(cnt) = i
            subIndxList(cnt) = j
         end do
      end do

      allocate (isActive(1:maxMol), stat=AllocationStatus)
      allocate (NeighborList(1:maxMol, 1:maxMol), stat=AllocationStatus)
      maxAtoms = maxval(nAtoms)

      allocate (newMol%x(1:maxAtoms), stat=AllocationStatus)
      allocate (newMol%y(1:maxAtoms), stat=AllocationStatus)
      allocate (newMol%z(1:maxAtoms), stat=AllocationStatus)

      allocate (newMol2%x(1:maxAtoms), stat=AllocationStatus)
      allocate (newMol2%y(1:maxAtoms), stat=AllocationStatus)
      allocate (newMol2%z(1:maxAtoms), stat=AllocationStatus)

      call CreateDistArrays

   end subroutine
!=============================================================================
   subroutine CreateJointArray
      use SimParameters
      use Coords
      implicit none
      integer :: i, j
      integer :: cnt

      cnt = 0
!      write(35,*) "Molecular Index List"
!      write(35,*) "Type, Mol, Index"
      do i = 1, nMolTypes
         do j = 1, NMAX(i)
            cnt = cnt + 1
!          JointArray(cnt)%molType = i
            MolArray(i)%mol(j)%indx = cnt
!          write(35,*) i, j, cnt
!          do jj=1,nAtoms(i)
!            JointArray(cnt)%x(jj) => MolArray(i)%mol(j)%x(jj)
!            JointArray(cnt)%y(jj) => MolArray(i)%mol(j)%y(jj)
!            JointArray(cnt)%z(jj) => MolArray(i)%mol(j)%z(jj)
!          enddo
         end do
      end do

   end subroutine
!=============================================================================
   subroutine ReadInitialConfiguration(fileName)
      use ParallelVar
      use SimParameters
      use ForceField
      use Coords
      implicit none
      character(len=maxLineLen), intent(in) :: fileName
      integer :: iType, iMol, iAtom, iIndx
      real(dp) :: x, y, z
      character(len=2) :: atmSymbol
      character(len=maxLineLen) :: format_string, fl_name, out1
      character(len=maxLineLen) :: test
      ! print *, "INITIAL CONFIGURATION FILENAME=", fileName 
      if (multipleInput) then
         if (myid .lt. 10) then
            format_string = "(A,I1,A)"
         elseif (myid .lt. 100) then
            format_string = "(A,I2,A)"
         elseif (myid .lt. 1000) then
            format_string = "(A,I3,A)"
         elseif (myid .lt. 1000) then
            format_string = "(A,I4,A)"
         else
            format_string = "(A,I5,A)"
         end if
         write (fl_name, format_string) fileName, myid, ".dat"
         open (unit=10, file=trim(adjustl(fl_name)), status="Old")
      else
         open (unit=10, file=fileName, status="Old")
      end if
      
      read (10, *)
      read (10, *) (NPART(iType), iType=1, nMolTypes)
      ! print *, "NPART", (NPART(iType), iType=1, nMolTypes)
!     This block ensures the initial configuration is within the min/max bounds set in the input parameter file
      do iType = 1, nMolTypes
         if (NPART(iType) .gt. NMAX(iType)) then
            write (*, *) "Error! Number of particles in the initial"
            write (*, *) "configuration are greater than the maxmium bounds set in the"
            write (*, *) "input parameters."
            write (*, *) "Molecule Type:", iType
            write (*, *) "Initial Size:", NPART(iType)
            write (*, *) "Maximum Allowed:", NMAX(iType)
            stop "Error! Configuration above maximum particle bounds!"
         end if
      end do
      do iType = 1, nMolTypes
         if (NPART(iType) .lt. NMIN(iType)) then
            write (*, *) "Error! Number of particles in the initial"
            write (*, *) "configuration are less than the minimum bounds set in the"
            write (*, *) "input parameters."
            write (*, *) "Molecule Type:", iType
            write (*, *) "Initial Size:", NPART(iType)
            write (*, *) "Minimum Allowed:", NMIN(iType)
            stop "Error! Configuration below minimum particle bounds!"
         end if
      end do

      do iType = 1, nMolTypes
         do iMol = 1, NPART(iType)
            do iAtom = 1, nAtoms(iType)
               read (10, *) atmSymbol, x, y, z
               MolArray(iType)%mol(iMol)%x(iAtom) = x
               MolArray(iType)%mol(iMol)%y(iAtom) = y
               MolArray(iType)%mol(iMol)%z(iAtom) = z
               ! print *, "atmSymbol, X, Y, Z | ", atmSymbol, x, y, z, NPART(iType), iAtom, iMol, iType 
            end do
         end do
      end do
      print *, "config file was properly closed"
      close (10)

      NTotal = sum(NPART)
      isActive = .false.
      do iType = 1, nMolTypes
         do iMol = 1, NPART(iType)
            iIndx = MolArray(iType)%mol(iMol)%indx
            isActive(iIndx) = .true.
         end do
      end do
   ! print *, "config file was properly closed"
   end subroutine
!========================================================
   subroutine ReadInitialGasPhase(fileName)
      use SimParameters
      use ForceField
      use Units
      use Coords
      implicit none
      character(len=maxLineLen), intent(in) :: fileName
      integer :: iType, iAtom
      character(len=2) :: atmSymbol
      real(dp) :: x1, y1, z1

      open (unit=10, file=fileName, status="Old")
      do iType = 1, nMolTypes
         read (10, *)
         do iAtom = 1, nAtoms(iType)
            read (10, *) atmSymbol, gasConfig(iType)%x(iAtom), gasConfig(iType)%y(iAtom), gasConfig(iType)%z(iAtom)
         end do
      end do
      close (10)
      ! print *, "GAS PHASE CONFIGURATION FILENAME=", fileName 
!     Shift the coordinates such that the first atom is located at (0,0,0) for convinence.
      do iType = 1, nMolTypes
         x1 = gasConfig(iType)%x(1)
         y1 = gasConfig(iType)%y(1)
         z1 = gasConfig(iType)%z(1)
         do iAtom = 1, nAtoms(iType)
            gasConfig(iType)%x(iAtom) = gasConfig(iType)%x(iAtom) - x1
            gasConfig(iType)%y(iAtom) = gasConfig(iType)%y(iAtom) - y1
            gasConfig(iType)%z(iAtom) = gasConfig(iType)%z(iAtom) - z1
         end do
      end do
   end subroutine
!========================================================
   subroutine RecenterCoordinates
      use SimParameters
      use ForceField
      use Coords
      implicit none
      integer :: iType, iMol, iAtom
      real(dp) :: xcm, ycm, zcm

      do iType = 1, nMolTypes
         if (nMax(iType) .gt. 0) then
            xcm = MolArray(iType)%mol(1)%x(1)
            ycm = MolArray(iType)%mol(1)%y(1)
            zcm = MolArray(iType)%mol(1)%z(1)
            exit
         end if
      end do

      do iType = 1, nMolTypes
         do iMol = 1, NPART(iType)
            do iAtom = 1, nAtoms(iType)
               MolArray(iType)%mol(iMol)%x(iAtom) = MolArray(iType)%mol(iMol)%x(iAtom) - xcm
               MolArray(iType)%mol(iMol)%y(iAtom) = MolArray(iType)%mol(iMol)%y(iAtom) - ycm
               MolArray(iType)%mol(iMol)%z(iAtom) = MolArray(iType)%mol(iMol)%z(iAtom) - zcm
            end do
         end do
      end do

   end subroutine
!========================================================
   subroutine Generate_UnitSphere(x, y, z)
      use VarPrecision
      implicit none
      real(dp), intent(out) :: x, y, z
      real(dp) :: u_12_sq, u1, u2, grnd

      u_12_sq = 2E0
      do while (u_12_sq .ge. 1)
         u1 = 2E0*grnd() - 1E0
         u2 = 2E0*grnd() - 1E0
         u_12_sq = u1*u1 + u2*u2
      end do

      x = 2E0*u1*sqrt(1E0 - u_12_sq)
      y = 2E0*u2*sqrt(1E0 - u_12_sq)
      z = (1E0 - 2E0*u_12_sq)

   end subroutine
!==========================================================================
!     The purpose of this function is to generate a random position for an atom
!     given a fixed bond angle (bond_ang) with a given vector (v1) and a given distance (r2).
!     The torsional rotational angle is chosen randomly from 0 to 2pi.
!     The coordinate is created using a relative orthonormal framework given by these vectors
!     w1=(x1,y1,z1)   w2=(-y1,x1,0)  w3=(-x1*z1, -y1*z1, x1^2 + y1^2)
!     Using these vectors, the new vector(v2) is calculated using a rotational matrix

   subroutine Generate_UnitCone(v1, r2, bond_ang, v2)
      use Constants
      use CoordinateTypes
      implicit none
      type(SimpleAtomCoords), intent(in) :: v1
      type(SimpleAtomCoords), intent(out) :: v2
      real(dp), intent(in) :: bond_ang, r2
      real(dp) :: tors_angle
      real(dp) :: r1
      real(dp) :: s_term, c_term
      real(dp) :: coeff1, coeff2, coeff3
      real(dp) :: r_proj, grnd

      r1 = v1%x*v1%x + v1%y*v1%y + v1%z*v1%z
      r1 = sqrt(r1)

      tors_angle = two_pi*grnd()

      s_term = sin(tors_angle)
      c_term = cos(tors_angle)
      r_proj = sqrt(v1%x*v1%x + v1%y*v1%y)

      coeff1 = (r2/r1)*cos(bond_ang)
      coeff2 = (r2/r_proj)*sin(bond_ang)
      coeff3 = coeff2/r1

      v2%x = coeff1*v1%x - coeff2*c_Term*v1%y - coeff3*s_Term*v1%x*v1%z
      v2%y = coeff1*v1%y + coeff2*c_term*v1%x - coeff3*s_term*v1%y*v1%z
      v2%z = coeff1*v1%z + coeff3*s_term*(r_proj*r_proj)

   end subroutine
   
!==========================================================================
!     The purpose of this function is to generate two random position for two atoms
!     given two fixed bond angles (bond_ang1 and bond_ang2), one dihedral angle (dihed) which
!     is the angle between two planes made by two angles with a given vector (v1) and two
!     given distances (r2 and r3).
!     The coordinate is created using a relative orthonormal framework given by these vectors
!     w1=(x1,y1,z1)   w2=(-y1,x1,0)  w3=(-x1*z1, -y1*z1, x1^2 + y1^2)
!     Using these vectors, the new vectors(v2 and v3) is calculated using a rotational matrix

   subroutine Generate_UnitPyramid(v1, r2, r3, bond_ang1, bond_ang2, dihed, v2, v3)
      use Constants
      use CoordinateTypes
      implicit none
      type(SimpleAtomCoords), intent(in) :: v1
      type(SimpleAtomCoords), intent(out) :: v2, v3
      real(dp), intent(in) :: bond_ang1, bond_ang2, dihed, r2, r3
      real(dp) :: tors_angle
      real(dp) :: r1
      real(dp) :: s_term, c_term
      real(dp) :: coeff1, coeff2, coeff3
      real(dp) :: r_proj, grnd

      r1 = v1%x*v1%x + v1%y*v1%y + v1%z*v1%z
      r1 = sqrt(r1)

      tors_angle = two_pi*grnd()

      s_term = sin(tors_angle)
      c_term = cos(tors_angle)
      r_proj = sqrt(v1%x*v1%x + v1%y*v1%y)

      coeff1 = (r2/r1)*cos(bond_ang1)
      coeff2 = (r2/r_proj)*sin(bond_ang1)
      coeff3 = coeff2/r1

      v2%x = coeff1*v1%x - coeff2*c_Term*v1%y - coeff3*s_Term*v1%x*v1%z
      v2%y = coeff1*v1%y + coeff2*c_term*v1%x - coeff3*s_term*v1%y*v1%z
      v2%z = coeff1*v1%z + coeff3*s_term*(r_proj*r_proj)

      tors_angle = tors_angle + dihed
      s_term = sin(tors_angle)
      c_term = cos(tors_angle)

      coeff1 = (r3/r1)*cos(bond_ang2)
      coeff2 = (r3/r_proj)*sin(bond_ang2)
      coeff3 = coeff2/r1

      v3%x = coeff1*v1%x - coeff2*c_Term*v1%y - coeff3*s_Term*v1%x*v1%z
      v3%y = coeff1*v1%y + coeff2*c_term*v1%x - coeff3*s_term*v1%y*v1%z
      v3%z = coeff1*v1%z + coeff3*s_term*(r_proj*r_proj)

   end subroutine

!==========================================================================
!     The purpose of this function is similar to the UnitCone function however
!     in this case the torsional angle is fixed in addition to the distance and
!     bond angle which means there is exactly a single point in space which
!     corresponds to these constraints.
!     The coordinate is created using a relative orthonormal framework given by these vectors
!     w1=(x2,y2,z2)   w2=(-y2,x2,0)  w3=(-x2*z2, -y2*z2, x2^2 + y2^2)
   subroutine Generate_UnitTorsion(v1, v2, r3, bond_ang, tors_angle, v3)
      use Constants
      use CoordinateTypes
      implicit none
      type(SimpleAtomCoords), intent(in) :: v1, v2
      type(SimpleAtomCoords), intent(out) :: v3
      real(dp), intent(in) :: bond_ang, tors_angle, r3
      real(dp) :: x1_u, y1_u, z1_u
      real(dp) :: x1_s, y1_s, z1_s
      real(dp) :: r2, rot_angle
      real(dp) :: s_term, c_term
      real(dp) :: coeff1, coeff2, coeff3
      real(dp) :: r_proj

      r2 = v2%x*v2%x + v2%y*v2%y + v2%z*v2%z
      r2 = sqrt(r2)
      r_proj = sqrt(v2%x*v2%x + v2%y*v2%y)

!
      x1_u = v1%x - v2%x
      y1_u = v1%y - v2%y
      z1_u = v1%z - v2%z

!     Calculate the v1 vector's w2 and w3 components in the new orthonormal framework.
!      x1_s =  ( v2%x*x1_u + v2%y*y1_u + v2%z * z1_u)/(r2)
      y1_s = (-v2%y*x1_u + v2%x*y1_u)/r_proj
      z1_s = (-v2%x*v2%z*x1_u - v2%y*v2%z*y1_u + r_proj*r_proj*z1_u)/(r_proj*r2)

!     Calculate the torsional rotation angle for the new v3 vector from the v1 components
      rot_angle = atan2(z1_s, y1_s)
!      write(2,*) rot_angle*180E0/pi, tors_angle*180E0/pi, (rot_angle + tors_angle)*180E0/pi
      rot_angle = rot_angle + tors_angle

!     Rescale the angle
!      do while(rot_angle .lt. 0E0)
!        rot_angle = rot_angle + two_pi
!      enddo
!      do while(rot_angle .gt. two_pi)
!        rot_angle = rot_angle - two_pi
!      enddo

      s_term = sin(rot_angle)
      c_term = cos(rot_angle)

      coeff1 = (r3/r2)*cos(bond_ang)
      coeff2 = (r3/r_proj)*sin(bond_ang)
      coeff3 = coeff2/r2

      v3%x = coeff1*v2%x - coeff2*c_Term*v2%y - coeff3*s_Term*v2%x*v2%z
      v3%y = coeff1*v2%y + coeff2*c_Term*v2%x - coeff3*s_Term*v2%y*v2%z
      v3%z = coeff1*v2%z + coeff3*s_Term*(r_proj*r_proj)

   end subroutine
!==========================================================================
!     The purpose of this function is to generate three random position for three atoms
!     given three fixed bond angles (bond_ang1 and bond_ang2 and bond_ang3),
!     two dihedral angle (dihed1 and dihed2) where dihed1 is the angle between two planes made by
!     first and second angles and dihed2 is the angle between two planes made by second and third angles
!     with a given vector (v1) and three given distances (r2 and r3 and r4).
!     The coordinate is created using a relative orthonormal framework given by these vectors
!     w1=(x1,y1,z1)   w2=(-y1,x1,0)  w3=(-x1*z1, -y1*z1, x1^2 + y1^2)
!     Using these vectors, the new vectors(v2 and v3) is calculated using a rotational matrix

   subroutine Generate_UnitTetrahedral(v1, r2, r3, r4, bond_ang1, bond_ang2, bond_ang3, dihed1, dihed2, v2, v3, v4)
      use Constants
      use CoordinateTypes
      implicit none
      type(SimpleAtomCoords), intent(in) :: v1
      type(SimpleAtomCoords), intent(out) :: v2, v3, v4
      real(dp), intent(in) :: bond_ang1, bond_ang2, bond_ang3, dihed1, dihed2, r2, r3, r4
      real(dp) :: tors_angle
      real(dp) :: r1
      real(dp) :: s_term, c_term
      real(dp) :: coeff1, coeff2, coeff3
      real(dp) :: r_proj, grnd

      r1 = v1%x*v1%x + v1%y*v1%y + v1%z*v1%z
      r1 = sqrt(r1)

      tors_angle = two_pi*grnd()

      s_term = sin(tors_angle)
      c_term = cos(tors_angle)
      r_proj = sqrt(v1%x*v1%x + v1%y*v1%y)

      coeff1 = (r2/r1)*cos(bond_ang1)
      coeff2 = (r2/r_proj)*sin(bond_ang1)
      coeff3 = coeff2/r1

      v2%x = coeff1*v1%x - coeff2*c_Term*v1%y - coeff3*s_Term*v1%x*v1%z
      v2%y = coeff1*v1%y + coeff2*c_term*v1%x - coeff3*s_term*v1%y*v1%z
      v2%z = coeff1*v1%z + coeff3*s_term*(r_proj*r_proj)

      tors_angle = tors_angle + dihed1
      s_term = sin(tors_angle)
      c_term = cos(tors_angle)

      coeff1 = (r3/r1)*cos(bond_ang2)
      coeff2 = (r3/r_proj)*sin(bond_ang2)
      coeff3 = coeff2/r1

      v3%x = coeff1*v1%x - coeff2*c_Term*v1%y - coeff3*s_Term*v1%x*v1%z
      v3%y = coeff1*v1%y + coeff2*c_term*v1%x - coeff3*s_term*v1%y*v1%z
      v3%z = coeff1*v1%z + coeff3*s_term*(r_proj*r_proj)

      tors_angle = tors_angle + dihed2
      s_term = sin(tors_angle)
      c_term = cos(tors_angle)

      coeff1 = (r4/r1)*cos(bond_ang3)
      coeff2 = (r4/r_proj)*sin(bond_ang3)
      coeff3 = coeff2/r1

      v4%x = coeff1*v1%x - coeff2*c_Term*v1%y - coeff3*s_Term*v1%x*v1%z
      v4%y = coeff1*v1%y + coeff2*c_term*v1%x - coeff3*s_term*v1%y*v1%z
      v4%z = coeff1*v1%z + coeff3*s_term*(r_proj*r_proj)

   end subroutine

!==========================================================================
   subroutine GenerateRotation(nTor, nToGrow, vTor, vPrev, rGrow, bendingAngles, dihedralAngles, torsTypes, vGrow, wTor)
      use SimParameters
      use Coords
!      use CoordinateFunctions
      use ForceField
      use ForceFieldFunctions
      use Constants
      use Rosenbluth_Functions_LJ_Q
      use CBMC_Variables
      use CBMC_Utility
      implicit none

      integer, intent(in) :: nTor, nToGrow, torsTypes(:, :)
      real(dp), intent(in) :: rGrow(:), bendingAngles(:), dihedralAngles(:)
      real(dp), intent(out) :: wTor
      type(SimpleAtomCoords), intent(in) :: vTor(:), vPrev
      type(SimpleAtomCoords), intent(out) :: vGrow(:)

      integer :: iGrow, iTor, iRosen, nTorsion, nSel
      real(dp) :: E_Trial(1:nRosenTorsion), ProbRosen(1:nRosenTorsion), rot_Trial(1:nRosenTorsion), rot_Angle, Prob, eng
      real(dp) :: rPrev, r_proj, rTor(1:3), azimuth_Grow(1:3), azimuth_Grow_0(1:3), azimuth_Tor(1:3), torsAngle
      real(dp) :: grnd, ranNum, sumInt, x1_u, y1_u, z1_u, x1_s, y1_s, z1_s, s_term, c_term, coeff1, coeff2, coeff3

      rPrev = vPrev%x*vPrev%x + vPrev%y*vPrev%y + vPrev%z*vPrev%z
      rPrev = sqrt(rPrev)
      r_proj = sqrt(vPrev%x*vPrev%x + vPrev%y*vPrev%y)
      azimuth_Grow_0(1) = 0d0
      if (nToGrow .ge. 2) then
         azimuth_Grow_0(2) = dihedralAngles(1)
         if (nToGrow .eq. 3) azimuth_Grow_0(3) = dihedralAngles(1) + dihedralAngles(2)
      end if
      nTorsion = nTor*nToGrow
      select case (nTorsion)
      case (0) !There is no torsion angle because nTor = 0, so a random rotation angle is generated
         rot_Angle = two_pi*grnd()
         wTor = dble(nRosenTorsion)
      case (1) !This part of the molecule is linear because nTor = nToGrow = 1, so the torsion angle is generated using rejection method
         x1_u = vTor(1)%x - vPrev%x
         y1_u = vTor(1)%y - vPrev%y
         z1_u = vTor(1)%z - vPrev%z

         y1_s = (-vPrev%y*x1_u + vPrev%x*y1_u)/r_proj
         z1_s = (-vPrev%x*vPrev%z*x1_u - vPrev%y*vPrev%z*y1_u + r_proj*r_proj*z1_u)/(r_proj*rPrev)

         azimuth_Tor(1) = atan2(z1_s, y1_s)
         call GenerateTorsAngle(torsAngle, torsTypes(1, 1), Prob)
         rot_Angle = azimuth_Tor(1) + torsAngle
         wTor = dble(nRosenTorsion)
      case default
         do iTor = 1, nTor
            x1_u = vTor(iTor)%x - vPrev%x
            y1_u = vTor(iTor)%y - vPrev%y
            z1_u = vTor(iTor)%z - vPrev%z

            y1_s = (-vPrev%y*x1_u + vPrev%x*y1_u)/r_proj
            z1_s = (-vPrev%x*vPrev%z*x1_u - vPrev%y*vPrev%z*y1_u + r_proj*r_proj*z1_u)/(r_proj*rPrev)

            azimuth_Tor(iTor) = atan2(z1_s, y1_s)
         end do
         E_Trial = 0d0
         do iRosen = 1, nRosenTorsion
            rot_Trial(iRosen) = two_pi*grnd()
            do iGrow = 1, nToGrow
               azimuth_Grow(iGrow) = azimuth_Grow_0(iGrow) + rot_Trial(iRosen)
               if (azimuth_Grow(iGrow) .ge. two_pi) then
                  iTor = floor(azimuth_Grow(iGrow)/two_pi)
                  azimuth_Grow(iGrow) = azimuth_Grow(iGrow) - dble(iTor)*two_pi
               end if
            end do
            do iGrow = 1, nToGrow
               do iTor = 1, nTor
                  torsAngle = azimuth_Grow(iGrow) - azimuth_Tor(iTor)
                  if (torsAngle .le. 0d0) torsAngle = torsAngle + two_pi
                  eng = Trappe_CosNx(torsAngle, torsData(torsTypes(iGrow, iTor))%a)
                  E_Trial(iRosen) = E_Trial(iRosen) + eng
               end do
            end do
            ProbRosen(iRosen) = exp(-beta*E_Trial(iRosen))
         end do
         wTor = sum(ProbRosen)
         ranNum = grnd()*wTor
         sumInt = ProbRosen(1)
         nSel = 1
         do while (sumInt .lt. ranNum)
            nSel = nSel + 1
            sumInt = sumInt + ProbRosen(nSel)
         end do
         rot_Angle = rot_Trial(nSel)
      end select

      do iGrow = 1, nToGrow
         azimuth_Grow(iGrow) = azimuth_Grow_0(iGrow) + rot_Angle

         s_term = sin(azimuth_Grow(iGrow))
         c_term = cos(azimuth_Grow(iGrow))

         coeff1 = (rGrow(iGrow)/rPrev)*cos(bendingAngles(iGrow))
         coeff2 = (rGrow(iGrow)/r_proj)*sin(bendingAngles(iGrow))
         coeff3 = coeff2/rPrev

         vGrow(iGrow)%x = coeff1*vPrev%x - coeff2*c_Term*vPrev%y - coeff3*s_Term*vPrev%x*vPrev%z
         vGrow(iGrow)%y = coeff1*vPrev%y + coeff2*c_Term*vPrev%x - coeff3*s_Term*vPrev%y*vPrev%z
         vGrow(iGrow)%z = coeff1*vPrev%z + coeff3*s_Term*(r_proj*r_proj)
      end do

   end subroutine
!==========================================================================

   subroutine GenerateRotation_Reverse(nTor, nToGrow, vTor, vPrev, torsTypes, vGrow, wTor)
      use SimParameters
      use Coords
      use ForceField
      use ForceFieldFunctions
      use Constants
      use Rosenbluth_Functions_LJ_Q
      use CBMC_Variables
      use CBMC_Utility
      implicit none

      integer, intent(in) :: nTor, nToGrow, torsTypes(:, :)
      real(dp), intent(out) :: wTor
      type(SimpleAtomCoords), intent(in) :: vTor(:), vPrev, vGrow(:)

      integer :: iGrow, iTor, iRosen, nTorsion, nSel
      real(dp) :: dihedralAngles(1:6)
      real(dp) :: E_Trial(1:nRosenTorsion), ProbRosen(1:nRosenTorsion), rot_Trial(1:nRosenTorsion), rot_Angle, Prob, eng
      real(dp) :: rPrev, r_proj, rTor(1:3), azimuth_Grow(1:3), azimuth_Grow_0(1:3), azimuth_Tor(1:3), torsAngle
      real(dp) :: grnd, ang1, ang2, x1_u, y1_u, z1_u, x1_s, y1_s, z1_s, s_term, c_term, coeff1, coeff2, coeff3

      nTorsion = nTor*nToGrow
      select case (nTorsion)
      case (0) !There is no torsion angle because nTor = 0
         wTor = dble(nRosenTorsion)
      case (1) !This part of the molecule is linear because nTor = nToGrow = 1
         wTor = dble(nRosenTorsion)
      case default
         rPrev = vPrev%x*vPrev%x + vPrev%y*vPrev%y + vPrev%z*vPrev%z
         rPrev = sqrt(rPrev)
         r_proj = sqrt(vPrev%x*vPrev%x + vPrev%y*vPrev%y)
         x1_u = vGrow(1)%x - vPrev%x
         y1_u = vGrow(1)%y - vPrev%y
         z1_u = vGrow(1)%z - vPrev%z
         y1_s = (-vPrev%y*x1_u + vPrev%x*y1_u)/r_proj
         z1_s = (-vPrev%x*vPrev%z*x1_u - vPrev%y*vPrev%z*y1_u + r_proj*r_proj*z1_u)/(r_proj*rPrev)
         rot_Trial(1) = atan2(z1_s, y1_s)
         if (nToGrow .ge. 2) then
            x1_u = vGrow(2)%x - vPrev%x
            y1_u = vGrow(2)%y - vPrev%y
            z1_u = vGrow(2)%z - vPrev%z
            y1_s = (-vPrev%y*x1_u + vPrev%x*y1_u)/r_proj
            z1_s = (-vPrev%x*vPrev%z*x1_u - vPrev%y*vPrev%z*y1_u + r_proj*r_proj*z1_u)/(r_proj*rPrev)
            ang1 = atan2(z1_s, y1_s)
            dihedralAngles(1) = ang1 - rot_Trial(1)
            if (dihedralAngles(1) .le. 0d0) dihedralAngles(1) = dihedralAngles(1) + two_pi
            if (nToGrow .eq. 3) then
               x1_u = vGrow(3)%x - vPrev%x
               y1_u = vGrow(3)%y - vPrev%y
               z1_u = vGrow(3)%z - vPrev%z
               y1_s = (-vPrev%y*x1_u + vPrev%x*y1_u)/r_proj
               z1_s = (-vPrev%x*vPrev%z*x1_u - vPrev%y*vPrev%z*y1_u + r_proj*r_proj*z1_u)/(r_proj*rPrev)
               ang2 = atan2(z1_s, y1_s)
               dihedralAngles(2) = ang2 - ang1
               if (dihedralAngles(2) .le. 0d0) dihedralAngles(2) = dihedralAngles(2) + two_pi
            end if
         end if
         azimuth_Grow_0(1) = 0d0
         if (nToGrow .ge. 2) then
            azimuth_Grow_0(2) = dihedralAngles(1)
            if (nToGrow .eq. 3) azimuth_Grow_0(3) = dihedralAngles(1) + dihedralAngles(2)
         end if
         do iTor = 1, nTor
            x1_u = vTor(iTor)%x - vPrev%x
            y1_u = vTor(iTor)%y - vPrev%y
            z1_u = vTor(iTor)%z - vPrev%z
            y1_s = (-vPrev%y*x1_u + vPrev%x*y1_u)/r_proj
            z1_s = (-vPrev%x*vPrev%z*x1_u - vPrev%y*vPrev%z*y1_u + r_proj*r_proj*z1_u)/(r_proj*rPrev)
            azimuth_Tor(iTor) = atan2(z1_s, y1_s)
         end do
         E_Trial = 0d0
         do iRosen = 2, nRosenTorsion
            rot_Trial(iRosen) = two_pi*grnd()
         end do
         do iRosen = 1, nRosenTorsion
            do iGrow = 1, nToGrow
               azimuth_Grow(iGrow) = azimuth_Grow_0(iGrow) + rot_Trial(iRosen)
               if (azimuth_Grow(iGrow) .ge. two_pi) then
                  iTor = floor(azimuth_Grow(iGrow)/two_pi)
                  azimuth_Grow(iGrow) = azimuth_Grow(iGrow) - dble(iTor)*two_pi
               end if
            end do
            do iGrow = 1, nToGrow
               do iTor = 1, nTor
                  torsAngle = azimuth_Grow(iGrow) - azimuth_Tor(iTor)
                  if (torsAngle .le. 0d0) torsAngle = torsAngle + two_pi
                  eng = Trappe_CosNx(torsAngle, torsData(torsTypes(iGrow, iTor))%a)
                  E_Trial(iRosen) = E_Trial(iRosen) + eng
               end do
            end do
            ProbRosen(iRosen) = exp(-beta*E_Trial(iRosen))
         end do
         wTor = sum(ProbRosen)
      end select

   end subroutine

!==========================================================================
   subroutine GenerateBondLength(dist, k_bond, r_eq, ProbGen)
      use Constants
      use SimParameters
      use ForceFieldFunctions
      use AcceptRates, only: distGen_accpt, distGen_atmp
      implicit none
      logical acpt
      real(dp), intent(in) :: k_bond, r_eq
      real(dp), intent(out) :: dist
      real(dp), intent(out), optional :: ProbGen
      real(dp) :: normalSigma, grnd, eng
!      real(dp) ::

      if (k_bond .eq. 0E0) then
         dist = r_eq
         ProbGen = 1E0
         return
      end if

      normalSigma = 1E0/(0.5E0*k_bond*beta)
      acpt = .false.
      do while (.not. acpt)
         distGen_atmp = distGen_atmp + 1E0
         dist = -log(1E0 - grnd())
!        dist = sqrt(-2E0*log(grnd())) * cos(two_pi*grnd())
!        dist = dist*normalSigma + r_eq
         eng = Harmonic(dist, k_bond, r_eq)
         ProbGen = dist*dist*exp(-beta*eng + dist)
!        ProbGen = dist * dist
         if (ProbGen/3E0 .gt. grnd()) then
            acpt = .true.
         end if
      end do
      distGen_accpt = distGen_accpt + 1E0
   end subroutine

   subroutine GenerateBondLength_Vashishta(dist, r_max, ProbGen)
      use Constants
      use ForceFieldPara_Vashishta
      use SimParameters, only : beta, rc, r1s, outputEConv
      use ForceFieldFunctions
      use ForceField
      use AcceptRates, only: distGen_accpt, distGen_atmp
      implicit none
      logical acpt
      real(dp), intent(in) :: r_max
      real(dp), intent(out) :: dist
      real(dp), intent(out), optional :: ProbGen
      real(dp) :: Ubond, grnd, eng

      integer(kind=atomIntType) :: atmType1, atmType2
      real(dp) :: H, eta, q, D, W, r4s, r, rmin_ij
      real(dp) :: steric, coulomb, chargeDipole, vanderWaals
      real(dp) :: tmp, rcheta, rcqr1s, rcDr4s, rcW, qO, qH 

      steric = 0E0_dp 
      coulomb = 0E0_dp 
      chargeDipole = 0E0_dp 
      vanderWaals = 0E0_dp 

      atmType1 = 1 
      atmtype2 = 2 
      qO = q_tab(atmType1, atmType2)
      qH = q_tab(atmType2, atmType2)
      H = H_tab(atmType1, atmType2)
      eta = eta_tab(atmType1, atmType2)
      r4s = r4s_tab(atmType1, atmType2)
      D = D_tab(atmType1, atmType2)
      W = W_tab(atmType1, atmType2)
      rmin_ij = r_min_tab(atmType1, atmType2)

      acpt = .false.
      do while (.not. acpt)
         distGen_atmp = distGen_atmp + 1E0
         dist = -log(1E0 - grnd())
         r = dist
         if (r*r > rmin_ij .and. r < r_max) then
            !============================================================
            rcheta = H/(rc**eta)                         
            tmp = -eta*rcheta/rc                         
            steric = H/(r**eta) - rcheta - (r-rc)*tmp
            rcqr1s = (qH*qO/rc)*exp(-rc/r1s) 
            tmp = -rcqr1s*(1/rc + 1/r1s)
            coulomb = (qH*qO/r)*exp(-r/r1s) - rcqr1s - (r-rc)*tmp
            rcDr4s = -D*exp(-rc/r4s)/(rc**4) 
            tmp = rcDr4s*(4/rc + 1/r4s)
            chargeDipole = - D*exp(-r/r4s)/(r**4) - rcDr4s + (r-rc)*tmp
            rcW = -W/(rc**6) 
            tmp = -rcW*6/rc
            vanderWaals = - W/(r**6) - rcW - (r - rc)*tmp 
            Ubond = steric + coulomb + chargeDipole + vanderWaals
            Ubond = Ubond*outputEConv
            !============================================================
            ProbGen = dist*dist*exp(-beta*Ubond  + dist)
            if (ProbGen/3 .gt. grnd()) then
               acpt = .true.
            end if
         end if
      end do
      distGen_accpt = distGen_accpt + 1E0
   end subroutine


   subroutine GenerateBendAngle_Vashishta(angle, r_ij, r_ik)
      use Constants
      use SimParameters
      use ForceField
      use ForceFieldFunctions
      use AcceptRates, only: angGen_accpt, angGen_atmp
      implicit none

      logical :: acpt
      integer :: nSel, dummy
      real(dp), intent(in) :: r_ij, r_ik
      real(dp), intent(out) :: angle
      real(dp) :: k_bend, theta_eq
      real(dp) :: grnd, eng, ranNum, sumInt
      real(dp) :: ProbSel, accptConstant
      real(dp) :: ang_eq, k_eq, st
      real(dp) :: angle1, angle2, angle3
      real(dp) :: eng1, eng2, eng3
      real(dp) :: prob1, prob2, prob3, probMax
      real(dp) :: totalProb, norm
      real(dp) :: maxima, accConstant, ProbGen
      real(dp) :: B, xi, r0, costheta, C, Ujik
      integer :: iAngle, iBin
      integer :: maxBin, startBin
      real(dp) :: Prob(1:nBendHistBins)
      character(len=100) :: format_string, fl_name
      !=============================================================================
      ! integrate bend angle probability distribution 
      !3-BODY INTERACTIONS parameters
      B  = atomDataV(1)%B
      xi = atomDataV(1)%xi
      r0 = atomDataV(1)%r0
      C = atomDataV(1)%C 
      costheta = atomDataV(1)%theta


      ! Using rejection sampling -- simplest way to sample 
      acpt = .false.
      do while (.not. acpt) 
         angle = (theta_max - theta_min)*grnd() + theta_min 
         eng1 = B*exp((xi/(r_ij-r0))+(xi/(r_ik-r0)))*(cos(angle)-costheta)**2 
         eng1 = eng1*outputEConv
         prob1 = sin(angle)*exp(-beta*eng1)
         if (prob1 > grnd()) then
            acpt = .true. 
         end if 
      end do 
      
      ! Prob = 0d0
      ! do iBin = 1, nBendHistBins
      !    angle1 = (iBin - 1)*bendBinWidth
      !    angle2 = (iBin - 0.5d0)*bendBinWidth
      !    angle3 = iBin*bendBinWidth
      !    eng1 = B*exp((xi/(r_ij-r0))+(xi/(r_ik-r0)))*(cos(angle1)-costheta)**2
      !    eng2 = B*exp((xi/(r_ij-r0))+(xi/(r_ik-r0)))*(cos(angle2)-costheta)**2
      !    eng3 = B*exp((xi/(r_ij-r0))+(xi/(r_ik-r0)))*(cos(angle3)-costheta)**2
      !    prob1 = sin(angle1)*exp(-beta*eng1)
      !    prob2 = sin(angle2)*exp(-beta*eng2)
      !    prob3 = sin(angle3)*exp(-beta*eng3)
      !    if (max(prob1, prob3)/probMax .lt. 1d-1) then
      !       totalProb = max(prob1, prob3)*bendBinWidth
      !    else
      !       totalProb = (bendBinWidth/6d0)*(prob1 + 4d0*prob2 + prob3)
      !    end if
      !    Prob(iBin)  = totalProb
      ! end do
      ! norm = sum(Prob)
      ! do iBin = 1, nBendHistBins
      !    Prob(iBin) = Prob(iBin)/norm
      ! end do
      ! accConstant = 0d0
      ! do iBin = 1, nBendHistBins
      !    angle1 = (iBin - 1)*bendBinWidth
      !    angle2 = iBin*bendBinWidth
      !    eng1 = B*exp((xi/(r_ij-r0))+(xi/(r_ik-r0)))*(cos(angle1)-costheta)**2
      !    eng2 = B*exp((xi/(r_ij-r0))+(xi/(r_ik-r0)))*(cos(angle2)-costheta)**2
      !    prob1 = sin(angle1)*exp(-beta*eng1)
      !    prob2 = sin(angle2)*exp(-beta*eng2)
      !    prob1 = max(prob1, prob2)
      !    totalProb = bendBinWidth*prob1/Prob(iBin)
      !    if (totalProb .ge. accConstant) then
      !       accConstant = totalProb
      !    end if
      ! end do

      ! ! eng1 = B*exp((xi/(r_ij-r0))+(xi/(r_ik-r0)))*(cos(angle1)-costheta)**2
      ! ! prob1 = sin(maxima)*!exp(-beta*eng1)
      ! ! totalProb = bendBinWidth*prob1/bendData(iAngle)%Prob(maxBin)

      ! ! find the peak of the distribution
      ! ! format_string = "(A,I1,A)"
      ! ! write (fl_name, format_string) "AngleBendProb.txt"
      ! ! open (unit=111, file=trim(adjustl(fl_name)))
      ! ! do iBin =1, nBendHistBins
      ! !    write (111, *) Prob(iBin) 
      ! ! end do
      ! ! angle = pi*grnd()
      ! accptConstant = bendBinWidth/accConstant
      ! startBin = 0
      ! do iBin = 2, nBendHistBins
      !    Prob(iBin) = Prob(iBin) + Prob(iBin - 1)
      !    if (startBin .eq. 0) then
      !       if (Prob(iBin) .ge. startProb) then
      !          startBin = iBin - 1
      !       end if
      !    end if
      ! end do
      ! !=============================================================================

      ! acpt = .false.
      ! do while (acpt .eqv. .false.)
      !    ranNum = grnd()
      !    !Since the majority of the probability density will be centered near the equilibrium angle
      !    !we can use this fact to find the correct bin with much fewer comparisons by simply starting
      !    !near the beginning of the gaussian curve unless the random number chosen is sufficiently small.
      !    if (ranNum .lt. startProb) then
      !       nSel = 1
      !    else
      !       nSel = startBin
      !    end if

      !    do while (Prob(nSel) .lt. ranNum)
      !       nSel = nSel + 1
      !    end do

      !    !Now that the bin has been chosen, select an angle uniformly from the bin and calculate the
      !    !acceptance probability.
      !    angle = (dble(nSel) - grnd())*bendBinWidth
      !    if (nSel .ne. 1) then
      !       ProbSel = Prob(nSel) - Prob(nSel - 1)
      !    else
      !       ProbSel =Prob(nSel)
      !    end if

      !    eng = B*exp((xi/(r_ij-r0))+(xi/(r_ik-r0)))*(cos(angle)-costheta)**2
      !    ProbGen = sin(angle)*exp(-beta*eng)/ProbSel
      !    ProbGen = ProbGen*accptConstant
      !    if (ProbGen .gt. grnd()) then
      !       acpt = .true.
      !    end if
      ! end do
      ! print *, angle, angle*180/pi
   end subroutine
!==========================================================================
   subroutine GenerateBendAngle(angle, bendType, ProbGen)
      use Constants
      use SimParameters
      use ForceField
      use ForceFieldFunctions
      use AcceptRates, only: angGen_accpt, angGen_atmp
      implicit none

      interface
         subroutine SelectBin_IntTable(integralTable, nBin)
            use VarPrecision
            implicit none
            real(dp), intent(in) :: integralTable(:)
            integer, intent(out) :: nBin
         end subroutine
      end interface

      logical :: acpt
      integer :: nSel, dummy
      integer, intent(in) :: bendType
      real(dp), intent(out) :: angle
      real(dp), intent(out), optional :: ProbGen
      real(dp) :: k_bend, theta_eq
      real(dp) :: grnd, eng, ranNum, sumInt
      real(dp) :: ProbSel

      acpt = .false.
      do while (acpt .eqv. .false.)
         angGen_atmp = angGen_atmp + 1E0

!         angle = pi*grnd()
!         angle = acos(1E0-2E0*grnd())

         ranNum = grnd()
         !Since the majority of the probability density will be centered near the equilibrium angle
         !we can use this fact to find the correct bin with much fewer comparisons by simply starting
         !near the beginning of the gaussian curve unless the random number chosen is sufficiently small.
         if (ranNum .lt. startProb) then
            nSel = 1
         else
            nSel = bendData(bendType)%startBin
         end if
         do while (bendData(bendType)%Prob(nSel) .lt. ranNum)
            nSel = nSel + 1
         end do

         !Now that the bin has been chosen, select an angle uniformly from the bin and calculate the
         !acceptance probability.
         angle = (dble(nSel) - grnd())*bendBinWidth
         if (nSel .ne. 1) then
            ProbSel = bendData(bendType)%Prob(nSel) - bendData(bendType)%Prob(nSel - 1)
         else
            ProbSel = bendData(bendType)%Prob(nSel)
         end if

         eng = Harmonic(angle, k_bend, theta_eq)
         ProbGen = sin(angle)*exp(-beta*eng)/ProbSel
         ProbGen = ProbGen*bendData(bendType)%accptConstant
         if (ProbGen .gt. grnd()) then
            acpt = .true.
         end if
      end do
      angGen_accpt = angGen_accpt + 1E0

   end subroutine
!==========================================================================
   subroutine GenerateTwoBranches(ang1, ang2, dihed, bendType1, bendType2, bendType3, wBending)
      use Constants
      use SimParameters
      use ForceFieldFunctions
      use ForceField
      use AcceptRates, only: dihedGen_accpt, dihedGen_atmp
      use CBMC_Variables
      use RandomTools
      implicit none

      logical acpt
      integer, intent(in) :: bendType1, bendType2, bendType3
      real(dp) :: ang1, ang2, dihed, wBending
      integer :: bin, nSel, cnt, iTrial
      real(dp) :: k_bend1, theta_eq1, k_bend2, theta_eq2, k_bend3, theta_eq3
      real(dp) :: ang1_val(1:nRosenTwoBranch), ang2_val(1:nRosenTwoBranch), dihed_val(1:nRosenTwoBranch)
      real(dp) :: wei(1:nRosenTwoBranch), ang3, std1, std2, std3, ang3Min, ang3Max
      real(dp) :: grnd, ranNum, sumInt, U1, U2, NG1, NG2, Gaussian
      real(dp) :: Sang1, Cang1, Sang2, Cang2, Sang3, Cang3, Sdihed, Cdihed

      k_bend1 = bendData(bendType1)%k_eq
      theta_eq1 = bendData(bendType1)%ang_eq
      k_bend2 = bendData(bendType2)%k_eq
      theta_eq2 = bendData(bendType2)%ang_eq
      k_bend3 = bendData(bendType3)%k_eq
      theta_eq3 = bendData(bendType3)%ang_eq
      std1 = 1d0/sqrt(beta*k_bend1)
      std2 = 1d0/sqrt(beta*k_bend2)
      std3 = 1d0/sqrt(beta*k_bend3)

      if ((k_bend1 .eq. 0E0) .or. (k_bend2 .eq. 0E0) .or. (k_bend3 .eq. 0E0)) then
         ang1 = theta_eq1
         ang2 = theta_eq2
         ang3 = theta_eq3
         dihed = acos((cos(ang3) - cos(ang1)*cos(ang2))/(sin(ang1)*sin(ang2)))
         wBending = 1E0
         return
      end if

      cnt = 0

      do iTrial = 1, nRosenTwoBranch
         acpt = .false.
         do while (acpt .eqv. .false.)
            ang1_val(iTrial) = Gaussian()*std1 + theta_eq1
            ang2_val(iTrial) = Gaussian()*std2 + theta_eq2
            ang3Min = abs(ang1_val(iTrial) - ang2_val(iTrial))
            ang3Max = two_pi - (ang1_val(iTrial) + ang2_val(iTrial))
            if (ang3Max .ge. pi) ang3Max = two_pi - ang3Max
            ang3 = Gaussian()*std3 + theta_eq3
            if ((ang3 .ge. ang3Min) .and. (ang3 .le. ang3Max)) then
               if ((ang1_val(iTrial) .gt. 0d0) .and. (ang1_val(iTrial) .lt. pi)) then
                  if ((ang2_val(iTrial) .gt. 0d0) .and. (ang2_val(iTrial) .lt. pi)) then
                     acpt = .true.
                  end if
               end if
            end if
            dihedGen_atmp = dihedGen_atmp + 1E0
         end do
         Sang1 = sin(ang1_val(iTrial))
         Cang1 = cos(ang1_val(iTrial))
         Sang2 = sin(ang2_val(iTrial))
         Cang2 = cos(ang2_val(iTrial))
         Sang3 = sin(ang3)
         Cang3 = cos(ang3)
         Cdihed = (Cang3 - Cang1*Cang2)/(Sang1*Sang2)
         if (Cdihed .le. -1d0) Cdihed = -1d0
         if (Cdihed .ge. 1d0) Cdihed = 1d0
         Sdihed = sqrt(1d0 - Cdihed*Cdihed)
         dihed_val(iTrial) = acos(Cdihed)
         if (grnd() .ge. 0.5d0) dihed_val(iTrial) = two_pi - dihed_val(iTrial)
         wei(iTrial) = abs(Sang3/Sdihed)
      end do
      wBending = sum(wei)
      ranNum = grnd()*wBending
      sumInt = 0d0
      nSel = 0
      do while (sumInt .lt. ranNum)
         nSel = nSel + 1
         ang1 = ang1_val(nSel)
         ang2 = ang2_val(nSel)
         dihed = dihed_val(nSel)
         sumInt = sumInt + wei(nSel)
      end do
      dihedGen_accpt = dihedGen_accpt + 1E0

   end subroutine
!==========================================================================

   subroutine GenerateTwoBranches_Reverse(v1, v3, v4, bendType1, bendType2, bendType3, wBending)
      use Constants
      use SimParameters
      use ForceFieldFunctions
      use ForceField
      use AcceptRates, only: dihedGen_accpt, dihedGen_atmp
      use CBMC_Variables
      use RandomTools
      use CoordinateTypes
      implicit none

      logical acpt
      type(SimpleAtomCoords), intent(in) :: v1, v3, v4
      integer, intent(in) :: bendType1, bendType2, bendType3
      real(dp) :: wBending
      integer :: bin, nSel, cnt, iTrial
      real(dp) :: k_bend1, theta_eq1, k_bend2, theta_eq2, k_bend3, theta_eq3
      real(dp) :: ang1_val(1:nRosenTwoBranch), ang2_val(1:nRosenTwoBranch), dihed_val(1:nRosenTwoBranch)
      real(dp) :: wei(1:nRosenTwoBranch), ang3, std1, std2, std3, ang3Min, ang3Max
      real(dp) :: grnd, ranNum, sumInt, U1, U2, NG1, NG2, Gaussian
      real(dp) :: r1, r3, r4, Sang1, Cang1, Sang2, Cang2, Sang3, Cang3, Sdihed, Cdihed

      k_bend1 = bendData(bendType1)%k_eq
      theta_eq1 = bendData(bendType1)%ang_eq
      k_bend2 = bendData(bendType2)%k_eq
      theta_eq2 = bendData(bendType2)%ang_eq
      k_bend3 = bendData(bendType3)%k_eq
      theta_eq3 = bendData(bendType3)%ang_eq
      std1 = 1d0/sqrt(beta*k_bend1)
      std2 = 1d0/sqrt(beta*k_bend2)
      std3 = 1d0/sqrt(beta*k_bend3)

      if ((k_bend1 .eq. 0E0) .or. (k_bend2 .eq. 0E0) .or. (k_bend3 .eq. 0E0)) then
         wBending = 1E0
         return
      end if

      r1 = sqrt(v1%x*v1%x + v1%y*v1%y + v1%z*v1%z)
      r3 = sqrt(v3%x*v3%x + v3%y*v3%y + v3%z*v3%z)
      r4 = sqrt(v4%x*v4%x + v4%y*v4%y + v4%z*v4%z)
      Cang1 = (v1%x*v3%x + v1%y*v3%y + v1%z*v3%z)/(r1*r3)
      Cang2 = (v1%x*v4%x + v1%y*v4%y + v1%z*v4%z)/(r1*r4)
      Cang3 = (v3%x*v4%x + v3%y*v4%y + v3%z*v4%z)/(r3*r4)
      if (Cang1 .le. -1d0) Cang1 = -1d0
      if (Cang1 .ge. 1d0) Cang1 = 1d0
      if (Cang2 .le. -1d0) Cang2 = -1d0
      if (Cang2 .ge. 1d0) Cang2 = 1d0
      if (Cang3 .le. -1d0) Cang3 = -1d0
      if (Cang3 .ge. 1d0) Cang3 = 1d0
      Sang1 = sqrt(1d0 - Cang1*Cang1)
      Sang2 = sqrt(1d0 - Cang2*Cang2)
      Sang3 = sqrt(1d0 - Cang3*Cang3)
      Cdihed = (Cang3 - Cang1*Cang2)/(Sang1*Sang2)
      if (Cdihed .le. -1d0) Cdihed = -1d0
      if (Cdihed .ge. 1d0) Cdihed = 1d0
      Sdihed = sqrt(1d0 - Cdihed*Cdihed)
      iTrial = 1
      wei(iTrial) = abs(Sang3/Sdihed)
      wBending = wBending + wei(iTrial)

      do iTrial = 2, nRosenTwoBranch
         acpt = .false.
         do while (acpt .eqv. .false.)
            ang1_val(iTrial) = Gaussian()*std1 + theta_eq1
            ang2_val(iTrial) = Gaussian()*std2 + theta_eq2
            ang3Min = abs(ang1_val(iTrial) - ang2_val(iTrial))
            ang3Max = two_pi - (ang1_val(iTrial) + ang2_val(iTrial))
            if (ang3Max .ge. pi) ang3Max = two_pi - ang3Max
            ang3 = Gaussian()*std3 + theta_eq3
            if ((ang3 .ge. ang3Min) .and. (ang3 .le. ang3Max)) then
               if ((ang1_val(iTrial) .gt. 0d0) .and. (ang1_val(iTrial) .lt. pi)) then
                  if ((ang2_val(iTrial) .gt. 0d0) .and. (ang2_val(iTrial) .lt. pi)) then
                     acpt = .true.
                  end if
               end if
            end if
         end do
         Sang1 = sin(ang1_val(iTrial))
         Cang1 = cos(ang1_val(iTrial))
         Sang2 = sin(ang2_val(iTrial))
         Cang2 = cos(ang2_val(iTrial))
         Sang3 = sin(ang3)
         Cang3 = cos(ang3)
         Cdihed = (Cang3 - Cang1*Cang2)/(Sang1*Sang2)
         if (Cdihed .le. -1d0) Cdihed = -1d0
         if (Cdihed .ge. 1d0) Cdihed = 1d0
         Sdihed = sqrt(1d0 - Cdihed*Cdihed)
         dihed_val(iTrial) = acos(Cdihed)
         if (grnd() .ge. 0.5d0) dihed_val(iTrial) = two_pi - dihed_val(iTrial)
         wei(iTrial) = abs(Sang3/Sdihed)
         wBending = wBending + wei(iTrial)
      end do

   end subroutine
!========================================================================================

   subroutine GenerateThreeBranches(ang1, ang2, ang3, dihed1, dihed2, &
                                    bendType1, bendType2, bendType3, bendType4, bendType5, bendType6, wBending)
      use Constants
      use SimParameters
      use ForceFieldFunctions
      use ForceField
      use AcceptRates, only: dihedGen_accpt, dihedGen_atmp
      use CBMC_Variables
      use RandomTools
      implicit none

      logical acpt
      integer, intent(in) :: bendType1, bendType2, bendType3, bendType4, bendType5, bendType6
      real(dp) :: ang1, ang2, ang3, dihed1, dihed2
      real(dp) :: wBending
      integer :: bin, nSel, cnt, iTrial
      real(dp) :: k_bend1, theta_eq1, k_bend2, theta_eq2, k_bend3, theta_eq3
      real(dp) :: k_bend4, theta_eq4, k_bend5, theta_eq5, k_bend6, theta_eq6, std1, std2, std3, std4, std5
      real(dp) :: grnd, eng, ang4, ang5, ang6, dihed3
      real(dp) :: ang1_val(1:nRosenThreeBranch), ang2_val(1:nRosenThreeBranch), ang3_val(1:nRosenThreeBranch)
      real(dp) :: dihed1_val(1:nRosenThreeBranch), dihed2_val(1:nRosenThreeBranch), wei(1:nRosenThreeBranch)
      real(dp) :: ranNum, sumInt, weight, U1, U2, NG1, NG2, ang4Min, ang4Max, ang5Min, ang5Max, Gaussian
      real(dp) :: Sang1, Cang1, Sang2, Cang2, Sang3, Cang3, Sang4, Cang4, Sang5, Cang5, Cang6, Sdihed1, Cdihed1, Sdihed2, Cdihed2

      k_bend1 = bendData(bendType1)%k_eq
      theta_eq1 = bendData(bendType1)%ang_eq
      k_bend2 = bendData(bendType2)%k_eq
      theta_eq2 = bendData(bendType2)%ang_eq
      k_bend3 = bendData(bendType3)%k_eq
      theta_eq3 = bendData(bendType3)%ang_eq
      k_bend4 = bendData(bendType4)%k_eq
      theta_eq4 = bendData(bendType4)%ang_eq
      k_bend5 = bendData(bendType5)%k_eq
      theta_eq5 = bendData(bendType5)%ang_eq
      k_bend6 = bendData(bendType6)%k_eq
      theta_eq6 = bendData(bendType6)%ang_eq
      std1 = 1d0/sqrt(beta*k_bend1)
      std2 = 1d0/sqrt(beta*k_bend2)
      std3 = 1d0/sqrt(beta*k_bend3)
      std4 = 1d0/sqrt(beta*k_bend4)
      std5 = 1d0/sqrt(beta*k_bend5)

      if ((k_bend1 .eq. 0E0) .or. (k_bend2 .eq. 0E0) .or. (k_bend3 .eq. 0E0) .or. (k_bend4 .eq. 0E0) &
          .or. (k_bend5 .eq. 0E0) .or. (k_bend6 .eq. 0E0)) then
         ang1 = theta_eq1
         ang2 = theta_eq2
         ang3 = theta_eq3
         ang4 = theta_eq4
         ang5 = theta_eq5
         ang6 = theta_eq6
         dihed1 = acos((cos(ang4) - cos(ang1)*cos(ang2))/(sin(ang1)*sin(ang2)))
         dihed2 = acos((cos(ang5) - cos(ang2)*cos(ang3))/(sin(ang2)*sin(ang3)))
         wBending = 1E0
         return
      end if

      cnt = 0
      do iTrial = 1, nRosenThreeBranch
         acpt = .false.
         do while (acpt .eqv. .false.)
            ang1_val(iTrial) = Gaussian()*std1 + theta_eq1
            ang2_val(iTrial) = Gaussian()*std2 + theta_eq2
            ang4Min = abs(ang1_val(iTrial) - ang2_val(iTrial))
            ang4Max = two_pi - (ang1_val(iTrial) + ang2_val(iTrial))
            if (ang4Max .ge. pi) ang4Max = two_pi - ang4Max
            ang4 = Gaussian()*std4 + theta_eq4
            if ((ang4 .ge. ang4Min) .and. (ang4 .le. ang4Max)) then
               if ((ang1_val(iTrial) .gt. 0d0) .and. (ang1_val(iTrial) .lt. pi)) then
                  if ((ang2_val(iTrial) .gt. 0d0) .and. (ang2_val(iTrial) .lt. pi)) then
                     ang3_val(iTrial) = Gaussian()*std3 + theta_eq3
                     ang5Min = abs(ang2_val(iTrial) - ang3_val(iTrial))
                     ang5Max = two_pi - (ang2_val(iTrial) + ang3_val(iTrial))
                     if (ang5Max .ge. pi) ang5Max = two_pi - ang5Max
                     ang5 = Gaussian()*std5 + theta_eq5
                     if ((ang5 .ge. ang5Min) .and. (ang5 .le. ang5Max)) then
                        if ((ang3_val(iTrial) .gt. 0d0) .and. (ang3_val(iTrial) .lt. pi)) then
                           Sang1 = sin(ang1_val(iTrial))
                           Cang1 = cos(ang1_val(iTrial))
                           Sang2 = sin(ang2_val(iTrial))
                           Cang2 = cos(ang2_val(iTrial))
                           Sang3 = sin(ang3_val(iTrial))
                           Cang3 = cos(ang3_val(iTrial))
                           Sang4 = sin(ang4)
                           Cang4 = cos(ang4)
                           Sang5 = sin(ang5)
                           Cang5 = cos(ang5)
                           Cdihed1 = (Cang4 - Cang1*Cang2)/(Sang1*Sang2)
                           if (Cdihed1 .le. -1d0) Cdihed1 = -1d0
                           if (Cdihed1 .ge. 1d0) Cdihed1 = 1d0
                           Sdihed1 = sqrt(1d0 - Cdihed1*Cdihed1)
                           dihed1_val(iTrial) = acos(Cdihed1)
                           if (grnd() .ge. 0.5d0) dihed1_val(iTrial) = two_pi - dihed1_val(iTrial)
                           Cdihed2 = (Cang5 - Cang2*Cang3)/(Sang2*Sang3)
                           if (Cdihed2 .le. -1d0) Cdihed2 = -1d0
                           if (Cdihed2 .ge. 1d0) Cdihed2 = 1d0
                           Sdihed2 = sqrt(1d0 - Cdihed2*Cdihed2)
                           dihed2_val(iTrial) = acos(Cdihed2)
                           if (grnd() .ge. 0.5d0) dihed2_val(iTrial) = two_pi - dihed2_val(iTrial)
                           dihed3 = two_pi - (dihed1_val(iTrial) + dihed2_val(iTrial))
                           if (dihed3 .le. 0d0) dihed3 = dihed3 + two_pi
                           Cang6 = Cang1*Cang3 + Sang1*Sang3*cos(dihed3)
                           if (Cang6 .le. -1d0) Cang6 = -1d0
                           if (Cang6 .ge. 1d0) Cang6 = 1d0
                           ang6 = acos(Cang6)
                           eng = Harmonic(ang6, k_bend6, theta_eq6)
                           if (exp(-beta*eng) .gt. grnd()) acpt = .true.
                        end if
                     end if
                  end if
               end if
            end if
            dihedGen_atmp = dihedGen_atmp + 1E0
         end do
         wei(iTrial) = abs((Sang4*Sang5)/(Sang2*Sdihed1*Sdihed2))
         wBending = wBending + wei(iTrial)
      end do
      ranNum = grnd()*wBending
      sumInt = 0d0
      nSel = 0
      do while (sumInt .lt. ranNum)
         nSel = nSel + 1
         ang1 = ang1_val(nSel)
         ang2 = ang2_val(nSel)
         ang3 = ang3_val(nSel)
         dihed1 = dihed1_val(nSel)
         dihed2 = dihed2_val(nSel)
         sumInt = sumInt + wei(nSel)
      end do
      dihedGen_accpt = dihedGen_accpt + 1E0

   end subroutine
!========================================================================================

subroutine GenerateThreeBranches_Reverse(v1, v3, v4, v5, bendType1, bendType2, bendType3, bendType4, bendType5, bendType6, wBending)
      use Constants
      use SimParameters
      use ForceFieldFunctions
      use ForceField
      use AcceptRates, only: dihedGen_accpt, dihedGen_atmp
      use CBMC_Variables
      use RandomTools
      use CoordinateTypes
      implicit none

      logical acpt
      type(SimpleAtomCoords), intent(in)  :: v1, v3, v4, v5
      integer, intent(in) :: bendType1, bendType2, bendType3, bendType4, bendType5, bendType6
      real(dp) :: wBending
      integer :: bin, nSel, cnt, iTrial
      real(dp) :: k_bend1, theta_eq1, k_bend2, theta_eq2, k_bend3, theta_eq3
      real(dp) :: k_bend4, theta_eq4, k_bend5, theta_eq5, k_bend6, theta_eq6, std1, std2, std3, std4, std5
      real(dp) :: grnd, eng, ang4, ang5, ang6, dihed3
      real(dp) :: ang1_val(1:nRosenThreeBranch), ang2_val(1:nRosenThreeBranch), ang3_val(1:nRosenThreeBranch)
      real(dp) :: dihed1_val(1:nRosenThreeBranch), dihed2_val(1:nRosenThreeBranch), wei(1:nRosenThreeBranch)
      real(dp) :: ranNum, sumInt, weight, U1, U2, NG1, NG2, ang4Min, ang4Max, ang5Min, ang5Max, Gaussian
      real(dp) :: Sang1, Cang1, Sang2, Cang2, Sang3, Cang3, Sang4, Cang4, Sang5, Cang5, Cang6, Sdihed1, Cdihed1, Sdihed2, Cdihed2
      real(dp) :: r1, r3, r4, r5

      k_bend1 = bendData(bendType1)%k_eq
      theta_eq1 = bendData(bendType1)%ang_eq
      k_bend2 = bendData(bendType2)%k_eq
      theta_eq2 = bendData(bendType2)%ang_eq
      k_bend3 = bendData(bendType3)%k_eq
      theta_eq3 = bendData(bendType3)%ang_eq
      k_bend4 = bendData(bendType4)%k_eq
      theta_eq4 = bendData(bendType4)%ang_eq
      k_bend5 = bendData(bendType5)%k_eq
      theta_eq5 = bendData(bendType5)%ang_eq
      k_bend6 = bendData(bendType6)%k_eq
      theta_eq6 = bendData(bendType6)%ang_eq
      std1 = 1d0/sqrt(beta*k_bend1)
      std2 = 1d0/sqrt(beta*k_bend2)
      std3 = 1d0/sqrt(beta*k_bend3)
      std4 = 1d0/sqrt(beta*k_bend4)
      std5 = 1d0/sqrt(beta*k_bend5)

      if ((k_bend1 .eq. 0E0) .or. (k_bend2 .eq. 0E0) .or. (k_bend3 .eq. 0E0) .or. (k_bend4 .eq. 0E0) &
         .or. (k_bend5 .eq. 0E0) .or. (k_bend6 .eq. 0E0)) then
         wBending = 1E0
         return
      end if

      r1 = sqrt(v1%x*v1%x + v1%y*v1%y + v1%z*v1%z)
      r3 = sqrt(v3%x*v3%x + v3%y*v3%y + v3%z*v3%z)
      r4 = sqrt(v4%x*v4%x + v4%y*v4%y + v4%z*v4%z)
      r5 = sqrt(v5%x*v5%x + v5%y*v5%y + v5%z*v5%z)
      Cang1 = (v1%x*v3%x + v1%y*v3%y + v1%z*v3%z)/(r1*r3)
      Cang2 = (v1%x*v4%x + v1%y*v4%y + v1%z*v4%z)/(r1*r4)
      Cang3 = (v1%x*v5%x + v1%y*v5%y + v1%z*v5%z)/(r1*r5)
      Cang4 = (v3%x*v4%x + v3%y*v4%y + v3%z*v4%z)/(r3*r4)
      Cang5 = (v4%x*v5%x + v4%y*v5%y + v4%z*v5%z)/(r4*r5)
      if (Cang1 .le. -1d0) Cang1 = -1d0
      if (Cang1 .ge. 1d0) Cang1 = 1d0
      if (Cang2 .le. -1d0) Cang2 = -1d0
      if (Cang2 .ge. 1d0) Cang2 = 1d0
      if (Cang3 .le. -1d0) Cang3 = -1d0
      if (Cang3 .ge. 1d0) Cang3 = 1d0
      if (Cang4 .le. -1d0) Cang4 = -1d0
      if (Cang4 .ge. 1d0) Cang4 = 1d0
      if (Cang5 .le. -1d0) Cang5 = -1d0
      if (Cang5 .ge. 1d0) Cang5 = 1d0
      Sang1 = sqrt(1d0 - Cang1*Cang1)
      Sang2 = sqrt(1d0 - Cang2*Cang2)
      Sang3 = sqrt(1d0 - Cang3*Cang3)
      Sang4 = sqrt(1d0 - Cang4*Cang4)
      Sang5 = sqrt(1d0 - Cang5*Cang5)
      Cdihed1 = (Cang4 - Cang1*Cang2)/(Sang1*Sang2)
      if (Cdihed1 .le. -1d0) Cdihed1 = -1d0
      if (Cdihed1 .ge. 1d0) Cdihed1 = 1d0
      Sdihed1 = sqrt(1d0 - Cdihed1*Cdihed1)
      Cdihed2 = (Cang5 - Cang2*Cang3)/(Sang2*Sang3)
      if (Cdihed2 .le. -1d0) Cdihed2 = -1d0
      if (Cdihed2 .ge. 1d0) Cdihed2 = 1d0
      Sdihed2 = sqrt(1d0 - Cdihed2*Cdihed2)
      iTrial = 1
      wei(iTrial) = abs((Sang4*Sang5)/(Sang2*Sdihed1*Sdihed2))
      wBending = wBending + wei(iTrial)

      do iTrial = 2, nRosenThreeBranch
         acpt = .false.
         do while (acpt .eqv. .false.)
            ang1_val(iTrial) = Gaussian()*std1 + theta_eq1
            ang2_val(iTrial) = Gaussian()*std2 + theta_eq2
            ang4Min = abs(ang1_val(iTrial) - ang2_val(iTrial))
            ang4Max = two_pi - (ang1_val(iTrial) + ang2_val(iTrial))
            if (ang4Max .ge. pi) ang4Max = two_pi - ang4Max
            ang4 = Gaussian()*std4 + theta_eq4
            if ((ang4 .ge. ang4Min) .and. (ang4 .le. ang4Max)) then
               if ((ang1_val(iTrial) .gt. 0d0) .and. (ang1_val(iTrial) .lt. pi)) then
                  if ((ang2_val(iTrial) .gt. 0d0) .and. (ang2_val(iTrial) .lt. pi)) then
                     ang3_val(iTrial) = Gaussian()*std3 + theta_eq3
                     ang5Min = abs(ang2_val(iTrial) - ang3_val(iTrial))
                     ang5Max = two_pi - (ang2_val(iTrial) + ang3_val(iTrial))
                     if (ang5Max .ge. pi) ang5Max = two_pi - ang5Max
                     ang5 = Gaussian()*std5 + theta_eq5
                     if ((ang5 .ge. ang5Min) .and. (ang5 .le. ang5Max)) then
                        if ((ang3_val(iTrial) .gt. 0d0) .and. (ang3_val(iTrial) .lt. pi)) then
                           Sang1 = sin(ang1_val(iTrial))
                           Cang1 = cos(ang1_val(iTrial))
                           Sang2 = sin(ang2_val(iTrial))
                           Cang2 = cos(ang2_val(iTrial))
                           Sang3 = sin(ang3_val(iTrial))
                           Cang3 = cos(ang3_val(iTrial))
                           Sang4 = sin(ang4)
                           Cang4 = cos(ang4)
                           Sang5 = sin(ang5)
                           Cang5 = cos(ang5)
                           Cdihed1 = (Cang4 - Cang1*Cang2)/(Sang1*Sang2)
                           if (Cdihed1 .le. -1d0) Cdihed1 = -1d0
                           if (Cdihed1 .ge. 1d0) Cdihed1 = 1d0
                           Sdihed1 = sqrt(1d0 - Cdihed1*Cdihed1)
                           dihed1_val(iTrial) = acos(Cdihed1)
                           if (grnd() .ge. 0.5d0) dihed1_val(iTrial) = two_pi - dihed1_val(iTrial)
                           Cdihed2 = (Cang5 - Cang2*Cang3)/(Sang2*Sang3)
                           if (Cdihed2 .le. -1d0) Cdihed2 = -1d0
                           if (Cdihed2 .ge. 1d0) Cdihed2 = 1d0
                           Sdihed2 = sqrt(1d0 - Cdihed2*Cdihed2)
                           dihed2_val(iTrial) = acos(Cdihed2)
                           if (grnd() .ge. 0.5d0) dihed2_val(iTrial) = two_pi - dihed2_val(iTrial)
                           dihed3 = two_pi - (dihed1_val(iTrial) + dihed2_val(iTrial))
                           if (dihed3 .le. 0d0) dihed3 = dihed3 + two_pi
                           Cang6 = Cang1*Cang3 + Sang1*Sang3*cos(dihed3)
                           if (Cang6 .le. -1d0) Cang6 = -1d0
                           if (Cang6 .ge. 1d0) Cang6 = 1d0
                           ang6 = acos(Cang6)
                           eng = Harmonic(ang6, k_bend6, theta_eq6)
                           if (exp(-beta*eng) .gt. grnd()) acpt = .true.
                        end if
                     end if
                  end if
               end if
            end if
         end do
         wei(iTrial) = abs((Sang4*Sang5)/(Sang2*Sdihed1*Sdihed2))
         wBending = wBending + wei(iTrial)
      end do
   end subroutine
!========================================================================================
   subroutine GenerateTorsAngle(angle, torsType, ProbGen)
      use Constants
      use SimParameters
      use ForceField
      use ForceFieldFunctions
      implicit none
      logical :: acpt
      integer, intent(in) :: torsType
      real(dp), intent(out) :: angle
      real(dp), intent(out), optional :: ProbGen
      real(dp) :: grnd, eng

      acpt = .false.
      do while (acpt .eqv. .false.)
         angle = two_pi*grnd()
         eng = Trappe_CosNx(angle, torsData(torsType)%a)
         ProbGen = exp(-beta*eng)
         if (0.9999E0*ProbGen .gt. grnd()) acpt = .true.
      end do

   end subroutine
!==========================================================================
   subroutine FindBond(nType, mem1, mem2, bondType)
      use Constants
      use SimParameters
      use ForceField
      implicit none
      integer, intent(in) :: nType, mem1, mem2
      integer, intent(out) :: bondType
      integer :: iBond

      do iBond = 1, nBonds(nType)
         if (any(bondArray(nType, iBond)%bondMembr(1:2) .eq. mem1)) then
            if (any(bondArray(nType, iBond)%bondMembr(1:2) .eq. mem2)) then
               bondType = bondArray(nType, iBond)%bondType
               return
            end if
         end if
      end do

      write (6, *) "Error! FindBond function unable to find a bond"
      write (6, *) "containing memebers:", mem1, mem2
      stop

   end subroutine
!==========================================================================
   subroutine FindAngle(nType, mem1, mem2, mem3, bendType)
      use Constants
      use SimParameters
      use ForceField
      implicit none
      integer, intent(in) :: nType, mem1, mem2, mem3
      integer, intent(out) :: bendType
      integer :: iBend

      do iBend = 1, nAngles(nType)
         if (any(bendArray(nType, iBend)%bendMembr(1:3) .eq. mem1)) then
            if (any(bendArray(nType, iBend)%bendMembr(1:3) .eq. mem2)) then
               if (any(bendArray(nType, iBend)%bendMembr(1:3) .eq. mem3)) then
                  bendType = bendArray(nType, iBend)%bendType
                  return
               end if
            end if
         end if
      end do

      write (*, *) "Error! FindAngle function unable to find an angle"
      write (*, *) "containing memebers:", mem1, mem2, mem3
      stop

   end subroutine
!==========================================================================
   subroutine FindTorsion(nType, mem1, mem2, mem3, mem4, torsType)
      use Constants
      use SimParameters
      use ForceField
      implicit none
      integer, intent(in) :: nType, mem1, mem2, mem3, mem4
      integer, intent(out) :: torsType
      integer :: iTors

      do iTors = 1, nTorsional(nType)
         if (any(torsArray(nType, iTors)%torsMembr(1:4) .eq. mem1)) then
            if (any(torsArray(nType, iTors)%torsMembr(1:4) .eq. mem2)) then
               if (any(torsArray(nType, iTors)%torsMembr(1:4) .eq. mem3)) then
                  if (any(torsArray(nType, iTors)%torsMembr(1:4) .eq. mem4)) then
                     torsType = torsArray(nType, iTors)%torsType
                     return
                  end if
               end if
            end if
         end if
      end do

      write (*, *) "Error! FindTorsional function unable to find a torsional angle"
      write (*, *) "containing memebers:", mem1, mem2, mem3, mem4
      stop

   end subroutine
!==========================================================================
   subroutine FindSingleDihedral(nType, hubIndx, dihedType)
      use CBMC_Variables
      implicit none
      integer, intent(in) :: nType, hubIndx
      integer, intent(out) :: dihedType
      integer :: iDihed

      do iDihed = 1, totalDihed
         if (dihedData(iDihed)%molType .eq. nType) then
            if (dihedData(iDihed)%hubIndx .eq. hubIndx) then
               dihedType = iDihed
               return
            end if
         end if
      end do

      write (*, *) "Error! FindSingleDihedral function unable to find a dihedral angle"
      write (*, *) nType, hubIndx
      stop

   end subroutine
!========================================================================
end module
