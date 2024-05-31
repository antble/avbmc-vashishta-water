



do iRosen = 1, nRosenTrials(nType):
    !..changeMolConfig..
    ...
    do i = 1, nAtoms(nType)
        rosenTrial(iRosen)%x = rosenTrial(iRosen)%x(i) + x1
        rosenTrial(iRosen)%y = rosenTrial(iRosen)%x(i) + y1
        rosenTrial(iRosen)%z = rosenTrial(iRosen)%x(i) + z1
        call Rosen_Mol_New(iRosen, ..., E_Trial(iRosen), overlap)
    end do 
    ...
end do 

E_Max = minval(E_Trial)
ProbRosen = 0d0
do iRosen = 1, nRosenTrials(nType)
    ProbRosen(iRosen) = exp(-beta*(E_Trial(iRosen) - E_Max))
end do

rosenNorm = sum(ProbRosen)

ranNum = grnd()*rosenNorm
sumInt = ProbRosen(1)
nSel = 1
do while (sumInt .lt. ranNum)
    nSel = nSel + 1
    sumInt = sumInt + ProbRosen(nSel)
end do

if (overlap(nSel) .eqv. .true.) then
    rejMove = .true.
    return
end if

! Update the coordinates
do i = 1, nAtoms(nType)
    newMol%x(i) = rosenTrial(nSel)%x(i)
    newMol%y(i) = rosenTrial(nSel)%y(i)
    newMol%z(i) = rosenTrial(nSel)%z(i)
end do

rosenRatio = (ProbRosen(nSel)*dble(nRosenTrials(nType)))/rosenNorm




module Rosenbluth_Functions_Vashishta_Q
   use VarPrecision
contains
!======================================================================================================
   pure subroutine two_body_interaction(H, eta, q, D, W, rc, r1s, r4s, r, steric, coulomb, chargeDipole, vanderWaals)
      use SimParameters, only : outputEConv 
      real(dp), intent(in) :: H, eta, q, D, W, rc, r1s, r4s, r
      real(dp), intent(out) :: steric, coulomb, chargeDipole, vanderWaals
      real(dp) :: tmp, rcheta, rcqr1s, rcDr4s, rcW
      rcheta = H/(rc**eta) ! 2nd term is the truncation
      tmp = -eta*rcheta/rc  ! 3rd term the necessary shift!  
      steric = H/(r**eta) - rcheta - (r-rc)*tmp  
      steric = steric*outputEConv ! eV * (T/eV)
      
      rcqr1s = (q/rc)*exp(-rc/r1s) 
      tmp = -rcqr1s*(1/rc + 1/r1s)
      coulomb = (q/r)*exp(-r/r1s) - rcqr1s - (r-rc)*tmp
      coulomb = coulomb*outputEConv 
      
      rcDr4s = -D*exp(-rc/r4s)/(rc**4) 
      tmp = rcDr4s*(4/rc + 1/r4s)
      chargeDipole = - D*exp(-r/r4s)/(r**4) - rcDr4s + (r-rc)*tmp
      chargeDipole = chargeDipole*outputEConv 

      rcW = -W/(rc**6) 
      tmp = -rcW*6/rc
      vanderWaals = - W/(r**6) - rcW - (r - rc)*tmp 
      vanderWaals = vanderWaals*outputEConv
   end subroutine 

   subroutine threebody_interaction(rx1, rx2, ry1, ry2, rz1, rz2, r_ij, r_ik, B, xi, r0, costheta, threebody)
      use SimParameters, only : outputEConv 
      real(dp), intent(in) :: rx1, rx2, ry1, ry2, rz1, rz2, r_ij, r_ik, B, xi, r0, costheta
      real(dp), intent(out) :: threebody
      real(dp) :: costheta_ijk
      costheta_ijk = ((rx1*rx2 + ry1*ry2 + rz1*rz2)/(r_ij*r_ik)) - costheta
      threebody = B*exp((xi/(r_ij - r0)) + (xi/(r_ik - r0)))*(costheta_ijk**2)
      ! print *, r_ij, r0
      ! print *, "CALCULATED THREEBODY", B, xi, r0, costheta, threebody
      threebody = threebody*outputEConv 
   end subroutine
!======================================================================================================
!      This subrotuine is intended to calculate the Rosenbluth weight for a single trial
!      in any method which regrows an entire molecule for the given trial.
   ! pure subroutine Rosen_Vashishta_Molecule_New(nRosen, nType, included, E_Trial, overlap)
   subroutine Rosen_Vashishta_Molecule_New(nRosen, nType, included, E_Trial, overlap)
      use ForceField, only: nAtoms, r_min_tab, atomArray, atomDataV
      use ForceFieldPara_Vashishta, only: q_tab, H_tab, eta_tab, W_tab, r1s_tab, D_Tab, r4s_tab
      use Coords , only: rosenTrial, MolArray
      use SimParameters, only: nMolTypes, NPART, r1s, rc
      implicit none

      logical, intent(in) :: included(:)
      integer, intent(in) :: nType, nRosen

      logical, intent(out) :: overlap
      real(dp), intent(out) :: E_Trial

      integer :: iAtom, iType, jType, iMol, jMol, kMol, jAtom, kAtom, ctr
      integer :: iIndx, jIndx, kIndx
      integer(kind=atomIntType) :: atmType1, atmType2
      real(dp) :: rx, ry, rz, r, rx1, ry1, rz1, rx2, ry2, rz2
      real(dp) :: rmin_ij
      real(dp) :: H, q, eta, r4s, D, W ! 2-body parameter
      real(dp) :: B, xi, r0, C, theta  ! 3-body parameters
      real(dp) :: O_x1, O_y1, O_z1
      real(dp) :: rcheta, rcqr1s, rcDr4s, rcW, tmp
      real(dp) :: r_ij, r_ik, costheta, costheta_ijk
      real(dp) :: E_steric, E_coulomb, E_chargeDipole, E_vanderWaals, E_twobody, E_threebody
      real(dp) :: steric, coulomb, chargeDipole, vanderWaals, twobody, threebody
      character(len=100) :: fl_name

      E_Trial = 0E0
      overlap = .false.
      
      steric = 0d0 
      coulomb = 0d0 
      chargeDipole = 0d0 
      vanderWaals = 0d0 
      twobody = 0d0
      threebody = 0d0

      E_steric = 0d0 
      E_coulomb = 0d0 
      E_chargeDipole = 0d0 
      E_vanderWaals = 0d0 
      E_twobody = 0d0 
      E_threebody = 0d0 

      
      ! write (fl_name,  "(A,I1,A)") "Energies_Vashishta-Rosen_Vashishta_Molecule_New_", nRosen,".out"
      ! open (unit=35, file=trim(adjustl(fl_name)))
      ! write (35, *) "INTRA/INTER, jMol, iAtom, jAtom, r, E_steric, E_coulomb, E_chargeDipole, E_vanderWaals"
      ! "intramolecular" contribution of the added atom
      do iAtom = 1, nAtoms(nType)-1
         atmType1 = atomArray(nType, iAtom)
         do jType = 1, nMolTypes
            do jAtom = iAtom+1, nAtoms(jType)
               atmType2 = atomArray(jType, jAtom)
               q = q_tab(atmType2, atmType1)
               H = H_tab(atmType2, atmType1)
               eta = eta_tab(atmType2, atmType1)
               r4s = r4s_tab(atmType2, atmType1)
               D = D_tab(atmType2, atmType1)
               W = W_tab(atmType2, atmType1)
               rmin_ij = r_min_tab(atmType2, atmType1)
               
               rx = rosenTrial(nRosen)%x(iAtom) - rosenTrial(nRosen)%x(jAtom)
               ry = rosenTrial(nRosen)%y(iAtom) - rosenTrial(nRosen)%y(jAtom)
               rz = rosenTrial(nRosen)%z(iAtom) - rosenTrial(nRosen)%z(jAtom)
               r = rx*rx + ry*ry + rz*rz
               if (r .lt. rmin_ij) then
                  E_Trial = huge(dp)
                  return
               end if
               r = sqrt(r)
               if (r < rc) then
                  call two_body_interaction(H, eta, q, D, W, rc, r1s, r4s, r, &
                                          steric, coulomb, chargeDipole, vanderWaals)
                  E_steric = E_steric + steric
                  E_coulomb = E_coulomb  + coulomb
                  E_chargeDipole = E_chargeDipole + chargeDipole
                  E_vanderWaals = E_vanderWaals + vanderWaals 
               end if 
            end do
         end do 
      end do 

      ! intermolecular energy contribution 
      do iAtom = 1, nAtoms(nType)
         atmType1 = atomArray(nType, iAtom)
         do jType = 1, nMolTypes
            do jAtom = 1, nAtoms(jType)
               atmType2 = atomArray(jType, jAtom)
               q = q_tab(atmType2, atmType1)
               H = H_tab(atmType2, atmType1)
               eta = eta_tab(atmType2, atmType1)
               r4s = r4s_tab(atmType2, atmType1)
               D = D_tab(atmType2, atmType1)
               W = W_tab(atmType2, atmType1)
               rmin_ij = r_min_tab(atmType2, atmType1)

               do jMol = 1, NPART(jType)
                  jIndx = molArray(jType)%mol(jMol)%indx
                  ! if (included(jIndx) .eqv. .false.) then
                  !    ! print *, "INCLUDED == FALSE", jIndx
                  !    cycle
                  ! end if

                  rx = rosenTrial(nRosen)%x(iAtom) - MolArray(jType)%mol(jMol)%x(jAtom)
                  ry = rosenTrial(nRosen)%y(iAtom) - MolArray(jType)%mol(jMol)%y(jAtom)
                  rz = rosenTrial(nRosen)%z(iAtom) - MolArray(jType)%mol(jMol)%z(jAtom)
                  r = rx*rx + ry*ry + rz*rz
                  if (r .lt. rmin_ij) then
                     E_Trial = huge(dp)
                     overlap = .true.
                     return
                  end if

                  r = sqrt(r)
                  if (r < rc) then
                     call two_body_interaction(H, eta, q, D, W, rc, r1s, r4s, r, &
                                          steric, coulomb, chargeDipole, vanderWaals)
                     E_steric = E_steric + steric
                     E_coulomb = E_coulomb  + coulomb
                     E_chargeDipole = E_chargeDipole + chargeDipole
                     E_vanderWaals = E_vanderWaals + vanderWaals
                  end if 
               end do
            end do
         end do
      end do


      ! 3-BODY INTERACTIONS 
      B  = atomDataV(1)%B
      xi = atomDataV(1)%xi
      r0 = atomDataV(1)%r0
      C = atomDataV(1)%C 
      costheta = atomDataV(1)%theta

      do iType = 1, nMolTypes
         do jType = iType, nMolTypes
            do iMol = 1, NPART(iType) + 1
                  ! Set Oxygen atom as the anchor
                  if (iMol == (NPART(1) + 1)) then 
                     O_x1 = rosenTrial(nRosen)%x(1)
                     O_y1 = rosenTrial(nRosen)%y(1)
                     O_z1 = rosenTrial(nRosen)%z(1)
                     ! print *, "1 Oxygen coordinate:", iMol, O_x1,  O_y1, O_z1
                  else 
                     iIndx = molArray(jType)%mol(iMol)%indx
                     ! if (included(iIndx) .eqv. .false.) then
                     !    cycle
                     ! end if
                     O_x1 = MolArray(iType)%mol(iMol)%x(1)
                     O_y1 = MolArray(iType)%mol(iMol)%y(1)
                     O_z1 = MolArray(iType)%mol(iMol)%z(1)
                     ! print *, "2 Oxygen coordinate:", iMol, O_x1,  O_y1, O_z1
                  end if 
                  do jMol = 1, NPART(1) + 1 
                     jIndx = molArray(jType)%mol(jMol)%indx
                     ! if (included(jIndx) .eqv. .false.) then
                     !    cycle
                     ! end if
                     ! find an O-H that is within r0
                     ! start at 2 since Oxygen=1, note: dont mix up the vashishta forcefield O,H sequence in atom definition
                     do jAtom = 2, nAtoms(jType)
                        if (jMol == NPART(1)+1) then 
                           rx1 = rosenTrial(nRosen)%x(jAtom) - O_x1
                           ry1 = rosenTrial(nRosen)%y(jAtom) - O_y1
                           rz1 = rosenTrial(nRosen)%z(jAtom) - O_z1 
                           r_ij = sqrt(rx1**2 + ry1**2 + rz1**2)
                        else 
                           rx1 = MolArray(iType)%mol(jMol)%x(jAtom) - O_x1
                           ry1 = MolArray(iType)%mol(jMol)%y(jAtom) - O_y1
                           rz1 = MolArray(iType)%mol(jMol)%z(jAtom) - O_z1 
                           r_ij = sqrt(rx1**2 + ry1**2 + rz1**2)
                        end if 
                        if (r_ij < r0) then
                           do kMol = 1, NPART(1) + 1
                              kIndx = molArray(jType)%mol(kMol)%indx
                              ! if (included(kIndx) .eqv. .false.) then
                              !    cycle
                              ! end if
                              if (iType .eq. jType) then
                                 if ((kMol .ne. (NPART(1)+1)) .and. (jMol .ne. (NPART(1)+1)) .and. (iMol .ne. (NPART(1)+1)))  then
                                    cycle
                                 end if
                                 if (jMol > kMol) then 
                                    cycle
                                 end if
                              end if
                              ! finder another O-H that is within r0
                              do kAtom = 2, nAtoms(jType)
                                 if (((jMol == kMol) .and. (jAtom < kAtom)) .or. (jMol < kMol)) then
                                    if (kMol == NPART(1) + 1) then 
                                       rx2 = rosenTrial(nRosen)%x(kAtom) - O_x1
                                       ry2 = rosenTrial(nRosen)%y(kAtom) - O_y1
                                       rz2 = rosenTrial(nRosen)%z(kAtom) - O_z1 
                                       r_ik = sqrt(rx2**2 + ry2**2 + rz2**2)
                                    else
                                       rx2 = MolArray(iType)%mol(kMol)%x(kAtom) - O_x1
                                       ry2 = MolArray(iType)%mol(kMol)%y(kAtom) - O_y1
                                       rz2 = MolArray(iType)%mol(kMol)%z(kAtom) - O_z1 
                                       r_ik = sqrt(rx2**2 + ry2**2 + rz2**2)
                                    end if 
                                    if (r_ik < r0) then
                                       call threebody_interaction(rx1, rx2, ry1, ry2, rz1, rz2, r_ij, r_ik, &
                                                                  B, xi, r0, costheta, threebody)
                                       E_threebody = E_threebody + threebody
                                    end if
                                 end if
                              end do 
                           end do 
                        end if 
                     end do 
                  end do
            end do 
         end do 
      end do

      E_Trial = (E_steric + E_coulomb + E_chargeDipole + E_vanderWaals)
      E_Trial = E_Trial + E_threebody
      print *, "REVERSE:", E_Trial 
      print *, "E_3body", E_threebody 
      print *, "E_2body", (E_steric + E_coulomb + E_chargeDipole + E_vanderWaals)
   end subroutine
!======================================================================================================
   ! pure subroutine Rosen_Vashishta_Molecule_Old(mol_x, mol_y, mol_z, nType, included, E_Trial, iMol)
   subroutine Rosen_Vashishta_Molecule_Old(mol_x, mol_y, mol_z, nType, included, E_Trial, iMol)
      use ForceField, only: nAtoms, r_min_tab, atomArray, atomDataV
      use ForceFieldPara_Vashishta, only: q_tab, H_tab, eta_tab, W_tab, r1s_tab, D_Tab, r4s_tab
      use Coords , only: rosenTrial, MolArray
      use SimParameters, only: nMolTypes, NPART, r1s, rc
      implicit none

      logical, intent(in) :: included(:)
      integer, intent(in) :: nType
      integer, intent(in)  :: iMol
      real(dp), intent(in) :: mol_x(:), mol_y(:), mol_z(:)
      real(dp), intent(out) :: E_Trial

      integer :: iAtom, iiIndx, iIndx, iType, jType, jIndx, iiMol, jMol, kMol, kIndx, jAtom, kAtom, ctr
      integer(kind=atomIntType) :: atmType1, atmType2
      real(dp) :: rx, ry, rz, r, rx1, ry1, rz1, rx2, ry2, rz2
      real(dp) :: rmin_ij
      real(dp) :: E_steric, E_coulomb, E_chargeDipole, E_vanderWaals, E_twobody, E_threebody
      real(dp) :: steric, coulomb, chargeDipole, vanderWaals, twobody, threebody
      real(dp) :: H, q, eta, r4s, D, W ! 2-body parameter
      real(dp) :: B, xi, r0, C, theta  ! 3-body parameters
      real(dp) :: O_x1, O_y1, O_z1
      real(dp) :: r_ij, r_ik, costheta, costheta_ijk
      real(dp) :: rcheta, rcqr1s, rcDr4s, rcW, tmp
      character(len=100) :: fl_name


      E_Trial = 0E0

      steric = 0d0 
      coulomb = 0d0 
      chargeDipole = 0d0 
      vanderWaals = 0d0 
      twobody = 0d0
      threebody = 0d0

      E_steric = 0d0 
      E_coulomb = 0d0 
      E_chargeDipole = 0d0 
      E_vanderWaals = 0d0 
      E_twobody = 0d0 
      E_threebody = 0d0

      print *, "iMol", iMol
      ! 2-body "intramolecular" contribution of the trial molecule
      do iAtom = 1, nAtoms(nType)-1
         atmType1 = atomArray(nType, iAtom)
         do jType = 1, nMolTypes
            do jAtom = iAtom+1, nAtoms(jType)
               atmType2 = atomArray(jType, jAtom)
               q = q_tab(atmType2, atmType1)
               H = H_tab(atmType2, atmType1)
               eta = eta_tab(atmType2, atmType1)
               r4s = r4s_tab(atmType2, atmType1)
               D = D_tab(atmType2, atmType1)
               W = W_tab(atmType2, atmType1)
               rmin_ij = r_min_tab(atmType2, atmType1)
         
               rx = mol_x(iAtom) - mol_x(jAtom)
               ry = mol_y(iAtom) - mol_y(jAtom)
               rz = mol_z(iAtom) - mol_z(jAtom)
               r = rx*rx + ry*ry + rz*rz
               if (r .lt. rmin_ij) then
                  E_Trial = huge(dp)
                  return
               end if
               r = sqrt(r)
               if (r < rc) then
                  call two_body_interaction(H, eta, q, D, W, rc, r1s, r4s, r, &
                                          steric, coulomb, chargeDipole, vanderWaals)
                  E_steric = E_steric + steric
                  E_coulomb = E_coulomb  + coulomb
                  E_chargeDipole = E_chargeDipole + chargeDipole
                  E_vanderWaals = E_vanderWaals + vanderWaals 
               end if 
            end do
         end do 
      end do 
      print *, "INTRAMOLECULAR ENERGY:=", E_steric, E_coulomb, E_chargeDipole,  E_vanderWaals, &
                                       "=", E_steric + E_coulomb + E_chargeDipole + E_vanderWaals

      ! calculate intramolecular energy of each molecule
      do jMol = 1, NPART(1)
         jIndx = molArray(1)%mol(jMol)%indx
         ! This check exclude the intramolecular contribution
         if (included(jIndx) .eqv. .false.) then
            cycle
         end if
         do iAtom = 1, nAtoms(nType)-1
            atmType1 = atomArray(nType, iAtom)
            do jType = 1, nMolTypes
               do jAtom = iAtom+1, nAtoms(jType)
                  atmType2 = atomArray(jType, jAtom)
                  q = q_tab(atmType2, atmType1)
                  H = H_tab(atmType2, atmType1)
                  eta = eta_tab(atmType2, atmType1)
                  r4s = r4s_tab(atmType2, atmType1)
                  D = D_tab(atmType2, atmType1)
                  W = W_tab(atmType2, atmType1)
                  rmin_ij = r_min_tab(atmType2, atmType1)
            
                  rx = MolArray(jType)%mol(jMol)%x(iAtom) - MolArray(jType)%mol(jMol)%x(jAtom)
                  ry = MolArray(jType)%mol(jMol)%y(iAtom) - MolArray(jType)%mol(jMol)%y(jAtom)
                  rz = MolArray(jType)%mol(jMol)%z(iAtom) - MolArray(jType)%mol(jMol)%z(jAtom)
                  r = rx*rx + ry*ry + rz*rz
                  if (r .lt. rmin_ij) then
                     E_Trial = huge(dp)
                     return
                  end if
                  r = sqrt(r)
                  if (r < rc) then
                     call two_body_interaction(H, eta, q, D, W, rc, r1s, r4s, r, &
                                             steric, coulomb, chargeDipole, vanderWaals)
                     E_steric = E_steric + steric
                     E_coulomb = E_coulomb  + coulomb
                     E_chargeDipole = E_chargeDipole + chargeDipole
                     E_vanderWaals = E_vanderWaals + vanderWaals 
                  end if 
               end do
            end do 
         end do 
      end do
      print *, "INTRAMOLECULAR ENERGY:=", E_steric, E_coulomb, E_chargeDipole,  E_vanderWaals, &
                                       "=", E_steric + E_coulomb + E_chargeDipole + E_vanderWaals
      ! 2-body intermolecular interaction 
      do iAtom = 1, nAtoms(nType)
         atmType1 = atomArray(nType, iAtom)
         do jType = 1, nMolTypes
            do jAtom = 1, nAtoms(jType)
               atmType2 = atomArray(jType, jAtom)

               q = q_tab(atmType1, atmType2)
               H = H_tab(atmType1, atmType2)
               eta = eta_tab(atmType1, atmType2)
               r4s = r4s_tab(atmType1, atmType2)
               D = D_tab(atmType1, atmType2)
               W = W_tab(atmType1, atmType2)
               rmin_ij = r_min_tab(atmType2, atmType1) 

               do jMol = 1, NPART(jType)
                  jIndx = molArray(jType)%mol(jMol)%indx
                  ! This check exclude the intramolecular contribution
                  if (included(jIndx) .eqv. .false.) then
                     cycle
                  end if
!             New Energy Calculation

                     rx = mol_x(iAtom) - MolArray(jType)%mol(jMol)%x(jAtom)
                     ry = mol_y(iAtom) - MolArray(jType)%mol(jMol)%y(jAtom)
                     rz = mol_z(iAtom) - MolArray(jType)%mol(jMol)%z(jAtom)
                     r = rx*rx + ry*ry + rz*rz
                     
                     if (r .lt. rmin_ij) then
                        E_Trial = huge(dp)
                        return
                     end if

                     r = sqrt(r)
                     if (r < rc) then
                        call two_body_interaction(H, eta, q, D, W, rc, r1s, r4s, r, &
                                          steric, coulomb, chargeDipole, vanderWaals)
                        E_steric = E_steric + steric
                        E_coulomb = E_coulomb  + coulomb
                        E_chargeDipole = E_chargeDipole + chargeDipole
                        E_vanderWaals = E_vanderWaals + vanderWaals 
                     end if
                  ! end if 
               end do
            end do
         end do
      end do
      print *, "INTRAMOLECULAR ENERGY:=", E_steric, E_coulomb, E_chargeDipole,  E_vanderWaals, &
                                       "=", E_steric + E_coulomb + E_chargeDipole + E_vanderWaals

      
      ! 3-BODY INTERACTIONS 
      B  = atomDataV(1)%B
      xi = atomDataV(1)%xi
      r0 = atomDataV(1)%r0
      C = atomDataV(1)%C 
      costheta = atomDataV(1)%theta
      
      ! Intramolecular 3-body contribution 
      O_x1 = mol_x(1)
      O_y1 = mol_y(1)
      O_z1 = mol_z(1)
      rx1 = mol_x(2) - O_x1
      ry1 = mol_y(2) - O_y1
      rz1 = mol_z(2) - O_z1 
      r_ij = sqrt(rx1**2 + ry1**2 + rz1**2)
      rx2 = mol_x(3) - O_x1
      ry2 = mol_y(3) - O_y1
      rz2 = mol_z(3) - O_z1 
      r_ik = sqrt(rx2**2 + ry2**2 + rz2**2)
      call threebody_interaction(rx1, rx2, ry1, ry2, rz1, rz2, r_ij, r_ik, &
                                 B, xi, r0, costheta, threebody)
      E_threebody = E_threebody + threebody
      print *, "INTRAMOLECULAR 3body", E_threebody


      do jType = nType, nMolTypes
         ! do iiMol = 1, NPART(nType)
         do iAtom = 1, nAtoms(nType)
            ! if (iMol == 0) then 
            O_x1 = mol_x(iAtom)
            O_y1 = mol_y(iAtom)
            O_z1 = mol_z(iAtom)
            ! else
            !    iiIndx = molArray(jType)%mol(iiMol)%indx
            !    if (included(iiIndx) .eqv. .false.) then
            !       cycle
            !    end if
            !    O_x1 = MolArray(nType)%mol(iiMol)%x(1)
            !    O_y1 = MolArray(nType)%mol(iiMol)%y(1)
            !    O_z1 = MolArray(nType)%mol(iiMol)%z(1)
            ! end if
            
            do jMol = 1, NPART(jType)
               ! find an O-H that is within r0
               ! start at 2 since Oxygen=1, note: dont mix up the vashishta forcefield O,H sequence in atom definition
               jIndx = molArray(jType)%mol(jMol)%indx
               if (included(jIndx) .eqv. .false.) then
                  cycle
               end if
               if (iAtom > 1) then ! if iAtom is not an oxygen 
                  rx1 = MolArray(jType)%mol(jMol)%x(1) - O_x1
                  ry1 = MolArray(jType)%mol(jMol)%y(1) - O_y1
                  rz1 = MolArray(jType)%mol(jMol)%z(1) - O_z1 
                  r_ij = sqrt(rx1**2 + ry1**2 + rz1**2)
                  print *, r_ij, r0
                  if (r_ij < r0) then
                     do kMol = 1, NPART(jType)
                        kIndx = molArray(jType)%mol(kMol)%indx
                        if (included(kIndx) .eqv. .false.) then
                           cycle
                        end if
                        do kAtom = 2, nAtoms(jType)
                           ! if (((jMol == kMol) .and. (jAtom < kAtom)) .or. (jMol < kMol)) then
                           !    if (iiMol == iMol) then
                           !       rx2 = mol_x(kAtom) - O_x1
                           !       ry2 = mol_y(kAtom) - O_y1
                           !       rz2 = mol_z(kAtom) - O_z1 
                           !    else 
                           rx2 = MolArray(nType)%mol(kMol)%x(kAtom) - O_x1
                           ry2 = MolArray(nType)%mol(kMol)%y(kAtom) - O_y1
                           rz2 = MolArray(nType)%mol(kMol)%z(kAtom) - O_z1 
                           ! end if
                           r_ik = sqrt(rx2**2 + ry2**2 + rz2**2)
                           if (r_ik < r0) then
                              print *,"==>",iAtom, jAtom, jMol, kAtom, kMol
                              call threebody_interaction(rx1, rx2, ry1, ry2, rz1, rz2, r_ij, r_ik, &
                                                   B, xi, r0, costheta, threebody)
                              E_threebody = E_threebody + threebody
                              ! print *, E_threebody
                           end if 
                           ! end if
                        end do 
                     end do
                  end if 
               else
                  do jAtom = 2, nAtoms(jType)  
                     ! if (iiMol == iMol) then
                     !    rx1 = mol_x(jAtom) - O_x1
                     !    ry1 = mol_y(jAtom) - O_y1
                     !    rz1 = mol_z(jAtom) - O_z1 
                     ! else 
                     rx1 = MolArray(jType)%mol(jMol)%x(jAtom) - O_x1
                     ry1 = MolArray(jType)%mol(jMol)%y(jAtom) - O_y1
                     rz1 = MolArray(jType)%mol(jMol)%z(jAtom) - O_z1 
                     ! end if
                     r_ij = sqrt(rx1**2 + ry1**2 + rz1**2)
                     ! print *, "r_ij", r_ij !, r0, B, C, xi, costheta 
                     if (r_ij < r0) then
                        do kMol = 1, NPART(jType)
                           kIndx = molArray(jType)%mol(kMol)%indx
                           if (included(kIndx) .eqv. .false.) then
                              cycle
                           end if
                           ! if (nType .eq. jType) then
                           !    if ((kMol .ne. iMol) .and. (jMol .ne. iMol) .and. (iiMol .ne. iMol))  then
                           !       cycle
                           !    end if
                           !    if (jMol > kMol) then
                           !       cycle
                           !    end if
                           ! end if
                           ! find another O-H that is within r0
                           do kAtom = 2, nAtoms(jType)
                              if (((jMol == kMol) .and. (jAtom < kAtom)) .or. (jMol < kMol)) then
                                 ! if (iiMol == iMol) then
                                 !    rx2 = mol_x(kAtom) - O_x1
                                 !    ry2 = mol_y(kAtom) - O_y1
                                 !    rz2 = mol_z(kAtom) - O_z1 
                                 ! else 
                                 rx2 = MolArray(nType)%mol(kMol)%x(kAtom) - O_x1
                                 ry2 = MolArray(nType)%mol(kMol)%y(kAtom) - O_y1
                                 rz2 = MolArray(nType)%mol(kMol)%z(kAtom) - O_z1 
                                 ! end if
                                 r_ik = sqrt(rx2**2 + ry2**2 + rz2**2)
                                 
                                 if (r_ik < r0) then
                                    print *,"==>",iAtom, jAtom, jMol, kAtom, kMol
                                    call threebody_interaction(rx1, rx2, ry1, ry2, rz1, rz2, r_ij, r_ik, &
                                                         B, xi, r0, costheta, threebody)
                                    E_threebody = E_threebody + threebody
                                    ! print *, E_threebody
                                 end if 
                              end if
                           end do 
                        end do 
                     end if 
                  end do 
               end if 
            end do
         end do 
      end do 
      ! print *, "RosenTrial 2body", E_steric + E_coulomb + E_chargeDipole + E_vanderWaals
      ! print *, "RosenTrial 3body", E_threebody
      E_Trial = (E_steric + E_coulomb + E_chargeDipole + E_vanderWaals)
      E_Trial = E_Trial + E_threebody 
      ! print *, "REVERSE:", E_Trial 
      ! print *, "E_3body", E_threebody 
      ! print *, "E_2body ", (E_steric + E_coulomb + E_chargeDipole + E_vanderWaals)
      ! print *, "E_2body -- steric", E_steric 
      ! print *, "E_2body -- coulomb", E_coulomb 
      ! print *, "E_2body -- charge Dipole",  E_chargeDipole
      ! print *, "E_2body -- vander Waals", E_vanderWaals
   end subroutine
end module





     
      ! 2-body "intramolecular" contribution of the trial molecule
      ! do iAtom = 1, nAtoms(nType)-1
      !    atmType1 = atomArray(nType, iAtom)
      !    do jType = 1, nMolTypes
      !       do jAtom = iAtom+1, nAtoms(jType)
      !          atmType2 = atomArray(jType, jAtom)
      !          q = q_tab(atmType2, atmType1)
      !          H = H_tab(atmType2, atmType1)
      !          eta = eta_tab(atmType2, atmType1)
      !          r4s = r4s_tab(atmType2, atmType1)
      !          D = D_tab(atmType2, atmType1)
      !          W = W_tab(atmType2, atmType1)
      !          rmin_ij = r_min_tab(atmType2, atmType1)
               
      !          rx = mol_x(iAtom) - mol_x(jAtom)
      !          ry = mol_y(iAtom) - mol_y(jAtom)
      !          rz = mol_z(iAtom) - mol_z(jAtom)
      !          r = rx*rx + ry*ry + rz*rz
      !          if (r .lt. rmin_ij) then
      !             E_Trial = huge(dp)
      !             return
      !          end if
      !          r = sqrt(r)
      !          if (r < rc) then
      !             call two_body_interaction(H, eta, q, D, W, rc, r1s, r4s, r, &
      !                                     steric, coulomb, chargeDipole, vanderWaals)
      !             E_steric = E_steric + steric
      !             E_coulomb = E_coulomb  + coulomb
      !             E_chargeDipole = E_chargeDipole + chargeDipole
      !             E_vanderWaals = E_vanderWaals + vanderWaals 
      !          end if 
      !       end do
      !    end do 
      ! end do 
      



      ! ========= INTRAMOLECULAR 3-BODY contribution =========================
      ! O_x1 = mol_x(1)
      ! O_y1 = mol_y(1)
      ! O_z1 = mol_z(1)
      ! rx1 = mol_x(2) - O_x1
      ! ry1 = mol_y(2) - O_y1
      ! rz1 = mol_z(2) - O_z1 
      ! r_ij = sqrt(rx1**2 + ry1**2 + rz1**2)
      ! rx2 = mol_x(3) - O_x1
      ! ry2 = mol_y(3) - O_y1
      ! rz2 = mol_z(3) - O_z1 
      ! r_ik = sqrt(rx2**2 + ry2**2 + rz2**2)
      ! print *, r_ik, r_ij
      ! if (r_ik < r0 .and. r_ij < r0) then
      !    call threebody_interaction(rx1, rx2, ry1, ry2, rz1, rz2, r_ij, r_ik, &
      !                               B, xi, r0, costheta, threebody)
      !    E_threebody = E_threebody + threebody
      ! end if 
      ! ========= INTRAMOLECULAR 3-BODY contribution =========================


      if ((iMol .ne. NPART(1) + 1) .and. jMol < NPART(1) + 1) then 
                     if ((jMol == iMol .and. iAtom  < jAtom) .or. iMol .ne. jMol) then
                        rx = mol_x(iAtom) - MolArray(jType)%mol(jMol)%x(jAtom)
                        ry = mol_y(iAtom) - MolArray(jType)%mol(jMol)%y(jAtom)
                        rz = mol_z(iAtom) - MolArray(jType)%mol(jMol)%z(jAtom)
                     end if 
                  else if ((iMol == NPART(1) + 1))
                     else if ((iMol == NPART(1) + 1) .and. iAtom  < jAtom) then 
                        rx = mol_x(iAtom) - mol_x(jAtom)
                        ry = mol_y(iAtom) - mol_y(jAtom)
                        rz = mol_z(iAtom) - mol_z(jAtom)
                     else if ((iMol == NPART(1) + 1) .and. jMol < NPART(1) + 1 .and. iAtom  < jAtom) then 
                  end if






       !==========================================
         ! TEST ITS REVERSIBILITY!! 
         !==========================================
            ! E_Inter = 0d0
            ! E_Intra = 0d0
            ! call Create_NeiETable(nType)
            ! call EBias_Remove_ChooseTarget(nTarget, nTargType, nTargMol, ProbTargOut, nTargMol)
            ! call EBias_Remove_ChooseNeighbor(nTarget, nType, nSel, ProbSel, NPART(1))
            ! nSel = NPART(1)
            ! nMol = subIndxList(nSel)
            ! call Simple_RosenConfigGen_Reverse_Vashishta(nType, nMol, nTarget, nTargType, rosenRatio)
            ! call SwapOut_ECalc(E_Inter, E_Intra, nType, nMol, dETable)
            ! call EBias_Remove_ReverseProbTarget(nTarget, nSel, nType, dETable, ProbTargIn)

            ! genProbRatio = (ProbTargIn*rosenRatio)/(ProbTargOut*ProbSel*avbmc_vol*gas_dens(nType))
            ! A21 = genProbRatio
            ! Boltzterm = exp(-beta*E_Inter)
            ! P21 = log(genProbRatio) - beta*E_Inter
            ! print *, "genProbRatio", genProbRatio, Boltzterm
            ! print *, "AVBMC_CheckReversibility: deletion", genProbRatio*Boltzterm, P21
            ! print *, "DELETION E_INTER", E_Inter

            ! ! molArray(nType)%mol(nMol)%x(1:nAtoms(nType)) = molArray(nType)%mol(NPART(nType))%x(1:nAtoms(nType))
            ! ! molArray(nType)%mol(nMol)%y(1:nAtoms(nType)) = molArray(nType)%mol(NPART(nType))%y(1:nAtoms(nType))
            ! ! molArray(nType)%mol(nMol)%z(1:nAtoms(nType)) = molArray(nType)%mol(NPART(nType))%z(1:nAtoms(nType))
            ! ! call NeighborUpdate_Delete(nIndx, molArray(nType)%mol(NPART(nType))%indx)
            ! isActive(molArray(nType)%mol(NPART(nType))%indx) = .false.
            ! E_T = E_T + E_Inter !+ E_Intra
            ! ETable = ETable - dETable
            ! ETable(nIndx) = ETable(molArray(nType)%mol(NPART(nType))%indx)
            ! ETable(molArray(nType)%mol(NPART(nType))%indx) = 0d0
            ! NPART(nType) = NPART(nType) - 1
            ! NTotal = NTotal - 1
            ! call Update_SubEnergies_Vashishta
         
            ! if( ((A12 == 0E0_dp) .or. (A21 == 0E0_dp)) .and. (A21 /= A12) ) then
            !    print *, "Zero Probability observed in move that should be reverisbile"
            !    print *, "This implies a calculation error in the detailed balance."
            !    print *, "A12, A21:", A12, A21
            !    error stop
            ! endif

            ! if(P12 > 0E0_dp) P12 = 0E0_dp
            ! if(P21 > 0E0_dp) P21 = 0E0_dp

            ! detailBalance = log(genProbRatio) + P12 - P21 - beta*E_Inter
            ! if (abs(detailBalance) > 1e-6_dp) then 
            !    print *, "Violation of Detailed Balance Detected!" 
            ! end if 
            ! print *, "E_inter", E_inter 
            ! print *, "A12, A21, A12*A21:", A12, A21, A12*A21 
            ! print *, "ln(P12), ln(P21):", P12, P21
            ! print *, "Detailed balance", detailBalance, log(genProbRatio), beta*E_Inter
         !==========================================
         ! TEST ITS REVERSIBILITY!! 
         !==========================================






!===================================================================================
!    subroutine AVBMC_CheckReversibility(E_T, arrayMax, accept)
!       use AVBMC_CBMC
!       use CBMC_Utility
!       use CBMC_Variables
!       use Constants
!       use Coords
!       use SimParameters
!       use ForceField
!       use EnergyPointers, only: SwapIn_ECalc, SwapOut_ECalc, Update_SubEnergies_Vashishta
!       use ForceField
!       use IndexingFunctions
!       use EnergyCriteria
!       use DistanceCriteria
!       use EnergyTables
!       use NeighborTable
!       use SwapBoundary
!       implicit none

!       integer, intent(in) :: arrayMax
!       real(dp), intent(inout) :: E_T
!       logical, intent(out) :: accept 
!       logical :: rejMove
!       logical :: isIncluded(1:arrayMax)
!       integer :: i, nTargType, nTargMol, nTargIndx, nTarget
!       integer :: nType, nIndx,  nSel, nMol   
!       real(dp) :: grnd
!       real(dp) :: genProbRatio, rosenRatio
!       real(dp) :: E_Inter, E_Intra, biasDiff
!       real(dp) :: PairList(1:arrayMax)
!       real(dp) :: dETable(1:arrayMax)
!       real(dp) :: newNeiETable(1:arrayMax)
!       real(dp) :: ProbTarg_In, ProbTarg_Out, ProbSel_Out  ! insertion
!       real(dp) :: ProbTargOut, ProbSel, ProbTargIn     ! deletion
!       real(dp) :: Boltzterm, P12, P21, A12, A21, detailBalance
!       nType = 1
!       !==========================================
!       ! GET THE INSERTION PROBABILITY 
!       !==========================================
! !      Choose the type of molecule to be inserted
!       E_Inter = 0d0
!       E_Intra = 0d0
!       print *, NPART(1), maxMol
!          call EBias_Insert_ChooseTarget(nType, nTarget, nTargType, nTargMol, ProbTarg_In)
!          nTargIndx = MolArray(nTargType)%mol(nTargMol)%indx
!          nIndx = MolArray(nType)%mol(NPART(nType) + 1)%indx
!          call Rosen_CreateSubset(nTarget, isIncluded)
!          call Simple_RosenConfigGen_Vashishta(nType, nIndx, nTarget, nTargType, isIncluded, rosenRatio, rejMove)
!          call SwapIn_ECalc(E_Inter, E_Intra, PairList, dETable, rejMove)
!          call Insert_NewNeiETable_Distance_V2(nType, PairList, dETable, newNeiETable) ! update the table
!          call EBias_Insert_ReverseProbTarget(nTarget, nType, newNeiETable, ProbTarg_Out) 
!          call EBias_Insert_ReverseProbSel(nTarget, nType, dETable, ProbSel_Out)
!          genProbRatio = (ProbTarg_Out*ProbSel_Out*avbmc_vol*gas_dens(nType))/(ProbTarg_In*rosenRatio)
!          A12 = genProbRatio
!          Boltzterm = exp(-beta*E_Inter)
!          P12 = log(genProbRatio) - beta*E_Inter
!          print *, "genProbRatio", genProbRatio, Boltzterm
!          print *, "AVBMC_CheckReversibility: insertion", genProbRatio*Boltzterm, P12
!          print *, "INSERTION E_INTER", E_Inter

!          do i = 1, nAtoms(nType)
!             molArray(nType)%mol(NPART(nType) + 1)%x(i) = newMol%x(i)
!             molArray(nType)%mol(NPART(nType) + 1)%y(i) = newMol%y(i)
!             molArray(nType)%mol(NPART(nType) + 1)%z(i) = newMol%z(i)
!          end do
!          isActive(nIndx) = .true.
!          call NeighborUpdate_SwapIn_Distance(PairList, nType)
!          NTotal = NTotal + 1
!          NPART(1) = NPART(1) + 1
!          E_T = E_T + E_Inter
!          ETable = ETable + dETable
!          call Update_SubEnergies_Vashishta

!          !==========================================
!          ! GET THE DELETION PROBABILITY
!          !==========================================
!          E_Inter = 0d0
!          E_Intra = 0d0
!          call Create_NeiETable(nType)
!          call EBias_Remove_ChooseTarget(nTarget, nTargType, nTargMol, ProbTargOut, nTargMol)
!          call EBias_Remove_ChooseNeighbor(nTarget, nType, nSel, ProbSel, NPART(1))
!          nSel = NPART(1)
!          nMol = subIndxList(nSel)
!          call Simple_RosenConfigGen_Reverse_Vashishta(nType, nMol, nTarget, nTargType, rosenRatio)
!          call SwapOut_ECalc(E_Inter, E_Intra, nType, nMol, dETable)
!          call EBias_Remove_ReverseProbTarget(nTarget, nSel, nType, dETable, ProbTargIn)

!          genProbRatio = (ProbTargIn*rosenRatio)/(ProbTargOut*ProbSel*avbmc_vol*gas_dens(nType))
!          A21 = genProbRatio
!          Boltzterm = exp(-beta*E_Inter)
!          P21 = log(genProbRatio) - beta*E_Inter
!          print *, "genProbRatio", genProbRatio, Boltzterm
!          print *, "AVBMC_CheckReversibility: deletion", genProbRatio*Boltzterm, P21
!          print *, "DELETION E_INTER", E_Inter

!          molArray(nType)%mol(nMol)%x(1:nAtoms(nType)) = molArray(nType)%mol(NPART(nType))%x(1:nAtoms(nType))
!          molArray(nType)%mol(nMol)%y(1:nAtoms(nType)) = molArray(nType)%mol(NPART(nType))%y(1:nAtoms(nType))
!          molArray(nType)%mol(nMol)%z(1:nAtoms(nType)) = molArray(nType)%mol(NPART(nType))%z(1:nAtoms(nType))
!          call NeighborUpdate_Delete(nIndx, molArray(nType)%mol(NPART(nType))%indx)
!          isActive(molArray(nType)%mol(NPART(nType))%indx) = .false.
!          E_T = E_T + E_Inter !+ E_Intra
!          ETable = ETable - dETable
!          ETable(nIndx) = ETable(molArray(nType)%mol(NPART(nType))%indx)
!          ETable(molArray(nType)%mol(NPART(nType))%indx) = 0d0
!          NPART(nType) = NPART(nType) - 1
!          NTotal = NTotal - 1
!          call Update_SubEnergies_Vashishta
      
!       if( ((A12 == 0E0_dp) .or. (A21 == 0E0_dp)) .and. (A21 /= A12) ) then
!          print *, "Zero Probability observed in move that should be reverisbile"
!          print *, "This implies a calculation error in the detailed balance."
!          print *, "A12, A21:", A12, A21
!          error stop
!       endif

!       if(P12 > 0E0_dp) P12 = 0E0_dp
!       if(P21 > 0E0_dp) P21 = 0E0_dp

!       detailBalance = log(genProbRatio) + P12 - P21 - beta*E_Inter
!       if (abs(detailBalance) > 1e-6_dp) then 
!          print *, "Violation of Detailed Balance Detected!" 
!       end if 
!       print *, "E_inter", E_inter 
!       print *, "A12, A21, A12*A21:", A12, A21, A12*A21 
!       print *, "ln(P12), ln(P21):", P12, P21
!       print *, "Detailed balance", detailBalance, log(genProbRatio), beta*E_Inter
!    end subroutine 






         ! endif

         ! if(P12 > 0E0_dp) P12 = 0E0_dp
         ! if(P21 > 0E0_dp) P21 = 0E0_dp

         ! detailBalance = log(A21) + P12 - P21 - beta*E_deletion
         ! ! detailBalance = log(A12) + P21 - P12 + beta*E_insertion
         ! if (abs(detailBalance) > 1e-6_dp) then 
         !    print *, "E_inter", E_inter 
         !    print *, "A12, A21, A12*A21:", A12, A21, A12*A21 
         !    print *, "ln(P12), ln(P21):", P12, P21
         !    print *, "Detailed balance", detailBalance, log(genProbRatio), beta*E_Inter
         !    stop "Violation of Detailed Balance Detected!" 
         ! end if 
         ! ==========================================
         ! TEST ITS REVERSIBILITY!! 
         ! ==========================================





!       nType = 1
!       !==========================================
!       ! GET THE INSERTION PROBABILITY 
!       !==========================================
! !      Choose the type of molecule to be inserted
!       E_Inter = 0d0
!       E_Intra = 0d0
!       print *, "CR start here=>",NPART(1), maxMol, regrowType(nType)
!          NDiff = 0
!          NDiff(nType) = 1
!          rejMove = boundaryFunction(NPART, NDiff)
!          call EBias_Insert_ChooseTarget(nType, nTarget, nTargType, nTargMol, ProbTarg_In)
!          nTargIndx = MolArray(nTargType)%mol(nTargMol)%indx
!          nIndx = MolArray(nType)%mol(NPART(nType) + 1)%indx
!          call Rosen_CreateSubset(nTarget, isIncluded)
!          select case (regrowType(nType))
!          case (0)
!             call Ridgid_RosenConfigGen(nType, nIndx, nTarget, nTargType, isIncluded, rosenRatio, rejMove)
!          case (1)
!             call Simple_RosenConfigGen(nType, nIndx, nTarget, nTargType, isIncluded, rosenRatio, rejMove)
!          case (2)
!             call StraightChain_RosenConfigGen(nType, nIndx, nTarget, nTargType, isIncluded, rosenRatio, rejMove)
!          case (3)
!             call BranchedMol_RosenConfigGen(nType, nIndx, nTarget, nTargType, isIncluded, rosenRatio, rejMove)
!          case (4)
!             call Simple_RosenConfigGen_Vashishta(nType, nIndx, nTarget, nTargType, isIncluded, rosenRatio, rejMove)
!          case (5)
!             call Ridgid_RosenConfigGen_Vashishta(nType, nIndx, nTarget, nTargType, isIncluded, rosenRatio, rejMove)
!          case default
!             write (*, *) "Error! EBias can not regrow a molecule of regrow type:", regrowType(nType)
!             stop
!          end select
!          call SwapIn_ECalc(E_Inter, E_Intra, PairList, dETable, rejMove)
!          call Insert_NewNeiETable_Distance_V2(nType, PairList, dETable, newNeiETable) ! update the table
!          call EBias_Insert_ReverseProbTarget(nTarget, nType, newNeiETable, ProbTarg_Out) 
!          call EBias_Insert_ReverseProbSel(nTarget, nType, dETable, ProbSel_Out)
!          genProbRatio = (ProbTarg_Out*ProbSel_Out*avbmc_vol*gas_dens(nType))/(ProbTarg_In*rosenRatio)
!          print *, "Insertion Probabilities: ",ProbTarg_In, rosenRatio, ProbTarg_Out, ProbSel_Out, avbmc_vol, gas_dens(nType)
!          A12 = genProbRatio
!          Boltzterm = exp(-beta*E_Inter)
!          P12 = log(genProbRatio) - beta*E_Inter
!          print *, "genProbRatio", genProbRatio, Boltzterm
!          print *, "AVBMC_CheckReversibility: insertion", genProbRatio*Boltzterm, P12
!          print *, "INSERTION E_INTER", E_Inter
!          print *, "TARGET during the insertion", nTarget
!          do i = 1, nAtoms(nType)
!             molArray(nType)%mol(NPART(nType) + 1)%x(i) = newMol%x(i)
!             molArray(nType)%mol(NPART(nType) + 1)%y(i) = newMol%y(i)
!             molArray(nType)%mol(NPART(nType) + 1)%z(i) = newMol%z(i)
!          end do
!          isActive(nIndx) = .true.
!          call NeighborUpdate_SwapIn_Distance(PairList, nType)
!          NTotal = NTotal + 1
!          NPART(1) = NPART(1) + 1
!          E_T = E_T + E_Inter
!          ETable = ETable + dETable
!          select case (trim(adjustl(ForceFieldName)))
!          case('Vashishta')
!             call Update_SubEnergies_Vashishta
!          case default
!             print *, "this is called!"
!             call Update_SubEnergies
!          end select

!       !    print *, ETable
!       ! !    !==========================================
!       ! !    ! GET THE DELETION PROBABILITY
!       ! !    !==========================================
!          ! E_Inter = 0d0
!          ! E_Intra = 0d0
!          call Create_NeiETable(nType)
!          call EBias_Remove_ChooseTarget(nTarget, nTargType, nTargMol, ProbTargOut, nTargMol)
!          call EBias_Remove_ChooseNeighbor(nTarget, nType, nSel, ProbSel, NPART(1))
!          nSel = NPART(1)
!          nMol = subIndxList(nSel)
!          select case (regrowType(nType))
!          case (0)
!             call Ridgid_RosenConfigGen_Reverse(nType, nMol, nTarget, nTargType, rosenRatio)
!          case (1)
!             call Simple_RosenConfigGen_Reverse(nType, nMol, nTarget, nTargType, rosenRatio)
!          case (2)
!             call StraightChain_RosenConfigGen_Reverse(nType, nMol, nTarget, nTargType, rosenRatio)
!          case (3)
!             call BranchedMol_RosenConfigGen_Reverse(nType, nMol, nTarget, nTargType, rosenRatio)
!          case (4)
!             call Simple_RosenConfigGen_Reverse_Vashishta(nType, nMol, nTarget, nTargType, rosenRatio)
!          case (5)
!             call Ridgid_RosenConfigGen_Reverse_Vashishta(nType, nMol, nTarget, nTargType, rosenRatio)
!          case default
!             write (*, *) "Error! EBias can not regrow a molecule of regrow type:", nType
!             stop
!          end select
!          call SwapOut_ECalc(E_Inter, E_Intra, nType, nMol, dETable)
!          call EBias_Remove_ReverseProbTarget(nTarget, nSel, nType, dETable, ProbTargIn)

!          genProbRatio = (ProbTargIn*rosenRatio)/(ProbTargOut*ProbSel*avbmc_vol*gas_dens(nType))
!          print *, "Deletion Probabilities: ",ProbTargIn, rosenRatio, ProbTargOut, ProbSel, avbmc_vol, gas_dens(nType)
!          A21 = genProbRatio
!          Boltzterm = exp(-beta*E_Inter)
!          P21 = log(genProbRatio) - beta*E_Inter
!          print *, "genProbRatio", genProbRatio, Boltzterm
!          print *, "AVBMC_CheckReversibility: deletion", genProbRatio*Boltzterm, P21
!          print *, "DELETION E_INTER", E_Inter, E_Intra

!          molArray(nType)%mol(nMol)%x(1:nAtoms(nType)) = molArray(nType)%mol(NPART(nType))%x(1:nAtoms(nType))
!          molArray(nType)%mol(nMol)%y(1:nAtoms(nType)) = molArray(nType)%mol(NPART(nType))%y(1:nAtoms(nType))
!          molArray(nType)%mol(nMol)%z(1:nAtoms(nType)) = molArray(nType)%mol(NPART(nType))%z(1:nAtoms(nType))
!          call NeighborUpdate_Delete(nIndx, molArray(nType)%mol(NPART(nType))%indx)
!          isActive(molArray(nType)%mol(NPART(nType))%indx) = .false.
!          E_T = E_T + E_Inter !+ E_Intra
!          ETable = ETable - dETable
!          ETable(nIndx) = ETable(molArray(nType)%mol(NPART(nType))%indx)
!          ETable(molArray(nType)%mol(NPART(nType))%indx) = 0d0
!          NPART(nType) = NPART(nType) - 1
!          NTotal = NTotal - 1
!          select case (trim(adjustl(ForceFieldName)))
!          case('Vashishta')
!             call Update_SubEnergies_Vashishta
!          case default
!             print *, "this is called!"
!             call Update_SubEnergies
!          end select
      
!       if( ((A12 == 0E0_dp) .or. (A21 == 0E0_dp)) .and. (A21 /= A12) ) then
!          print *, "Zero Probability observed in move that should be reverisbile"
!          print *, "This implies a calculation error in the detailed balance."
!          print *, "A12, A21:", A12, A21
!          error stop
!       endif

!       if(P12 > 0E0_dp) P12 = 0E0_dp
!       if(P21 > 0E0_dp) P21 = 0E0_dp

!       detailBalance = log(A21) + P12 - P21 - beta*E_Inter
!       print *, "E_inter", E_inter 
!       print *, "A12, A21, A12*A21:", A12, A21, A12*A21 
!       print *, "ln(P12), ln(P21):", P12, P21
!       print *, "Detailed balance", detailBalance, log(genProbRatio), beta*E_Inter
!       ! detailBalance = log(A12) + P21 - P12 + beta*E_insertion
!       if (abs(detailBalance) > 1e-6_dp) then 
!          print *, "E_inter", E_inter 
!          print *, "A12, A21, A12*A21:", A12, A21, A12*A21 
!          print *, "ln(P12), ln(P21):", P12, P21
!          print *, "Detailed balance", detailBalance, log(genProbRatio), beta*E_Inter
!          stop "Violation of Detailed Balance Detected!" 
!       end if 