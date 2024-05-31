module Rosenbluth_Functions_Vashishta_Q
   use VarPrecision
contains
!======================================================================================================
   pure subroutine two_body_interaction(H, eta, q, D, W, rc, r1s, r4s, r, steric, coulomb, chargeDipole, vanderWaals)
      use SimParameters, only : outputEConv 
      implicit none
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
      implicit none
      real(dp), intent(in) :: rx1, rx2, ry1, ry2, rz1, rz2, r_ij, r_ik, B, xi, r0, costheta
      real(dp), intent(out) :: threebody
      real(dp) :: costheta_ijk
      costheta_ijk = ((rx1*rx2 + ry1*ry2 + rz1*rz2)/(r_ij*r_ik)) - costheta
      threebody = B*exp((xi/(r_ij - r0)) + (xi/(r_ik - r0)))*(costheta_ijk**2)
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
      use SimParameters, only : outputEConv 
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
               if (r <= rc) then
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
                  if (r <= rc) then
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
                  else 
                     iIndx = molArray(jType)%mol(iMol)%indx
                     O_x1 = MolArray(iType)%mol(iMol)%x(1)
                     O_y1 = MolArray(iType)%mol(iMol)%y(1)
                     O_z1 = MolArray(iType)%mol(iMol)%z(1)
                  end if 
                  do jMol = 1, NPART(1) + 1 
                     jIndx = molArray(jType)%mol(jMol)%indx
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
   end subroutine
!======================================================================================================
   ! pure subroutine Rosen_Vashishta_Molecule_Old(mol_x, mol_y, mol_z, nType, included, E_Trial, iMol)
   subroutine Rosen_Vashishta_Molecule_Old(mol_x, mol_y, mol_z, nType, included, E_Trial, iMol, X)
      use ForceField, only: nAtoms, r_min_tab, atomArray, atomDataV
      use ForceFieldPara_Vashishta, only: q_tab, H_tab, eta_tab, W_tab, r1s_tab, D_Tab, r4s_tab
      use Coords , only: rosenTrial, MolArray
      use SimParameters, only: nMolTypes, NPART, r1s, rc
      implicit none

      logical, intent(in) :: included(:)
      integer, intent(in) :: nType
      integer, intent(in)  :: iMol, X
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

               do jMol = 1, NPART(jType) + 1
                  if (iMol == X .and. jMol <= NPART(jType)) then 
                     if (iMol == jMol .and. iAtom == jAtom) then ! .or. iAtom > jAtom) then
                        cycle
                     else if (iMol==jMol .and. iAtom > jAtom) then
                        cycle
                     end if
                     if ((jMol == iMol .and. iAtom  < jAtom) .or. iMol .ne. jMol) then
                        rx = mol_x(iAtom) - MolArray(jType)%mol(jMol)%x(jAtom)
                        ry = mol_y(iAtom) - MolArray(jType)%mol(jMol)%y(jAtom)
                        rz = mol_z(iAtom) - MolArray(jType)%mol(jMol)%z(jAtom)
                     end if
                     r = rx*rx + ry*ry + rz*rz

                     if (r .lt. rmin_ij) then
                        E_Trial = huge(dp)
                        return
                     end if

                     r = sqrt(r)
                     if (r <= rc) then
                        call two_body_interaction(H, eta, q, D, W, rc, r1s, r4s, r, &
                                          steric, coulomb, chargeDipole, vanderWaals)
                        E_steric = E_steric + steric
                        E_coulomb = E_coulomb  + coulomb
                        E_chargeDipole = E_chargeDipole + chargeDipole
                        E_vanderWaals = E_vanderWaals + vanderWaals 
                     end if
                  else if (iMol .ne. X) then 
                     if (jMol == X .and. iAtom == jAtom) then 
                        cycle
                     else if (jMol == X .and. iAtom > jAtom) then
                        cycle
                     end if

                     if (jMol == iMol) then  ! skip the selected molecule in the calculation
                        cycle  
                     end if 

                     if (jMol <= NPART(jType)) then 
                        rx = mol_x(iAtom) - MolArray(jType)%mol(jMol)%x(jAtom)
                        ry = mol_y(iAtom) - MolArray(jType)%mol(jMol)%y(jAtom)
                        rz = mol_z(iAtom) - MolArray(jType)%mol(jMol)%z(jAtom)
                     else if (jMol == X .and. iAtom  < jAtom) then 
                        rx = mol_x(iAtom) - mol_x(jAtom)
                        ry = mol_y(iAtom) - mol_y(jAtom)
                        rz = mol_z(iAtom) - mol_z(jAtom)
                     end if
                     r = rx*rx + ry*ry + rz*rz

                     if (r .lt. rmin_ij) then
                        E_Trial = huge(dp)
                        return
                     end if

                     r = sqrt(r)
                     if (r <= rc) then
                        
                        call two_body_interaction(H, eta, q, D, W, rc, r1s, r4s, r, &
                                          steric, coulomb, chargeDipole, vanderWaals)
                        E_steric = E_steric + steric
                        E_coulomb = E_coulomb  + coulomb
                        E_chargeDipole = E_chargeDipole + chargeDipole
                        E_vanderWaals = E_vanderWaals + vanderWaals 
                     end if
                  else 
                     cycle
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

      do jType = nType, nMolTypes
         do iiMol = 1, NPART(nType) + 1
            if (iiMol == NPART(1) + 1 .and. iMol .ne. NPART(1) + 1 ) then 
               cycle 
            end if 
            if ((iiMol == NPART(1) + 1) .and. (iMol == NPART(1) + 1)) then 
               O_x1 = mol_x(1)
               O_y1 = mol_y(1)
               O_z1 = mol_z(1)
            else
               O_x1 = MolArray(nType)%mol(iiMol)%x(1)
               O_y1 = MolArray(nType)%mol(iiMol)%y(1)
               O_z1 = MolArray(nType)%mol(iiMol)%z(1)
            end if
            
            do jMol = 1, NPART(jType) + 1
               ! find an O-H that is within r0
               ! start at 2 since Oxygen=1, note: dont mix up the vashishta forcefield O,H sequence in atom definition
               if (jMol == NPART(1) + 1 .and. iMol .ne. NPART(1) + 1 ) then 
                  cycle 
               end if 
               do jAtom = 2, nAtoms(jType)  
                  if ((jMol == NPART(1)+1) .and. (iMol == NPART(1) + 1)) then
                     rx1 = mol_x(jAtom) - O_x1
                     ry1 = mol_y(jAtom) - O_y1
                     rz1 = mol_z(jAtom) - O_z1 
                  else 
                     rx1 = MolArray(jType)%mol(jMol)%x(jAtom) - O_x1
                     ry1 = MolArray(jType)%mol(jMol)%y(jAtom) - O_y1
                     rz1 = MolArray(jType)%mol(jMol)%z(jAtom) - O_z1 
                  end if
                  r_ij = sqrt(rx1**2 + ry1**2 + rz1**2)
                  if (r_ij < r0) then
                     do kMol = 1, NPART(jType) + 1
                        if (kMol == NPART(1) + 1 .and. iMol .ne. NPART(1) + 1 ) then 
                           cycle 
                        end if 
                        if (nType .eq. jType) then
                           if ((kMol .ne. iMol) .and. (jMol .ne. iMol) .and. (iiMol .ne. iMol))  then
                              cycle
                           end if
                           if (jMol > kMol) then
                              cycle
                           end if
                        end if
                        ! find another O-H that is within r0
                        do kAtom = 2, nAtoms(jType)
                           if (((jMol == kMol) .and. (jAtom < kAtom)) .or. (jMol < kMol)) then
                              if ((kMol ==  NPART(1) + 1) .and. (iMol == NPART(1) + 1)) then
                                 rx2 = mol_x(kAtom) - O_x1
                                 ry2 = mol_y(kAtom) - O_y1
                                 rz2 = mol_z(kAtom) - O_z1 
                              else 
                                 rx2 = MolArray(nType)%mol(kMol)%x(kAtom) - O_x1
                                 ry2 = MolArray(nType)%mol(kMol)%y(kAtom) - O_y1
                                 rz2 = MolArray(nType)%mol(kMol)%z(kAtom) - O_z1 
                              end if
                              r_ik = sqrt(rx2**2 + ry2**2 + rz2**2)
                              
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

      E_Trial = (E_steric + E_coulomb + E_chargeDipole + E_vanderWaals)
      E_Trial = E_Trial + E_threebody 
   end subroutine
end module

