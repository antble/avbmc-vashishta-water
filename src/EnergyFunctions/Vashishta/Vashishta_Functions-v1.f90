!*********************************************************************************************************************
!     This file contains the energy functions that work for Vashishta style forcefields
!     these functions are enclosed inside of the module "InterMolecularEnergy" so that
!     the energy functions can be freely exchanged from the simulation.
!     The prefix naming scheme implies the following:
!           Detailed - Complete energy calculation intended for use at the beginning and end
!                      of the simulation.  This function is not inteded for use mid-simulation.
!             Shift  - Calculates the energy difference for any move that does not result
!                      in molecules being added or removed from the cluster. This function
!                      receives any number of Displacement vectors from the parent function as input.
!              Mol   - Calculates the energy for a molecule already present in the system. For
!                      use in moves such as particle deletion moves.
!             NewMol - Calculates the energy for a molecule that has been freshly inserted into the system.
!                      Intended for use in Rosenbluth Sampling, Swap In, etc.
!*********************************************************************************************************************
module InterEnergy_Vashishta
   use VarPrecision
contains
!======================================================================================
   subroutine two_body_interaction(H, eta, q, D, W, rc, r1s, r4s, r, steric, coulomb, chargeDipole, vanderWaals)
      
      use SimParameters, only : outputEConv 
      implicit none
      real(dp), intent(inout) :: H, eta, q, D, W, rc, r1s, r4s, r
      real(dp), intent(out) :: steric, coulomb, chargeDipole, vanderWaals
      real(dp) :: tmp, rcheta, rcqr1s, rcDr4s, rcW

      rcheta = H/(rc**eta)                         ! 2nd term is the truncation
      tmp = -eta*rcheta/rc                         ! 3rd term the necessary shift!  
      steric = H/(r**eta) - rcheta - (r-rc)*tmp    ! Full truncated (rc) and shifted (r<rc) potential
      steric = steric*outputEConv                  ! eV * (T/eV)

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
!======================================================================================
   subroutine Detailed_ECalc_Inter(E_T, PairList)
      use ParallelVar
      use ForceField
      use ForceFieldPara_Vashishta
      use Coords
      use SimParameters
      use EnergyTables
      implicit none
      
      real(dp), intent(inOut) :: E_T
      real(dp), intent(inOut) :: PairList(:, :)
      

      integer :: iType, jType, iMol, jMol, kMol, iAtom, jAtom, kAtom, iPair
      integer :: iIndx, jIndx, kIndx, globIndx1, globIndx2, globIndx3
      integer(kind=atomIntType) :: atmType1, atmType2
      real(dp) :: rx, ry, rz, r, rx1, ry1, rz1, rx2, ry2, rz2
      real(dp) :: rmin_ij
      real(dp) :: O_x1, O_y1, O_z1
      real(dp) :: H, q, eta, r4s, D, W    ! 2-body parameter
      real(dp) :: B, xi, r0, C, costheta  ! 3-body parameters
      real(dp) :: rcheta, rcqr1s, rcDr4s, rcW, tmp
      real(dp) :: r_ij, r_ik, costheta_ijk, test
      real(dp) :: E_steric, E_coulomb, E_chargeDipole, E_vanderWaals, E_twobody, E_threebody
      real(dp) :: steric, coulomb, chargeDipole, vanderWaals, twobody, threebody, E_Intra


      steric = 0E0_dp 
      coulomb = 0E0_dp 
      chargeDipole = 0E0_dp 
      vanderWaals = 0E0_dp 
      twobody = 0E0_dp
      threebody = 0E0_dp

      E_steric = 0E0 
      E_coulomb = 0E0 
      E_chargeDipole = 0E0 
      E_vanderWaals = 0E0 
      E_twobody = 0E0 
      E_threebody = 0E0 

      E_Inter_T = 0E0
      
      PairList = 0E0 
      ETable = 0E0_dp

      E_Intra = 0E0
      do iType = 1, nMolTypes
         do jType = iType, nMolTypes
            do iMol = 1, NPART(iType)
               do jMol = 1, NPART(jType)
                  iIndx = MolArray(iType)%mol(iMol)%indx
                  jIndx = MolArray(jType)%mol(jMol)%indx
                  do iAtom = 1, nAtoms(iType)
                     atmType1 = atomArray(1, iAtom)
                     do jAtom = 1, nAtoms(jType)
                        atmType2 = atomArray(1, jAtom)
                        if (iMol==jMol .and. (jAtom == iAtom .or. iAtom > jAtom)) then 
                           cycle 
                        end if

                        if (iMol < jMol .or. iMol == jMol) then
                           rmin_ij = r_min_tab(atmType1, atmType2)

                           q = q_tab(atmType1, atmType2)
                           H = H_tab(atmType1, atmType2)
                           eta = eta_tab(atmType1, atmType2)
                           r4s = r4s_tab(atmType1, atmType2)
                           D = D_tab(atmType1, atmType2)
                           W = W_tab(atmType1, atmType2)

                           rx = MolArray(iType)%mol(iMol)%x(iAtom) - MolArray(jType)%mol(jMol)%x(jAtom)
                           ry = MolArray(iType)%mol(iMol)%y(iAtom) - MolArray(jType)%mol(jMol)%y(jAtom)
                           rz = MolArray(iType)%mol(iMol)%z(iAtom) - MolArray(jType)%mol(jMol)%z(jAtom)
                           r = rx**2 + ry**2 + rz**2
                           if (distCriteria) then
                              if (iAtom .eq. 1) then
                                 if (jAtom .eq. 1) then
                                    PairList(iIndx, jIndx) = r
                                    PairList(jIndx, iIndx) = PairList(iIndx, jIndx)
                                 end if
                              end if
                           end if
                           if(r .lt. rmin_ij) then
                              print *, iAtom, iMol, jMol
                              stop "ERROR! Overlaping atoms found in the current configuration!"
                           endif 

                           r = sqrt(r)
                           if (r <= rc) then
                              call two_body_interaction(H, eta, q, D, W, rc, r1s, r4s, r, &
                                                      steric, coulomb, chargeDipole, vanderWaals)
                              E_steric = E_steric + steric
                              E_coulomb = E_coulomb  + coulomb
                              E_chargeDipole = E_chargeDipole + chargeDipole
                              E_vanderWaals = E_vanderWaals + vanderWaals

                              ETable(iIndx) = ETable(iIndx) + steric + coulomb + chargeDipole + vanderWaals
                              ETable(jIndx) = ETable(jIndx) + steric + coulomb + chargeDipole + vanderWaals
                           end if
                        end if
                     end do
                  end do
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
            do iMol = 1, NPART(iType)
               O_x1 = MolArray(iType)%mol(iMol)%x(1)
               O_y1 = MolArray(iType)%mol(iMol)%y(1)
               O_z1 = MolArray(iType)%mol(iMol)%z(1)
               iIndx = MolArray(iType)%mol(iMol)%indx
               do jMol = 1, NPART(jType)
                  jIndx = MolArray(jType)%mol(jMol)%indx
                  ! find an O-H that is within r0
                  ! start at 2 since Oxygen=1, note: dont mix up the vashishta forcefield O,H sequence in atom definition
                  do jAtom = 2, nAtoms(jType)  
                     rx1 = MolArray(iType)%mol(jMol)%x(jAtom) - O_x1
                     ry1 = MolArray(iType)%mol(jMol)%y(jAtom) - O_y1
                     rz1 = MolArray(iType)%mol(jMol)%z(jAtom) - O_z1 
                     r_ij = sqrt(rx1**2 + ry1**2 + rz1**2)
                     if (r_ij < r0) then
                        do kMol = 1, NPART(jType)
                           kIndx = MolArray(jType)%mol(kMol)%indx
                           ! find another O-H that is within r0
                           do kAtom = 2, nAtoms(jType)
                              if (((jMol == kMol) .and. (jAtom < kAtom)) .or. (jMol < kMol)) then
                                 rx2 = MolArray(iType)%mol(kMol)%x(kAtom) - O_x1
                                 ry2 = MolArray(iType)%mol(kMol)%y(kAtom) - O_y1
                                 rz2 = MolArray(iType)%mol(kMol)%z(kAtom) - O_z1 
                                 r_ik = sqrt(rx2**2 + ry2**2 + rz2**2)
                                 if (r_ik < r0) then
                                    call threebody_interaction(rx1, rx2, ry1, ry2, rz1, rz2, r_ij, r_ik, &
                                                               B, xi, r0, costheta, threebody)
                                    E_threebody = E_threebody + threebody
                                    ETable(iIndx) = ETable(iIndx) + threebody
                                    ETable(jIndx) = ETable(jIndx) + threebody
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

      !=============
      !use this to check the detailed energy values
      !=============
      ! write (nout, *) "Steric repulsion energy(eV):",E_steric/outputEConv
      ! write (nout, *) "Electrostatic Energy(eV):",E_coulomb/outputEConv
      ! write (nout, *) "Charge-dipole Energy(eV):",E_chargeDipole/outputEConv
      ! write (nout, *) "Vander Waals Energy(eV):",E_vanderWaals/outputEConv
      ! write (nout, *) "Threebody Energy(eV):",E_threebody/outputEConv
      ! write (nout, *) "Twobody Energy(eV):",(E_steric + E_coulomb + E_chargeDipole + E_vanderWaals)/outputEConv
      
      E_T = E_T + (E_steric + E_coulomb + E_chargeDipole + E_vanderWaals) 
      E_T = E_T + E_threebody
      E_Inter_T =  E_threebody + (E_steric + E_coulomb + E_chargeDipole + E_vanderWaals)
   end subroutine
!======================================================================================
   subroutine Shift_ECalc_Inter(E_Trial, disp, PairList, dETable, rejMove, update)
      use ForceField
      use ForceFieldPara_Vashishta
      use Coords
      use SimParameters
      implicit none

      type(Displacement), intent(in) :: disp(:)
      real(dp), intent(out) :: E_Trial
      real(dp), intent(inout) :: PairList(:), dETable(:)
      logical, intent(in) :: update
      logical, intent(out) :: rejMove

      integer :: iType, jType, iMol, iiMol, jMol, kMol, iAtom, jAtom, kAtom, iDisp, jDisp
      integer(kind=atomIntType) :: atmType1, atmType2, iIndx, jIndx, kIndx
      integer :: sizeDisp, dispMol, shifted_ctr
      real(dp) :: rx, ry, rz, r, rx1, ry1, rz1, rx2, ry2, rz2
      real(dp) :: r_ij, r_ik, r_new, r_old, rmin_ij, r_min1_sq
      real(dp) :: O_x1, O_y1, O_z1,O_x1_old, O_y1_old, O_z1_old, O_x1_new, O_y1_new, O_z1_new
      real(dp) :: H, eta, q, r4s, D, W                   ! 2-body parameter
      real(dp) :: B, xi, r0, C, costheta, costheta_ijk   ! 3-body parameters
      real(dp) :: E_steric, E_coulomb, E_chargeDipole, E_vanderWaals, E_twobody, E_threebody
      real(dp) :: E_threebody_new, E_threebody_old
      real(dp) :: steric, coulomb, chargeDipole, vanderWaals, twobody, threebody
      real(dp) :: rcheta, rcqr1s, rcDr4s, rcW, tmp, r1_dot_r2
      logical :: test 

      sizeDisp = size(disp)
      E_Trial = 0E0
      
      steric = 0E0 
      coulomb = 0E0 
      chargeDipole = 0E0 
      vanderWaals = 0E0 
      twobody = 0E0
      threebody = 0E0

      E_steric = 0E0 
      E_coulomb = 0E0 
      E_chargeDipole = 0E0 
      E_vanderWaals = 0E0 
      E_twobody = 0E0
      E_threebody = 0E0
      E_threebody_old = 0E0
      E_threebody_new  = 0E0

      dETable = 0E0
      E_Trial = 0E0 
      PairList = 0E0
      iType = disp(1)%molType
      iMol = disp(1)%molIndx
      iIndx = MolArray(iType)%mol(iMol)%Indx

      test = .false.
!  This section calculates the "Intramolecular" interaction between the atoms that
!  have been modified in this trial move with the atoms that have remained stationary
      do iDisp = 1, sizeDisp-1
         iAtom = disp(iDisp)%atmIndx
         atmType1 = atomArray(iType, iAtom)

         do jType = 1, nMolTypes
            do jDisp = iDisp+1, sizeDisp
               jAtom = disp(jDisp)%atmIndx
               atmType2 = atomArray(jType, jAtom)
               rmin_ij = r_min_tab(atmType2, atmType1)
               q = q_tab(atmType2, atmType1)
               H = H_tab(atmType2, atmType1)
               eta = eta_tab(atmType2, atmType1)
               r4s = r4s_tab(atmType2, atmType1)
               D = D_tab(atmType2, atmType1)
               W = W_tab(atmType2, atmType1)

               ! Distance for the New position
               rx = disp(iDisp)%x_new - disp(jDisp)%x_new
               ry = disp(iDisp)%y_new - disp(jDisp)%y_new
               rz = disp(iDisp)%z_new - disp(jDisp)%z_new
               r_new = rx*rx + ry*ry + rz*rz

               ! If r_new is less than r_min reject the move.
               if (r_new .lt. rmin_ij) then
                  rejMove = .true.
                  return
               end if
               if(r_new .lt. rmin_ij) then
                  print *, iAtom, iMol, jMol
                  stop "SHIFT ERROR! Overlaping atoms found in the current configuration!"
               endif 
               
               jIndx = MolArray(jType)%mol(iMol)%indx
               if (distCriteria) then
                  if (iAtom .eq. 1) then
                     if (jAtom .eq. 1) then
                        PairList(jIndx) = r_new
                     end if
                  end if
               end if

               ! Distance for the Old position
               rx = disp(iDisp)%x_old - disp(jDisp)%x_old
               ry = disp(iDisp)%y_old - disp(jDisp)%y_old
               rz = disp(iDisp)%z_old - disp(jDisp)%z_old
               r_old = rx*rx + ry*ry + rz*rz

               ! new position Energy calculations 
               r_new = sqrt(r_new)
               if (r_new <= rc) then
                  call two_body_interaction(H, eta, q, D, W, rc, r1s, r4s, r_new, &
                                          steric, coulomb, chargeDipole, vanderWaals)
                  E_steric = E_steric + steric
                  E_coulomb = E_coulomb  + coulomb
                  E_chargeDipole = E_chargeDipole + chargeDipole
                  E_vanderWaals = E_vanderWaals + vanderWaals 
                  
                  dETable(iIndx) = dETable(iIndx) + steric + coulomb  + chargeDipole + vanderWaals
                  dETable(jIndx) = dETable(jIndx) + (steric + coulomb  + chargeDipole + vanderWaals)
                  
               end if
               
               ! old position Energy calculations 
               r_old = sqrt(r_old)
               if (r_old <= rc) then
                  call two_body_interaction(H, eta, q, D, W, rc, r1s, r4s, r_old, &
                                          steric, coulomb, chargeDipole, vanderWaals)
                  E_steric = E_steric - steric
                  E_coulomb = E_coulomb  - coulomb
                  E_chargeDipole = E_chargeDipole - chargeDipole
                  E_vanderWaals = E_vanderWaals - vanderWaals 
                  
                  dETable(iIndx) = dETable(iIndx) - (steric + coulomb  + chargeDipole + vanderWaals)
                  dETable(jIndx) = dETable(jIndx) - (steric + coulomb  + chargeDipole + vanderWaals)
                  
               end if
            end do 
         end do 
      end do 
!  This section calculates the Intermolecular interaction between the atoms that
!  have been modified in this trial move with the atoms that have remained stationary
      do iDisp = 1, sizeDisp
         iAtom = disp(iDisp)%atmIndx
         
         atmType1 = atomArray(iType, iAtom)
         do jType = 1, nMolTypes
            do jAtom = 1, nAtoms(jType)
               atmType2 = atomArray(jType, jAtom)
               rmin_ij = r_min_tab(atmType2, atmType1)
               q = q_tab(atmType2, atmType1)
               H = H_tab(atmType2, atmType1)
               eta = eta_tab(atmType2, atmType1)
               r4s = r4s_tab(atmType2, atmType1)
               D = D_tab(atmType2, atmType1)
               W = W_tab(atmType2, atmType1)
               do jMol = 1, NPART(jType)
                  if (iType .eq. jType) then
                     if (iMol .eq. jMol) then
                        cycle
                     end if
                  end if
                  ! Distance for the New position
                  rx = disp(iDisp)%x_new - MolArray(jType)%mol(jMol)%x(jAtom)
                  ry = disp(iDisp)%y_new - MolArray(jType)%mol(jMol)%y(jAtom)
                  rz = disp(iDisp)%z_new - MolArray(jType)%mol(jMol)%z(jAtom)
                  r_new = rx*rx + ry*ry + rz*rz

                  ! If r_new is less than r_min reject the move.
                  if (r_new .lt. rmin_ij) then
                     rejMove = .true.
                     return
                  end if
                  jIndx = MolArray(jType)%mol(jMol)%indx
                  if (distCriteria) then
                     if (iAtom .eq. 1) then
                        if (jAtom .eq. 1) then
                           PairList(jIndx) = r_new
                        end if
                     end if
                  end if

                  ! Distance for the Old position
                  rx = disp(iDisp)%x_old - MolArray(jType)%mol(jMol)%x(jAtom)
                  ry = disp(iDisp)%y_old - MolArray(jType)%mol(jMol)%y(jAtom)
                  rz = disp(iDisp)%z_old - MolArray(jType)%mol(jMol)%z(jAtom)
                  r_old = rx*rx + ry*ry + rz*rz
                  
                  ! new position Energy calculations 
                  r_new = sqrt(r_new)
                  if (r_new <= rc) then
                     call two_body_interaction(H, eta, q, D, W, rc, r1s, r4s, r_new, &
                                          steric, coulomb, chargeDipole, vanderWaals)
                     E_steric = E_steric + steric
                     E_coulomb = E_coulomb  + coulomb
                     E_chargeDipole = E_chargeDipole + chargeDipole
                     E_vanderWaals = E_vanderWaals + vanderWaals 
                     
                     dETable(iIndx) = dETable(iIndx) + steric + coulomb  + chargeDipole + vanderWaals
                     dETable(jIndx) = dETable(jIndx) + steric + coulomb  + chargeDipole + vanderWaals
                  end if
                  
                  ! old position Energy calculations 
                  r_old = sqrt(r_old)
                  if (r_old <= rc) then
                     call two_body_interaction(H, eta, q, D, W, rc, r1s, r4s, r_old, &
                                          steric, coulomb, chargeDipole, vanderWaals)
                     E_steric = E_steric - steric
                     E_coulomb = E_coulomb  - coulomb
                     E_chargeDipole = E_chargeDipole - chargeDipole
                     E_vanderWaals = E_vanderWaals - vanderWaals
                     
                     dETable(iIndx) = dETable(iIndx) - (steric + coulomb  + chargeDipole + vanderWaals)
                     dETable(jIndx) = dETable(jIndx) - (steric + coulomb  + chargeDipole + vanderWaals)
                     
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

! Determines the change of the 3-body energy if you shift a molecule
      dispMol = disp(1)%molIndx  ! index of the displaced molecule

      ! NEW CONFIGURATION - 3BODY INTRXN ENERGY
      do iType = 1, nMolTypes
         do jType = iType, nMolTypes
            do iiMol = 1, NPART(iType)
               iIndx = MolArray(iType)%mol(iiMol)%indx
               if (iiMol == dispMol) then
                  O_x1 = disp(1)%x_new
                  O_y1 = disp(1)%y_new
                  O_z1 = disp(1)%z_new
               else
                  O_x1 = MolArray(iType)%mol(iiMol)%x(1)
                  O_y1 = MolArray(iType)%mol(iiMol)%y(1)
                  O_z1 = MolArray(iType)%mol(iiMol)%z(1)
               end if 
               
               
               do jMol = 1, NPART(jType)
                  jIndx = MolArray(jType)%mol(jMol)%indx !?? 
                  ! find an O-H that is within r0
                  ! start at 2 since Oxygen=1, note: dont mix up the vashishta forcefield O,H sequence in atom definition
                  do jAtom = 2, nAtoms(jType)  
                     if (jMol == dispMol) then 
                        rx1 = O_x1 - disp(jAtom)%x_new 
                        ry1 = O_y1 - disp(jAtom)%y_new
                        rz1 = O_z1 - disp(jAtom)%z_new 
                     else
                        rx1 = O_x1 - MolArray(jType)%mol(jMol)%x(jAtom) 
                        ry1 = O_y1 - MolArray(jType)%mol(jMol)%y(jAtom) 
                        rz1 = O_z1 - MolArray(jType)%mol(jMol)%z(jAtom) 
                     end if 
                     r_ij = sqrt(rx1**2 + ry1**2 + rz1**2)

                     if (r_ij < r0) then
                        do kMol = 1, NPART(jType)
                           kIndx = MolArray(jType)%mol(kMol)%indx !?? 
                           ! find another O-H that is within r0
                           do kAtom = 2, nAtoms(jType)
                              if (((jMol == kMol) .and. (jAtom < kAtom)) .or. (jMol < kMol)) then
                                 if (kMol == dispMol) then 
                                    rx2 = O_x1 - disp(kAtom)%x_new
                                    ry2 = O_y1 - disp(kAtom)%y_new
                                    rz2 = O_z1 - disp(kAtom)%z_new
                                 else 
                                    rx2 = O_x1 - MolArray(jType)%mol(kMol)%x(kAtom)
                                    ry2 = O_y1 - MolArray(jType)%mol(kMol)%y(kAtom)
                                    rz2 = O_z1 - MolArray(jType)%mol(kMol)%z(kAtom)
                                 end if 
                                 r_ik = sqrt(rx2**2 + ry2**2 + rz2**2)
                                 if (r_ik < r0) then
                                    call threebody_interaction(rx1, rx2, ry1, ry2, rz1, rz2, r_ij, r_ik, &
                                                               B, xi, r0, costheta, threebody)
                                    E_threebody_new = E_threebody_new + threebody 
                                    dETable(iIndx) = dETable(iIndx) + threebody
                                    dETable(jIndx) = dETable(jIndx) + threebody
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
      ! OLD CONFIGURATION - 3BODY INTRXN ENERGY
      do iType = 1, nMolTypes
         do jType = iType, nMolTypes
            do iMol = 1, NPART(iType)
               iIndx = MolArray(jType)%mol(iMol)%indx
               O_x1 = MolArray(iType)%mol(iMol)%x(1)
               O_y1 = MolArray(iType)%mol(iMol)%y(1)
               O_z1 = MolArray(iType)%mol(iMol)%z(1)

               do jMol = 1, NPART(jType)
                  jIndx = MolArray(jType)%mol(jMol)%indx !?? 
                  ! find an O-H that is within r0
                  ! start at 2 since Oxygen=1, note: dont mix up the vashishta forcefield O,H sequence in atom definition
                  do jAtom = 2, nAtoms(jType)  
                     rx1 = O_x1 - MolArray(jType)%mol(jMol)%x(jAtom)
                     ry1 = O_y1 - MolArray(jType)%mol(jMol)%y(jAtom)
                     rz1 = O_z1 - MolArray(jType)%mol(jMol)%z(jAtom)
   
                     r_ij = sqrt(rx1**2 + ry1**2 + rz1**2)
                     if (r_ij < r0) then
                        do kMol = 1, NPART(jType)
                           kIndx = MolArray(jType)%mol(kMol)%indx !??
                           ! find another O-H that is within r0
                           do kAtom = 2, nAtoms(jType)
                              if (((jMol == kMol) .and. (jAtom < kAtom)) .or. (jMol < kMol)) then
                                 rx2 = O_x1 - MolArray(jType)%mol(kMol)%x(kAtom)
                                 ry2 = O_y1 - MolArray(jType)%mol(kMol)%y(kAtom)
                                 rz2 = O_z1 - MolArray(jType)%mol(kMol)%z(kAtom) 
                                 r_ik = sqrt(rx2**2 + ry2**2 + rz2**2)

                                 if (r_ik < r0) then
                                    call threebody_interaction(rx1, rx2, ry1, ry2, rz1, rz2, r_ij, r_ik, &
                                                               B, xi, r0, costheta, threebody)
                                    E_threebody_old = E_threebody_old  + threebody
                                    
                                    dETable(iIndx) = dETable(iIndx) - threebody
                                    dETable(jIndx) = dETable(jIndx) - threebody
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
      E_Trial = (E_threebody_new - E_threebody_old) + E_Trial 
      
   end subroutine
!======================================================================================
   subroutine Mol_ECalc_Inter(iType, iMol, dETable, E_Trial, E_intra)
      use ForceField
      use ForceFieldPara_Vashishta
      use Coords
      use SimParameters
      implicit none
      integer, intent(in) :: iType, iMol
      real(dp), intent(out) :: E_Trial
      real(dp), intent(inout) :: dETable(:), E_intra

      integer :: iAtom, iIndx, jType, jIndx, iiMol, jMol, kMol, jAtom, kAtom, kIndx
      integer :: jAtom11, iAtom22
      integer(kind=atomIntType)  :: atmType1, atmType2
      real(dp) :: rx, ry, rz, r
      real(dp) :: E_steric, E_coulomb, E_chargeDipole, E_vanderWaals, E_twobody, E_threebody
      real(dp) :: steric, coulomb, chargeDipole, vanderWaals, twobody, threebody
      real(dp) :: H, eta, r4s, D, W ! 2-body parameter
      real(dp) :: B, xi, r0, q, C, costheta, costheta_ijk  ! 3-body parameters
      real(dp) :: rcheta, rcqr1s, rcDr4s, rcW, tmp
      real(dp) :: rx1, ry1, rz1, rx2, ry2, rz2, r_ij, r_ik
      real(dp) :: O_x1, O_y1, O_z1

      steric = 0E0 
      coulomb = 0E0 
      chargeDipole = 0E0 
      vanderWaals = 0E0 
      twobody = 0E0
      threebody = 0E0

      E_steric = 0E0 
      E_coulomb = 0E0 
      E_chargeDipole = 0E0 
      E_vanderWaals = 0E0 
      E_twobody = 0E0
      E_threebody = 0E0 

      
      E_Trial = 0E0
      dETable = 0E0
      
      iIndx = MolArray(iType)%mol(iMol)%indx
      ! 2-body intermolecular interaction 
      do iAtom = 1, nAtoms(iType)
         atmType1 = atomArray(iType, iAtom)
         do jType = 1, nMolTypes
            do jAtom = 1, nAtoms(jType)
               atmType2 = atomArray(jType, jAtom)

               q = q_tab(atmType2, atmType1)
               H = H_tab(atmType2, atmType1)
               eta = eta_tab(atmType2, atmType1)
               r4s = r4s_tab(atmType2, atmType1)
               D = D_tab(atmType2, atmType1)
               W = W_tab(atmType2, atmType1)

               do jMol = 1, NPART(jType)
                  
!             New Energy Calculation
                  if (iMol == jMol .and. iAtom == jAtom) then ! .or. iAtom > jAtom) then
                     cycle
                  end if
                  jIndx = MolArray(jType)%mol(jMol)%indx
                  if ((iMol == jMol .and. iAtom  < jAtom) .or. iMol .ne. jMol) then
                     rx = MolArray(iType)%mol(iMol)%x(iAtom) - MolArray(jType)%mol(jMol)%x(jAtom)
                     ry = MolArray(iType)%mol(iMol)%y(iAtom) - MolArray(jType)%mol(jMol)%y(jAtom)
                     rz = MolArray(iType)%mol(iMol)%z(iAtom) - MolArray(jType)%mol(jMol)%z(jAtom)
                     r = rx*rx + ry*ry + rz*rz
                     
                     r = sqrt(r)
                     if (r <= rc) then
                        call two_body_interaction(H, eta, q, D, W, rc, r1s, r4s, r, &
                                          steric, coulomb, chargeDipole, vanderWaals)
                        E_steric = E_steric + steric
                        E_coulomb = E_coulomb  + coulomb
                        E_chargeDipole = E_chargeDipole + chargeDipole
                        E_vanderWaals = E_vanderWaals + vanderWaals
                        ! record the intramolecular energy of the molecule to be removed
                        if (iMol == jMol) then 
                           E_intra = E_intra + steric + coulomb  + chargeDipole + vanderWaals
                        end if 
                        dETable(iIndx) = dETable(iIndx) + steric + coulomb  + chargeDipole + vanderWaals
                        dETable(jIndx) = dETable(jIndx) + steric + coulomb  + chargeDipole + vanderWaals 
                     end if
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

      do jType = iType, nMolTypes
         do iiMol = 1, NPART(iType)
            O_x1 = MolArray(iType)%mol(iiMol)%x(1)
            O_y1 = MolArray(iType)%mol(iiMol)%y(1)
            O_z1 = MolArray(iType)%mol(iiMol)%z(1)
            iIndx = MolArray(jType)%mol(iiMol)%indx
            do jMol = 1, NPART(jType)
               jIndx = MolArray(jType)%mol(jMol)%indx
               ! find an O-H that is within r0
               ! start at 2 since Oxygen=1, note: dont mix up the vashishta forcefield O,H sequence in atom definition
               do jAtom = 2, nAtoms(jType)  
                  rx1 = MolArray(iType)%mol(jMol)%x(jAtom) - O_x1
                  ry1 = MolArray(iType)%mol(jMol)%y(jAtom) - O_y1
                  rz1 = MolArray(iType)%mol(jMol)%z(jAtom) - O_z1 
                  r_ij = sqrt(rx1**2 + ry1**2 + rz1**2)
                  if (r_ij < r0) then
                     do kMol = 1, NPART(jType)
                        kIndx = MolArray(jType)%mol(kMol)%indx
                        if (iType .eq. jType) then
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
                              rx2 = MolArray(iType)%mol(kMol)%x(kAtom) - O_x1
                              ry2 = MolArray(iType)%mol(kMol)%y(kAtom) - O_y1
                              rz2 = MolArray(iType)%mol(kMol)%z(kAtom) - O_z1 
                              r_ik = sqrt(rx2**2 + ry2**2 + rz2**2)
                              
                              if (r_ik < r0) then
                                 call threebody_interaction(rx1, rx2, ry1, ry2, rz1, rz2, r_ij, r_ik, &
                                                            B, xi, r0, costheta, threebody)
                                 E_threebody = E_threebody + threebody
                                 if (iMol == jMol .and. jMol == kMol) then 
                                    E_intra = E_intra + E_threebody
                                 end if 
                                 dETable(iIndx) = dETable(iIndx) + threebody
                                 dETable(jIndx) = dETable(jIndx) + threebody
                              end if 
                           end if 
                        end do 
                     end do 
                  end if 
               end do 
            end do
         end do 
      end do 

      E_Trial = (E_steric + E_coulomb + E_chargeDipole + E_vanderWaals) + E_threebody
   end subroutine
!======================================================================================
   subroutine NewMol_ECalc_Inter(E_Trial, PairList, dETable, rejMove, E_intra)
      use ForceField
      use ForceFieldPara_Vashishta
      use Coords
      use SimParameters
      implicit none
      logical, intent(out) :: rejMove
      real(dp), intent(out) :: E_Trial
      real(dp), intent(inout) :: PairList(:), dETable(:), E_intra

      integer :: iAtom, iIndx, iType, jType, jIndx, kIndx, iMol, jMol, kMol, jAtom, kAtom, ctr
      integer(kind=atomIntType) :: atmType1, atmType2
      real(dp) :: rx, ry, rz, r, rx1, ry1, rz1, rx2, ry2, rz2, r_ij, r_ik
      real(dp) :: rmin_ij
      real(dp) :: E_steric, E_coulomb, E_chargeDipole, E_vanderWaals, E_twobody, E_threebody
      real(dp) :: steric, coulomb, chargeDipole, vanderWaals, twobody, threebody
      real(dp) :: O_x1, O_y1, O_z1
      real(dp) :: H, eta, r4s, D, W                         ! 2-body parameters
      real(dp) :: B, xi, r0, q, C, costheta, costheta_ijk   ! 3-body parameters
      real(dp) :: rcheta, rcqr1s, rcDr4s, rcW, tmp

      steric = 0E0 
      coulomb = 0E0
      chargeDipole = 0E0 
      vanderWaals = 0E0 
      twobody = 0E0
      threebody = 0E0 

      E_steric = 0E0 
      E_coulomb = 0E0 
      E_chargeDipole = 0E0 
      E_vanderWaals = 0E0 
      E_twobody = 0E0
      E_threebody = 0E0

      dETable = 0E0
      PairList = 0E0
      E_Trial = 0E0
      rejMove = .false.

      iIndx = molArray(newMol%molType)%mol(NPART(newMol%molType) + 1)%indx
      jIndx = molArray(newMol%molType)%mol(NPART(newMol%molType) + 1)%indx
      ! "intramolecular" contribution of the added atom
      do iAtom = 1, nAtoms(newMol%molType)-1
         atmType1 = atomArray(newMol%molType, iAtom)
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

               rx = newMol%x(iAtom) - newMol%x(jAtom)
               ry = newMol%y(iAtom) - newMol%y(jAtom)
               rz = newMol%z(iAtom) - newMol%z(jAtom)
               r = rx*rx + ry*ry + rz*rz

               r = sqrt(r)
               if (r <= rc) then
                  call two_body_interaction(H, eta, q, D, W, rc, r1s, r4s, r, &
                                          steric, coulomb, chargeDipole, vanderWaals)
                  E_steric = E_steric + steric
                  E_coulomb = E_coulomb  + coulomb
                  E_chargeDipole = E_chargeDipole + chargeDipole
                  E_vanderWaals = E_vanderWaals + vanderWaals
                  dETable(iIndx) = dETable(iIndx) + steric + coulomb  + chargeDipole + vanderWaals 
                  dETable(jIndx) = dETable(jIndx) + steric + coulomb  + chargeDipole + vanderWaals
               end if 
            end do
         end do 
      end do 

      ! record two-body intramolecular energy
      E_intra =  E_steric + E_coulomb + E_chargeDipole + E_vanderWaals

      ! intermolecular energy contribution 
      do iAtom = 1, nAtoms(newMol%molType)
         atmType1 = atomArray(newMol%molType, iAtom)
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

                  rx = newMol%x(iAtom) - MolArray(jType)%mol(jMol)%x(jAtom)
                  ry = newMol%y(iAtom) - MolArray(jType)%mol(jMol)%y(jAtom)
                  rz = newMol%z(iAtom) - MolArray(jType)%mol(jMol)%z(jAtom)
                  r = rx*rx + ry*ry + rz*rz

                  if (r .lt. rmin_ij) then
                     rejMove = .true.
                     return
                  end if
                  if (distCriteria) then
                     if (iAtom .eq. 1) then
                        if (jAtom .eq. 1) then
                           PairList(jIndx) = r
                        end if
                     end if
                  end if

                  r = sqrt(r)
                  if (r <= rc) then
                     call two_body_interaction(H, eta, q, D, W, rc, r1s, r4s, r, &
                                          steric, coulomb, chargeDipole, vanderWaals)
                     E_steric = E_steric + steric
                     E_coulomb = E_coulomb  + coulomb
                     E_chargeDipole = E_chargeDipole + chargeDipole
                     E_vanderWaals = E_vanderWaals + vanderWaals 
                     dETable(iIndx) = dETable(iIndx) + steric + coulomb  + chargeDipole + vanderWaals
                     dETable(jIndx) = dETable(jIndx) + steric + coulomb  + chargeDipole + vanderWaals
                     
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
                  iIndx = molArray(jType)%mol(iMol)%indx
                  ! Set Oxygen atom as the anchor
                  if (iMol == NPART(1) + 1) then 
                     O_x1 = newMol%x(1)
                     O_y1 = newMol%y(1)
                     O_z1 = newMol%z(1)
                  else 
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
                           rx1 = newMol%x(jAtom) - O_x1
                           ry1 = newMol%y(jAtom) - O_y1
                           rz1 = newMol%z(jAtom) - O_z1 
                           r_ij = sqrt(rx1**2 + ry1**2 + rz1**2 + 0.0_dp)
                        else 
                           rx1 = MolArray(iType)%mol(jMol)%x(jAtom) - O_x1
                           ry1 = MolArray(iType)%mol(jMol)%y(jAtom) - O_y1
                           rz1 = MolArray(iType)%mol(jMol)%z(jAtom) - O_z1 
                           r_ij = sqrt(rx1**2 + ry1**2 + rz1**2 + 0.0_dp)
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
                                       rx2 = newMol%x(kAtom) - O_x1
                                       ry2 = newMol%y(kAtom) - O_y1
                                       rz2 = newMol%z(kAtom) - O_z1 
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
                                       ! record the 3body intramolecular energy contribution
                                       if (iMol == jMol .and. jMol == kMol) then 
                                          E_intra = E_intra + E_threebody
                                       end if 

                                       dETable(iIndx) = dETable(iIndx) + threebody
                                       dETable(jIndx) = dETable(jIndx) + threebody
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
      E_Trial = (E_steric + E_coulomb + E_chargeDipole + E_vanderWaals) + E_threebody
   end subroutine
!======================================================================================
end module

