!*********************************************************************************************************************
!     This file contains the pressure functions that work for Lennard-Jones w/ Columbic style forcefields
!     these functions are enclosed inside of the module "InterMolecularEnergy" so that
!     the energy functions can be freely exchanged from the simulation.
!     The prefix naming scheme implies the following:
!           Detailed - Complete energy calculation inteded for use at the beginning and end
!                      of the simulation.  This function is not inteded for use mid-simulation.
!             Shift  - Calculates the energy difference for any move that does not result
!                      in molecules being added or removed from the cluster. This function
!                      receives any number of Displacement vectors from the parent function as input.
!              Mol   - Calculates the energy for a molecule already present in the system. For
!                      use in moves such as particle deletion moves.
!             NewMol - Calculates the energy for a molecule that has been freshly inserted into the system.
!                      Intended for use in Rosenbluth Sampling, Swap In, etc.
!*********************************************************************************************************************
module Pressure_LJ_Electro
   use VarPrecision
   use InterEnergy_LJ_Electro, only: lj_cut, lj_Cut_sq, q_cut, q_cut_sq

contains
!======================================================================================
   pure function LJ_Func(r_sq, ep, sig) result(LJ)
      implicit none
      real(dp), intent(in) :: r_sq, ep, sig
      real(dp) :: LJ

      LJ = (sig/r_sq)
      LJ = LJ*LJ*LJ
      LJ = 12d0*ep*LJ*(LJ - 0.5E0_dp)

   end function
!======================================================================================
   pure function Ele_Func(r, q) result(Ele)
      implicit none
      real(dp), intent(in) :: r, q
      real(dp) :: Ele

!      r = sqrt(r_sq)
      Ele = q/r

   end function
!======================================================================================
   subroutine Detailed_PressCalc_Inter(P_T)
!      use ParallelVar
      use ForceField, only: nAtoms, atomArray
      use ForceFieldPara_LJ_Q, only: ep_tab, q_tab, sig_tab
      use Coords, only: MolArray
      use SimParameters, only: nMolTypes, NPART, distCriteria
      use PairStorage, only: rPair, distStorage
      use ParallelVar, only: nout
      implicit none
      real(dp), intent(out) :: P_T

      integer :: iType, jType, iMol, jMol, iAtom, jAtom
      integer(kind=atomIntType) :: atmType1, atmType2
      integer :: iIndx, jIndx, globIndx1, globIndx2, jMolMin
      real(dp) :: r_sq, r
      real(dp) :: ep, sig_sq, q
      real(dp) :: LJ, Ele
      real(dp) :: P_Ele, P_LJ

      P_LJ = 0E0_dp
      P_Ele = 0E0_dp
      P_T = 0E0_dp
      do iType = 1, nMolTypes
         do jType = iType, nMolTypes
            do iMol = 1, NPART(iType)
               if (iType .eq. jType) then
                  jMolMin = iMol + 1
               else
                  jMolMin = 1
               end if
               do jMol = jMolMin, NPART(jType)
                  iIndx = MolArray(iType)%mol(iMol)%indx
                  jIndx = MolArray(jType)%mol(jMol)%indx
                  do iAtom = 1, nAtoms(iType)
                     atmType1 = atomArray(iType, iAtom)
                     globIndx1 = MolArray(iType)%mol(iMol)%globalIndx(iAtom)
                     do jAtom = 1, nAtoms(jType)
                        atmType2 = atomArray(jType, jAtom)
                        ep = ep_tab(atmType1, atmType2)
                        q = q_tab(atmType1, atmType2)
                        sig_sq = sig_tab(atmType1, atmType2)
                        globIndx2 = MolArray(jType)%mol(jMol)%globalIndx(jAtom)
                        r_sq = rPair(globIndx1, globIndx2)%p%r_sq
!                 if(r_sq .gt. lj_Cut_sq) then
!                   cycle
!                 endif
                        if (r_sq .lt. lj_cut_sq) then
                           LJ = LJ_Func(r_sq, ep, sig_sq)
                           P_LJ = P_LJ + LJ
                        end if
                        r = sqrt(r_sq)
                        Ele = q/r
                        P_Ele = P_Ele + Ele
                     end do
                  end do
               end do
            end do
         end do
      end do

      write (nout, *) "Lennard-Jones Pressure:", P_LJ
      write (nout, *) "Eletrostatic Pressure:", P_Ele
      P_T = P_Ele + P_LJ
      write (nout, *) "Total Pressure:", P_T
   end subroutine
!======================================================================================
   subroutine Shift_PressCalc_Inter(P_Trial, disp)
      use Coords, only: Displacement, atomIndicies, molArray
      use ForceField, only: nAtoms, atomArray
      use ForceFieldPara_LJ_Q, only: ep_tab, q_tab, sig_tab
      use SimParameters, only: distCriteria
      use PairStorage, only: distStorage, rPair, rPairNew, newDist, DistArrayNew, nNewDist, oldIndxArray
      implicit none

      type(Displacement), intent(in) :: disp(:)
      real(dp), intent(out) :: P_Trial

      integer :: iType, jType, iMol, jMol, iAtom, jAtom, iPair
      integer(kind=atomIntType) :: atmType1, atmType2, iIndx, jIndx
      integer :: sizeDisp
!      integer, pointer :: oldIndx
      integer :: globIndx1, globIndx2
      real(dp) :: r_new, r_new_sq, r, r_sq
      real(dp) :: ep, sig_sq, q
      real(dp) :: LJ, Ele, P_New, P_PairOld, P_Old
      real(dp) :: P_Ele, P_LJ

      P_LJ = 0E0_dp
      P_Ele = 0E0_dp
      P_Trial = 0E0_dp
      P_Old = 0E0_dp
      iType = disp(1)%molType
      iMol = disp(1)%molIndx
      iIndx = MolArray(iType)%mol(iMol)%indx

!      !This section calculates the Intermolecular interaction between the atoms that
!      !have been modified in this trial move with the atoms that have remained stationary

      do iPair = 1, nNewDist

         globIndx1 = newDist(iPair)%indx1
         globIndx2 = newDist(iPair)%indx2
         if (.not. rPair(globIndx1, globIndx2)%p%usePair) then
            cycle
         end if
         jType = atomIndicies(globIndx2)%nType
         jMol = atomIndicies(globIndx2)%nMol
         jIndx = MolArray(jType)%mol(jMol)%indx
         jAtom = atomIndicies(globIndx2)%nAtom
         iAtom = atomIndicies(globIndx1)%nAtom

         atmType1 = atomArray(iType, iAtom)
         atmType2 = atomArray(jType, jAtom)

         ep = ep_tab(atmType2, atmType1)
         q = q_tab(atmType2, atmType1)
         LJ = 0E0_dp
         if (ep .ne. 0E0_dp) then
            sig_sq = sig_tab(atmType2, atmType1)
            r_new_sq = rPairNew(globIndx1, globIndx2)%p%r_sq
            if (r_new_sq .lt. lj_cut_sq) then
               LJ = LJ + LJ_Func(r_new_sq, ep, sig_sq)
            end if
            r_sq = rPair(globIndx1, globIndx2)%p%r_sq
            if (r_sq .lt. lj_cut_sq) then
               LJ = LJ - LJ_Func(r_sq, ep, sig_sq)
            end if
         end if
         Ele = 0E0_dp
         if (q .ne. 0E0_dp) then
            r_new = rPairNew(globIndx1, globIndx2)%p%r
            Ele = q/r_new
            r = rPair(globIndx1, globIndx2)%p%r
            Ele = Ele - q/r
         end if
         P_Trial = P_Trial + Ele + LJ
         if (abs(P_Trial) > 1E6) then
            write (*, *) r_new, r, r_sq, LJ, Ele, q, ep, sig_sq, rPair(globIndx1, globIndx2)%p%storeRValue
         end if
      end do

   end subroutine
!======================================================================================
   pure subroutine Mol_PressCalc_Inter(iType, iMol, P_Trial)
      use ForceField, only: nAtoms, atomArray
      use ForceFieldPara_LJ_Q, only: ep_tab, q_tab, sig_tab
      use Coords, only: Displacement, atomIndicies, molArray
      use SimParameters, only: distCriteria, nMolTypes, NPART
      use PairStorage, only: distStorage, rPair
      implicit none
      integer, intent(in) :: iType, iMol
      real(dp), intent(out) :: P_Trial

      integer :: iAtom, iIndx, jType, jIndx, jMol, jAtom
      integer  :: globIndx1, globIndx2
      integer  :: atmType1, atmType2
      real(dp) :: ep, sig_sq, q
      real(dp) :: Ele, LJ, P
      real(dp) :: r_sq, r

      P_Trial = 0E0_dp
      iIndx = MolArray(iType)%mol(iMol)%indx

      do iAtom = 1, nAtoms(iType)
         atmType1 = atomArray(iType, iAtom)
         globIndx1 = MolArray(iType)%mol(iMol)%globalIndx(iAtom)
         do jType = 1, nMolTypes
            do jMol = 1, NPART(jType)
               jIndx = MolArray(jType)%mol(jMol)%indx
               if (iIndx .eq. jIndx) then
                  cycle
               end if
               do jAtom = 1, nAtoms(jType)
                  globIndx2 = MolArray(jType)%mol(jMol)%globalIndx(jAtom)
                  if (.not. rPair(globIndx1, globIndx2)%p%usePair) then
                     cycle
                  end if
                  atmType2 = atomArray(jType, jAtom)
                  ep = ep_tab(atmType2, atmType1)
                  q = q_tab(atmType2, atmType1)
                  LJ = 0E0_dp
!              r_sq = rPair(globIndx1, globIndx2)%p%r_sq
                  if (ep .ne. 0E0_dp) then
                     r_sq = rPair(globIndx1, globIndx2)%p%r_sq
                     if (r_sq .lt. lj_Cut_sq) then
                        sig_sq = sig_tab(atmType2, atmType1)
                        LJ = LJ_Func(r_sq, ep, sig_sq)
                     end if
                  end if

                  Ele = 0E0_dp
                  if (q .ne. 0E0_dp) then
                     r = rPair(globIndx1, globIndx2)%p%r
                     Ele = q/r
                  end if
                  P_Trial = P_Trial + Ele + LJ
               end do
            end do
         end do
      end do

   end subroutine
!======================================================================================
   subroutine NewMol_PressCalc_Inter(P_Trial)
      use ForceField, only: nAtoms, atomArray
      use ForceFieldPara_LJ_Q, only: ep_tab, q_tab, sig_tab
      use Coords, only: Displacement, atomIndicies, molArray, newmol
      use SimParameters, only: distCriteria, nMolTypes, NPART
      use PairStorage, only: distStorage, rPair, rPairNew, newDist, nNewDist
      implicit none
      real(dp), intent(out) :: P_Trial

      integer :: iPair
      integer :: iType, iMol, iAtom, iIndx, jType, jIndx, jMol, jAtom
      integer(kind=atomIntType) :: atmType1, atmType2
      integer :: globIndx1, globIndx2
      real(dp) :: r_sq, r
      real(dp) :: ep, sig_sq, q
      real(dp) :: LJ, Ele

      real(dp) :: P_Ele, P_LJ

      P_LJ = 0E0_dp
      P_Ele = 0E0_dp
      P_Trial = 0E0_dp

      iType = newMol%molType
      iMol = NPART(iType) + 1
      iIndx = molArray(iType)%mol(iMol)%indx
      do iPair = 1, nNewDist
         globIndx1 = newDist(iPair)%indx1
         globIndx2 = newDist(iPair)%indx2
         if (.not. rPair(globIndx1, globIndx2)%p%usePair) then
            cycle
         end if
         jType = atomIndicies(globIndx2)%nType
         jMol = atomIndicies(globIndx2)%nMol
         jIndx = MolArray(jType)%mol(jMol)%indx
         if (iIndx .ne. jIndx) then
            iAtom = atomIndicies(globIndx1)%nAtom
            jAtom = atomIndicies(globIndx2)%nAtom

            atmType1 = atomArray(iType, iAtom)
            atmType2 = atomArray(jType, jAtom)

            ep = ep_tab(atmType2, atmType1)
            q = q_tab(atmType2, atmType1)
            if (ep .ne. 0E0_dp) then
               r_sq = newDist(iPair)%r_sq
               if (r_sq .lt. lj_Cut_sq) then
                  sig_sq = sig_tab(atmType2, atmType1)
                  LJ = LJ_Func(r_sq, ep, sig_sq)
                  P_LJ = P_LJ + LJ
               end if
            end if

            if (q .ne. 0E0_dp) then
               r = newDist(iPair)%r
               Ele = q/r
               P_Ele = P_Ele + Ele
            end if
         end if
      end do

      P_Trial = P_LJ + P_Ele

   end subroutine
!======================================================================================
end module

