!==============================================================
!  This file contains a set of small functions that are either used in creating forcefield
!  tables such as mixing rules or contain functional forms that are used in the energy calculations
!  themselves.
!==============================================================
module ForceFieldVariableType
   use VarPrecision
!     This block defines custom variable types for forcefield
   type AtomDef
      character(len=5) :: Symb
      character(len=20) :: atmName
      real(dp) :: sig, ep, q, mass
   end type

   type AtomDefV
      character(len=20) :: atmName
      integer(kind=atomIntType) :: Symb1, Symb2
      real(dp) :: B, xi, r0, C, theta
   end type

   type BondDef
      character(len=20) :: bondName
      real(dp) :: k_eq, r_eq
   end type
   
   integer, parameter :: nBendHistBins = 200
   real(dp), parameter :: bendBinWidth = 4e0*atan(1e0)/real(nBendHistBins, dp)
   real(dp), parameter :: startProb = 0.05e0
   type BendAngleDef
      character(len=20) :: angleName
      real(dp) :: k_eq, ang_eq

      integer :: startBin
      real(dp) :: accptConstant
      real(dp) :: Prob(1:nBendHistBins)
   end type

   type TorsionAngleDef
      character(len=20) :: torsName
      integer(kind=atomIntType) :: nPara
      real(dp), allocatable :: a(:)
   end type

   type NonBondedIndex
!        integer :: nonMembr(1:2)
      integer(kind=atomIntType) :: nonMembr(1:2)
   end type

   type BondIndex
      integer(kind=atomIntType) :: bondType
      integer(kind=atomIntType) :: bondMembr(1:2)
!        integer :: bondType
!        integer :: bondMembr(1:2)
   end type

   type BendAngleIndex
      integer(kind=atomIntType) :: bendType
      integer(kind=atomIntType) :: bendMembr(1:3)
!        integer :: bendType
!        integer :: bendMembr(1:3)

   end type

   type TorsionAngleIndex
      integer(kind=atomIntType) :: TorsType
      integer(kind=atomIntType) :: torsMembr(1:4)
!        integer :: TorsType
!        integer :: torsMembr(1:4)
   end type

   type ImproperAngleIndex
      integer(kind=atomIntType) :: ImpropType
      integer(kind=atomIntType) :: impropMembr(1:4)
   end type

end module
!==============================================================
module ForceField
   use ForceFieldVariableType
   use VarPrecision

   character(len=20) :: ForceFieldName

   integer :: nAtomTypes = 0
   integer :: nAtomTypes2 = 0
   integer :: nAtomTypes3 = 0
   integer :: nBondTypes = 0
   integer :: nAngleTypes = 0
   integer :: nTorsionalTypes = 0
   integer :: nImproperTypes = 0

   type(AtomDef), allocatable :: atomData(:)
   type(AtomDefV), allocatable :: atomDataV(:)
   type(BondDef), allocatable :: bondData(:)
   type(BendAngleDef), allocatable :: bendData(:)
   type(TorsionAngleDef), allocatable :: torsData(:)
   type(TorsionAngleDef), allocatable :: impropData(:)

   real(dp), allocatable :: r_min(:), r_min_sq(:)
   real(dp), allocatable :: r_min_tab(:, :)

   integer, allocatable :: nAtoms(:)
   integer, allocatable :: nIntraNonBond(:)
   integer, allocatable :: nBonds(:)
   integer, allocatable :: nAngles(:)
   integer, allocatable :: nTorsional(:)
   integer, allocatable :: nImproper(:)

   real(dp), allocatable :: totalMass(:)

   integer(kind=atomIntType), allocatable :: atomArray(:, :)
   type(NonBondedIndex), allocatable :: nonBondArray(:, :)
   type(BondIndex), allocatable :: bondArray(:, :)
   type(BendAngleIndex), allocatable :: bendArray(:, :)
   type(TorsionAngleIndex), allocatable :: torsArray(:, :)
   type(ImproperAngleIndex), allocatable :: impropArray(:, :)

!      character(len=16), allocatable :: atomUseByBond(:,:)
!      character(len=66), allocatable :: atomUseByBend(:,:)
!      character(len=66), allocatable :: atomUseByTors(:,:)
!      character(len=66), allocatable :: atomUseByImprop(:,:)

end module
!==============================================================

module ForceFieldFunctions
   interface
      real(dp) function MixRule(par1, par2)
         use VarPrecision
         real(dp), Intent(IN) :: par1, par2
      end function
   end interface
   interface
      real(dp) function TorsEng(angle, coeffArray)
         use VarPrecision
         real(dp), Intent(IN) :: angle
         real(dp), Intent(IN) :: coeffArray(:)
      end function
   end interface
contains
!        !----------------------------------------------------------------
!       !Geometric Mean Mixing rule
   real(dp) function GeoMean_MixingFunc(par1, par2)
      use VarPrecision
      real(dp), Intent(IN) :: par1, par2

      if ((par1 .eq. 0) .or. (par2 .eq. 0)) then
         GeoMean_MixingFunc = 0e0
      else
         GeoMean_MixingFunc = sqrt(par1*par2)
      end if
   end function
!       !----------------------------------------------------------------
!      !Standard Mean(average) Mixing rule
   real(dp) function Mean_MixingFunc(par1, par2)
      use VarPrecision
      real(dp), Intent(IN) :: par1, par2
      if ((par1 .eq. 0e0) .or. (par2 .eq. 0e0)) then
         Mean_MixingFunc = 0e0
      else
         Mean_MixingFunc = (par1 + par2)*0.5e0
      end if

   end function
!       !----------------------------------------------------------------
!      !(Fill in Later) Mixing rule
   real(dp) function MeanMax_MixingFunc(par1, par2)
      use VarPrecision
      real(dp), Intent(IN) :: par1, par2
      if ((par1 .gt. 0e0) .and. (par2 .gt. 0e0)) then
         MeanMax_MixingFunc = abs((par1 + par2)*0.5e0)
      else
         MeanMax_MixingFunc = max(abs(par1), abs(par2))
      end if

   end function
!       !----------------------------------------------------------------
!      !Minimum Mixing rule
   real(dp) function Min_MixingFunc(par1, par2)
      use VarPrecision
      real(dp), Intent(IN) :: par1, par2
      Min_MixingFunc = min(par1, par2)
   end function
!       !----------------------------------------------------------------
!      !Maximum Mixing rule
   real(dp) function Max_MixingFunc(par1, par2)
      use VarPrecision
      real(dp), Intent(IN) :: par1, par2
      Max_MixingFunc = max(par1, par2)
   end function
!       !----------------------------------------------------------------
!       ! Standard Harmonic Potential Energy Function
   pure real(dp) function Harmonic(angle, k_eq, ang_eq)
      use VarPrecision
      real(dp), Intent(IN) :: angle, k_eq, ang_eq
      Harmonic = 0.5e0*k_eq*(angle - ang_eq)**2
   end function
!       !----------------------------------------------------------------
!       ! Standard Harmonic Potential Energy Function that accepts a Torsional Parameter Array
   pure real(dp) function TorsHarmonic(angle, coeffArray)
      use VarPrecision
      real(dp), Intent(IN) :: angle
      real(dp), Intent(IN) :: coeffArray(:)

      TorsHarmonic = 0.5e0*coeffArray(1)*(angle - coeffArray(2))**2
   end function
!       !----------------------------------------------------------------
!       ! Standard Torsional Potential using a series of Cos^n functionals
   pure real(dp) function Torsion_Cos_N(angle, coeffArray)
      use VarPrecision
      integer :: nCoeff, i
      real(dp), Intent(IN) :: angle
      real(dp), Intent(IN) :: coeffArray(:)
      real(dp) :: eng

      nCoeff = size(coeffArray)
      eng = coeffArray(1)
      do i = 2, nCoeff
         if (coeffArray(i) .ne. 0e0) then
            eng = eng + coeffArray(i)*cos(angle)**(i - 1)
         end if
      end do
      Torsion_Cos_N = eng
   end function
!       !----------------------------------------------------------------
!       ! Torsional Potential using a series of Cos(nx) functions
   pure real(dp) function Torsion_CosNx(angle, coeffArray)
      use VarPrecision
      integer :: nCoeff, i
      real(dp), Intent(IN) :: angle
      real(dp), Intent(IN) :: coeffArray(:)
      real(dp) :: eng

      nCoeff = size(coeffArray)
      eng = coeffArray(1)
      do i = 2, nCoeff
         if (coeffArray(i) .ne. 0e0) then
            eng = eng + coeffArray(i)*cos(dble(i - 1)*angle)
         end if
      end do
      Torsion_CosNx = eng
   end function
!       !----------------------------------------------------------------
!       ! Trappe Torsional Potential using a series of Cos(nx) functions
   pure real(dp) function Trappe_CosNx(angle, coeffArray)
      use VarPrecision
      integer :: nCoeff
      real(dp), Intent(IN) :: angle
      real(dp), Intent(IN) :: coeffArray(:)
      real(dp) :: eng

      nCoeff = size(coeffArray)
      eng = coeffArray(1)
      eng = eng + coeffArray(2)*(1e0 + cos(angle))
      eng = eng + coeffArray(3)*(1e0 - cos(2e0*angle))
      eng = eng + coeffArray(4)*(1e0 + cos(3e0*angle))
      Trappe_CosNx = eng
   end function
!       !----------------------------------------------------------------

end module
!==============================================================
module ForceFieldPara_Vashishta
   use ForceFieldFunctions
   use ForceFieldVariableType
   use VarPrecision

   procedure(MixRule), pointer :: ep_func => null()
   procedure(MixRule), pointer :: sig_func => null()
   real(dp), allocatable :: ep_tab(:, :), sig_tab(:, :)
   real(dp), allocatable :: q_tab(:, :)

   real(dp), allocatable :: H_tab(:, :), eta_tab(:, :), W_tab(:, :)
   real(dp), allocatable :: r1s_tab(:, :), D_tab(:, :), r4s_tab(:, :)
   ! real(dp), allocatable :: rc_tab(:, :)
end module
!==============================================================
module ForceFieldPara_LJ_Q
   use ForceFieldFunctions
   use ForceFieldVariableType
   use VarPrecision

   procedure(MixRule), pointer :: ep_func => null()
   procedure(MixRule), pointer :: sig_func => null()
   real(dp), allocatable :: ep_tab(:, :), sig_tab(:, :)
   real(dp), allocatable :: q_tab(:, :)

end module
!==============================================================
module ForceFieldPara_Pedone
   use VarPrecision

   type AtomDefPedone
      character(len=20) :: atmName
      character(len=5) :: Symb
      real(dp) :: repul, rEq, q, alpha, delta, mass
   end type

   logical :: implcSolvent = .false.

   type(AtomDefPedone), allocatable :: pedoneData(:)
   real(dp), allocatable :: repul_tab(:, :), rEq_tab(:, :)
   real(dp), allocatable :: q_tab(:, :), alpha_Tab(:, :), D_Tab(:, :)
   real(dp), allocatable :: bornRad(:)

end module
!==============================================================
module ForceFieldPara_Tersoff
   use VarPrecision

   type AtomDefTersoff
      character(len=20) :: atmName
      character(len=5) :: Symb
      real(dp) :: A, B, c, d, n, lam1, lam2, h, R, D2, beta, mass
   end type

   type(AtomDefTersoff), allocatable :: tersoffData(:)

end module
!==============================================================
