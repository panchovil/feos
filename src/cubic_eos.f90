module cubic_eos
   !! Module that encompass the calculations of the residual Helmholtz energy
   !! and related properties like fugacity coefficents.
   use constants
   implicit none

   private
   public :: pure_compound
   public :: pr
   public :: srk
   public :: rkpr

   type :: pure_compound
      !! Generic EoS
      character(len=:), allocatable :: name !! Compound name
      real(8) :: ac !! Critical atractive parameter
      real(8) :: b !! Repulsive parameter
      real(8) :: tc !! Critical temperature
      real(8) :: pc !! Critical pressure
      real(8) :: w !! Accentric factor
      real(8) :: k !! Atractive parameter constant
      real(8) :: a = 0 !! Atractive parameter valuated at temperature
      real(8) :: dadt = 0 !! Atractive parameter first derivative with tempetarue
      real(8) :: da2dt2 = 0 !! Atractive parameter second derivative with tempetarue
   contains
      procedure :: a_t => a_parameter
   end type pure_compound

   type, extends(pure_compound) :: pr
      !! Peng-Robinson EoS
      real(8) :: del1 = 1 !! \[\delta_1\] parameter.
   end type pr

   type, extends(pure_compound) :: srk
      !! Soave-Redlich-Kwong EoS
      real(8) :: del1 = 1 + sqrt(2.d0) !! \[\delta_1\] parameter.
   end type srk

   type, extends(pure_compound) :: rkpr
      real(8) :: del1
   end type rkpr

contains

   subroutine a_parameter(self, T)
      !! Calculate the atractive parameter at T temperature.
      !! the subroutine will read the mixture's model and based on that
      !! will use the corresponding rule.
      class(pure_compound) :: self
      real(8), intent(in) :: T !! Temperature where to calculate

      real(8) :: Tr
      real(8) :: ac, k, Tc

      Tc = self%Tc
      ac = self%ac
      k = self%k
      Tr = T/Tc

      select type (self)
      class is (rkpr)
         self%a = ac*(3/(2 + Tr))**k
         self%dadT = -k*self%a/Tc/(2 + Tr)
         self%da2dT2 = -(k + 1)*self%dadT/Tc/(2 + Tr)

      class default
         self%a = ac*(1 + k*(1 - sqrt(Tr)))**2
         self%dadT = ac*k*(k - (k + 1)/sqrt(Tr))/Tc
         self%da2dT2 = ac*k*(k + 1)/(2*Tc**2*Tr**1.5)
      end select
   end subroutine a_parameter

   !subroutine d1nder(nc, n, d1, d1_mix, dD1dni, dD12dnij2)
   !   !! delta_1 of the mixture and it's compositional derivatives.
   !   integer, intent(in) :: nc !! Number of components.
   !   real(wp), intent(in) :: n(nc) !! Array of mole numbers for each component.
   !   real(wp), intent(in) :: d1(nc) !! Array of delta 1 parameters for each component.
   !   real(wp), intent(out) :: d1_mix !! Mixture's delta 1
   !   real(wp), intent(out) :: dD1dni(nc) !! delta1 parameter first derivative with composition.
   !   real(wp), intent(out) :: dD12dnij2(nc, nc) !! delta1 parameter second derivative with composition.

   !   integer :: i, j
   !   real(wp) :: totn ! Total number of moles

   !   d1_mix = 0.0_wp
   !   d1_mix = sum(n*d1)

   !   totn = sum(n)
   !   D1_mix = D1_mix/totn
   !   dD1dni = (d1 - d1_mix)/totn

   !   do i = 1, nc
   !      do j = 1, nc
   !         dD12dnij2(i, j) = (2.0_wp*d1_mix - d1(i) - d1(j))/totn**2
   !      end do
   !   end do
   !end subroutine d1nder

   !subroutine Bmixnder(nc, n, bij, B_mix, dBdni, dB2dnij2)
   !   !! Repulsive parameter of the mixture and it's compositional derivatives.
   !   integer, intent(in) :: nc !! Number of components
   !   real(wp), intent(in) :: n(nc) !! Number of moles of each component
   !   real(wp), intent(in) :: bij(nc, nc) !! Repulsive parameter matrix

   !   real(wp), intent(out) :: B_mix !! Mixture repulsive parameter
   !   real(wp), intent(out) :: dBdni(nc) !! Repulsive parameter derivatives wrt number of moles
   !   real(wp), intent(out) :: dB2dnij2(nc, nc) !! Repulsive parameter second derivatives wrt number of moles

   !   real(wp) :: totn, aux(nc)
   !   integer :: i, j

   !   totn = sum(n)
   !   B_mix = 0.0_wp
   !   aux = 0.0_wp

   !   do i = 1, nc
   !      do j = 1, nc
   !         aux(i) = aux(i) + n(j)*bij(i, j)
   !      end do
   !      B_mix = B_mix + n(i)*aux(i)
   !   end do

   !   B_mix = B_mix/totn

   !   do i = 1, nc
   !      dBdni(i) = (2*aux(i) - B_mix)/totn
   !      do j = 1, i
   !         dB2dnij2(i, j) = (2*bij(i, j) - dBdni(i) - dBdni(j))/totn
   !         dB2dnij2(j, i) = dB2dnij2(i, j)
   !      end do
   !   end do
   !end subroutine Bmixnder
end module cubic_eos
