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

   type, abstract :: pure_compound
      !! Generic EoS
      character(len=:), allocatable :: name !! Compound name
      real(wp) :: ac !! Critical atractive parameter
      real(wp) :: b !! Repulsive parameter
      real(wp) :: tc !! Critical temperature
      real(wp) :: pc ! Critical pressure
      real(wp) :: w !! Accentric factor
      real(wp) :: k !! Atractive parameter constant
      real(wp) :: a = 0 !! Atractive parameter valuated at temperature
      real(wp) :: dadt = 0 !! Atractive parameter first derivative with tempetarue
      real(wp) :: da2dt2 = 0 !! Atractive parameter second derivative with tempetarue
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

end module cubic_eos
