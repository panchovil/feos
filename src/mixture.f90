module multicomponent_eos
   use constants, only: wp
   use mixing_rules, only: ConstantBIP
   use cubic_eos, only: CubicEoS
   
   implicit none

   private

   public :: mixture

   type :: mixture(n)
      integer, len :: n
      real(wp) :: P
      real(wp) :: V
      real(wp) :: T
      type(CubicEoS) :: as_a_fluid !! As a fluid representation of the mixture
      type(CubicEoS) :: compounds(n)
      !real(wp), allocatable :: concentrations(:)
      type(ConstantBIP) :: mixing_rule(n)
   end type mixture
end module multicomponent_eos
