module multicomponent_eos
   use constants
   use mixing_rules
   use cubic_eos

   implicit none
   private

   public :: mixture

   type :: mixture
      class(pure_compound), allocatable :: compounds(:)
      real(wp), allocatable :: concentrations(:)
      class(quadratic), allocatable :: mixing_rule
   end type mixture
end module multicomponent_eos
