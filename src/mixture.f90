module multicomponent_eos
   use constants
   use mixing_rules
   use cubic_eos
   implicit none

   type :: mixture
      class(pure_compound), allocatable :: compounds(:)
      class(quadratic), allocatable :: mixing_rule
   end type mixture

end module multicomponent_eos
