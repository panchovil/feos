module multicomponent_eos
   use constants
   use mixing_rules
   use cubic_eos
   implicit none

   type :: mixture
      class(pure_compound), allocatable :: compounds(:)
      class(kij_constant), allocatable :: kij_calculator
      class(lij_constant), allocatable :: lij_calculator
      real(wp), allocatable :: kij_matrix(:, :)
      real(wp), allocatable :: lij_matrix(:, :)
      real(wp), allocatable :: aij(:, :)
      real(wp), allocatable :: bij(:, :)
         contains
            procedure :: mix => mix
   end type mixture

contains

