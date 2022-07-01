module multicomponent_eos
   use constants
   use mixing_rules
   use cubic_eos
   implicit none

   type :: mixture
      class(pure_compound), allocatable :: compounds(:)
      class(quadratic), allocatable :: mixing_rule
   end type mixture

contains

   subroutine mix(self, T)
      class(mixture) :: self !! Mixture object
      real(wp), intent(in) :: T !! Temperature [K]

      class(kij_constant), pointer :: kij_rule
      class(lij_constant), pointer :: lij_rule

      real(wp), allocatable :: a(:) 
      real(wp), allocatable :: dadt(:)
      real(wp), allocatable :: da2dt2(:)
      real(wp), allocatable :: b(:), bij(:, :), lij(:, :)
      real(wp), allocatable :: kij(:, :), dkijdt(:, :), dkij2dt2(:, :)
      integer :: i, j, nc

      nc = size(self%compounds)
      ! =======================================================================
      !   Calculate all the pure compounds atractive and repulsive 
      !    parameters at T
      ! -----------------------------------------------------------------------
      allocate(a(nc))
      allocate(dadt(nc))
      allocate(da2dt2(nc))
      allocate(b(nc))
      do i=1,nc
         call self%compounds(i)%a_t(T)
         a(i) = self%compounds(i)%a
         dadt(i) = self%compounds(i)%dadt
         da2dt2(i) = self%compounds(i)%da2dt2
         
         b(i) = self%compounds(i)%b
      end do
      ! =======================================================================

      ! =======================================================================
      !  Calculate the kij and lij matrices
      ! -----------------------------------------------------------------------
      associate(kij_rule => self%kij_calculator)
      select type (kij_rule)
      type is (kij_constant)
         call kij_rule%get_kij(T, kij, dkijdt, dkij2dt2)

      type is (kij_exp_t)
         call kij_rule%get_kij(T, kij, dkijdt, dkij2dt2)
      
      class default
         print *, "Not implemented Kij rule"
      end select
      end associate

      associate(lij_rule => self%lij_calculator)
      select type(lij_rule)
      type is (lij_constant)
         call lij_rule%get_lij(T, lij)

      class default
         print *, "Not implemented Kij rule"
      end select

      self%kij_matrix = kij
      self%lij_matrix = lij
      end associate
      ! =======================================================================

      ! =======================================================================
      ! Calculate each aij
      ! -----------------------------------------------------------------------
      allocate(self%aij(nc, nc))
      allocate(self%bij(nc, nc))
      do i = 1, nc
         self%aij(i, :) = sqrt(a(i) * a)*(1 - kij(i, :))
         self%bij(i, :) = (1 - lij(i, :))*(b(i) + b(:))/2
      end do
      ! =======================================================================
   end subroutine mix
end module multicomponent_eos
