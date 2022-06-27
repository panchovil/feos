module mixture
   use constants
   use mixing_rules
   use cubic_eos
   implicit none

   type :: mix
      class(pure_compound), allocatable :: compounds(:)
      class(kij_constant), allocatable :: kij_calculator
      class(lij_constant), allocatable :: lij_calculator
      real(wp), allocatable :: kij_matrix(:, :)
      real(wp), allocatable :: aij(:, :)
      real(wp), allocatable :: bij(:, :)
         contains
            procedure :: get_aij => aij
            !procedure :: get_kij => kij
   end type mix

contains

   subroutine aij(self, T)
      class(mix) :: self !! Mixture object
      real(wp), intent(in) :: T !! Temperature [K]

      class(kij_constant), pointer :: kij_rule
      real(wp), allocatable :: a(:)
      real(wp), allocatable :: dadt(:)
      real(wp), allocatable :: da2dt2(:)
      real(wp), allocatable :: b(:)
      real(wp), allocatable :: kij(:, :), dkijdt(:, :), dkij2dt2(:, :)

      integer :: i, j, nc

      nc = size(self%compounds)
      allocate(a(nc))
      allocate(dadt(nc))
      allocate(da2dt2(nc))
      allocate(b(nc))
      
      ! Calculate all the pure compounds atractive parameters at T
      do i=1,nc
         call self%compounds(i)%a_t(T)
         a(i) = self%compounds(i)%a
         dadt(i) = self%compounds(i)%dadt
         da2dt2(i) = self%compounds(i)%da2dt2
         
         b(i) = self%compounds(i)%b
      end do

      ! Calculate the kij matrix
      associate(kij_rule => self%kij_calculator)
      
      select type (kij_rule)
      type is (kij_constant)
         call kij_rule%get_kij(T, kij, dkijdt, dkij2dt2)

      type is (kij_exp_t)
         call kij_rule%get_kij(T, kij, dkijdt, dkij2dt2)
      
      class default
         print *, "Error"
      end select
      
      self%kij_matrix = kij
      end associate

      allocate(self%aij(nc, nc))
      do i = 1, nc
         self%aij(i, :) = sqrt(a(i) * a)*(1 - kij(i, :))
      end do
   end subroutine aij

end module mixture
