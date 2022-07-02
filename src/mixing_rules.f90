module mixing_rules
   !! Module that contains the available mixing rules to be used.
   use constants
   use cubic_eos

   implicit none

   type :: quadratic
      !! Basic Mixing rule with constant \[K_ij\] and \[l_ij\]
      real(wp), allocatable :: kij(:, :) !! \[K_ij\] matrix
      real(wp), allocatable :: lij(:, :) !! \[l_ij\] matrix
      real(wp), allocatable :: aij(:, :) !! \[a_ij\] matrix
      real(wp), allocatable :: bij(:, :) !! \[b_ij\] matrix
    contains 
        procedure :: mix => quadratic_mix
   end type

   type, extends(quadratic) :: kij_exp_t
      !! Kij with temperature dependance according to the equation:
      !! \[ K_{ij}(T) = K_{ij\infty} + K_{ij0} e^{-T/T^*} \]
      !! The parameters of the equation are obtained from the mixture module.
      real(wp), allocatable :: dkijdt(:, :) !! \[K_ij\] matrix
      real(wp), allocatable :: dkij2dt2(:, :) !! \[K_ij\] matrix
      real(wp), allocatable :: kij_0(:, :) !! Exponential term
      real(wp), allocatable :: kij_inf(:, :) !! K_ij at infinite Temperature
      real(wp), allocatable :: T_star(:, :) !! Reference temperature
   contains
      procedure :: get_kij => kij_tdep
   end type

contains

   subroutine kij_tdep(self, T, kij, dkijdt, dkij2dt2)
      implicit none
      class(kij_exp_t) :: self
      real(wp), intent(in) :: T !! Temperature

      real(wp), allocatable, intent(out) :: kij(:, :) !! Binary interaction parameter matrix
      real(wp), allocatable, intent(out) :: dkijdt(:, :) !! Binary interaction parameter first derivative with T matrix
      real(wp), allocatable, intent(out) :: dkij2dt2(:, :) !! Binary interaction parameter second derivative with T matrix

      real(wp), allocatable :: kij_inf(:, :), kij_0(:, :), T_star(:, :)
      integer :: nc, sp(2)

      kij_0 = self%kij_0
      kij_inf = self%kij_inf
      T_star = self%T_star

      sp = shape(self%kij_0)
      nc = sp(1)

      allocate (kij(nc, nc))
      allocate (dkijdt(nc, nc))
      allocate (dkij2dt2(nc, nc))

      kij = kij_inf + kij_0*exp(-T/T_star)
      dkijdt = -kij_0/T_star*exp(-T/T_star)
      dkij2dt2 = kij_0/T_star**2*exp(-T/T_star)

      self%kij = kij
      self%dkijdt = dkijdt
      self%dkij2dt2 = dkij2dt2
   end subroutine kij_tdep

   subroutine quadratic_mix(self, T, compounds)
      class(quadratic) :: self !! Mixing rule
      class(pure_compound), allocatable, intent(in) :: compounds(:) !! Compounds to mix
      real(wp), intent(in) :: T !! Temperature [K]

      real(wp), allocatable :: a(:), aij(:, :)
      real(wp), allocatable :: dadt(:), daijdt(:, :)
      real(wp), allocatable :: da2dt2(:), daij2dt2(:, :)
      real(wp), allocatable :: b(:), bij(:, :), lij(:, :)
      real(wp), allocatable :: kij(:, :), dkijdt(:, :), dkij2dt2(:, :)
      integer :: i, nc

      nc = size(compounds)
      ! =======================================================================
      !   Calculate all the pure compounds atractive and repulsive
      !    parameters at T
      ! -----------------------------------------------------------------------
      allocate (a(nc))
      allocate (dadt(nc))
      allocate (da2dt2(nc))
      allocate (b(nc))

      do i = 1, nc
         call compounds(i)%a_t(T)
         a(i) = compounds(i)%a
         dadt(i) = compounds(i)%dadt
         da2dt2(i) = compounds(i)%da2dt2

         b(i) = compounds(i)%b
      end do
      ! =======================================================================

      ! =======================================================================
      !  Calculate the kij and lij matrices
      ! -----------------------------------------------------------------------
      select type(self)
         class is (kij_exp_t)
            call self%get_kij(T, kij, dkijdt, dkij2dt2)
            lij = self%lij

         class default
             kij = self%kij
             lij = self%lij
             dkijdt = 0*kij
             dkij2dt2 = 0*kij
      end select

      self%kij = kij
      self%lij = lij
      ! =======================================================================

      ! =======================================================================
      ! Calculate each aij
      ! -----------------------------------------------------------------------
      allocate(aij(nc, nc))
      allocate(bij(nc, nc))

      do i = 1, nc
         aij(i, :) = sqrt(a(i)*a)*(1 - kij(i, :))
         bij(i, :) = (1 - lij(i, :))*(b(i) + b(:))/2

         
      end do

      self%aij = aij
      self%bij = bij
      ! =======================================================================
   end subroutine quadratic_mix
end module mixing_rules
