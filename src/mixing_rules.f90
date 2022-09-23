module mixing_rules
   !! Module that contains the available mixing rules to be used.
   use constants, only: wp, R
   use properties
   use cubic_eos

   implicit none

   private

   public :: ConstantBIP
   public :: ClassicVdW

   type :: ConstantBIP(n)
      !! Classic Mixing rules with constant K_ij and l_ij
      integer, len :: n
      type(scalar_property(n)) :: D
      type(scalar_property(n)) :: B

      type(binary_property(n)) :: kij !! Binary interaction parameter
      type(binary_property(n)) :: lij !! Binary repulsive parameter
      type(binary_property(n)) :: aij !! Binary atractive parameter
      type(binary_property(n)) :: bij !! Binary repulsive parameter
    !contains
      !procedure :: mix => quadratic_mix    !! Atractive and repulsive matrices calculation
      !procedure :: d_mix => attractive_mix !! Mixture's Attractive parameter calculation
      !procedure :: b_mix => repulsive_mix  !! Mixture's Repulsive parameter calculation
      !procedure :: d1_mix => delta1_mix    !! Mixture's Delta 1 calculation
   end type ConstantBIP
contains

   ! ==========================================================================
   !  Constant kij and lij
   ! --------------------------------------------------------------------------
   function ClassicVdW(n, kij, lij) result(mixing_rule)
      !! Classic "Van der Waals" Mixing Rules, with constant Binary interaction parameters
      integer, intent(in) :: n
      real(wp), intent(in) :: kij(n, n)
      real(wp), intent(in) :: lij(n, n)
      type(ConstantBIP(n)) :: mixing_rule

      type(scalar_property(n)) :: D_in
      type(scalar_property(n)) :: B_in
      type(binary_property(n)) :: kij_in !! Binary interaction parameter
      type(binary_property(n)) :: lij_in !! Binary repulsive parameter
      type(binary_property(n)) :: aij_in !! Binary atractive parameter
      type(binary_property(n)) :: bij_in !! Binary repulsive parameter

      kij_in%val = kij
      lij_in%val = lij

      mixing_rule = ConstantBIP(&
         n=n, D=D_in, B=B_in, kij=kij_in, lij=lij_in, aij=aij_in, bij=bij_in&
      )
   end function ClassicVdW

   function aij_classic(n, a, kij) result(aij)
      integer, intent(in) :: n
      type(scalar_property(n)), intent(in) :: a(n)  !! Atractive parameters
      type(binary_property(n)), intent(in) :: kij  !! kij

      type(binary_property(n)) :: aij !! aij

      integer :: j

      do j=1,n
         aij%val(j, j) = a(n)%val
         aij%val(:, j) = sqrt(a(:)%val * a(j)%val) * (1 - kij%val(:, j))
      end do

   end function aij_classic

   function bij_classic(n, b, lij) result(bij)
      integer, intent(in) :: n
      type(scalar_property(n)), intent(in) :: b(n)        !! Repulsive parameters
      type(binary_property(n)), intent(in) :: lij         !! lij matrix
      type(binary_property(n)) :: bij

      integer :: j

      do j=1, n
         bij%val(j, j) = b(j)%val
         bij%val(:, j) = (b(:)%val + b(j)%val)/2 * (1 - lij%val(:, j))
      end do
   end function bij_classic
   ! ==========================================================================

   ! =============================================================================
   !  Kij Exponentialy dependant of temperature and constant lij
   ! -----------------------------------------------------------------------------
   
   ! =============================================================================

!   subroutine kij_tdep(self, T)
!      !! Kij with temperature dependance according to the equation:
!      !! \[ K_{ij}(T) = K_{ij\infty} + K_{ij0} e^{-T/T^*} \]
!      !! The parameters of the equation are obtained from the mixture module.
!      implicit none
!      class(kij_exp_t) :: self
!      real(wp), intent(in) :: T !! Temperature
!
!      real(wp), allocatable :: kij_inf(:, :), kij_0(:, :), T_star(:, :)
!      integer :: nc, sp(2)
!
!      kij_0 = self%kij_0
!      kij_inf = self%kij_inf
!      T_star = self%T_star
!
!      sp = shape(self%kij_0)
!      nc = sp(1)
!
!      if (allocated(self%kij))      deallocate(self%kij)
!      if (allocated(self%dkijdt))   deallocate(self%dkijdt)
!      if (allocated(self%dkij2dt2)) deallocate(self%dkij2dt2)
!
!      allocate(self%kij(nc, nc))
!      allocate(self%dkijdt(nc, nc))
!      allocate(self%dkij2dt2(nc, nc))
!
!      self%kij = kij_inf + kij_0*exp(-T/T_star)
!      self%dkijdt = -kij_0/T_star*exp(-T/T_star)
!      self%dkij2dt2 = kij_0/T_star**2*exp(-T/T_star)
!   end subroutine kij_tdep
!
!   subroutine aij_classic(self, temperature)
!       type(quadratic) :: self
!       real(wp), intent(in) :: temperature
!       aij(i, :) = sqrt(a(i)*a)*(1 - kij(i, :))
!       daijdt(i, :) = - (kij(i, :) - 1) * (a * dadt(i) + a(i) * dadt) * (a(i)*a)**(-1._wp/2._wp)/2._wp
!       daij2dt2(i, :) = (1 - kij(i, :)) * ((a*da2dt2(i) + 2._wp * dadt(i)*dadt + a(i)) * (a(i)*a)**(-1._wp/2._wp)/2._wp) - (a*dadt(i) + a(i)*dadt)**(3._wp/4._wp)
!
!
!       bij(i, :) = (1 - lij(i, :))*(b(i) + b(:))/2
!   end subroutine
!
!   subroutine quadratic_mix(self, T, compounds, concentrations)
!      !! Calculate aij, bij and kij matrices.
!      class(quadratic) :: self !! Mixing rule
!      real(wp), intent(in) :: T !! Temperature [K]
!      class(pure_compound), allocatable, intent(in) :: compounds(:) !! Compounds to mix
!      real(wp), allocatable, intent(in) :: concentrations(:) !! Number of moles of each component
!
!      real(wp), allocatable :: a(:), aij(:, :)
!      real(wp), allocatable :: dadt(:), daijdt(:, :)
!      real(wp), allocatable :: da2dt2(:), daij2dt2(:, :)
!
!      real(wp), allocatable :: b(:), bij(:, :), lij(:, :)
!      
!      real(wp), allocatable :: kij(:, :), dkijdt(:, :), dkij2dt2(:, :)
!      
!      real(wp), allocatable :: d1(:)
!
!      real(wp), allocatable :: aux(:)
!
!      integer :: i, nc
!
!      nc = size(compounds)
!      
!      ! =======================================================================
!      !   Calculate all the pure compounds atractive and repulsive
!      !    parameters at T
!      ! -----------------------------------------------------------------------
!      allocate (a(nc))
!      allocate (dadt(nc))
!      allocate (da2dt2(nc))
!      allocate (b(nc))
!
!      do i = 1, nc
!         call compounds(i)%a_t(T)
!         a(i) = compounds(i)%a
!         dadt(i) = compounds(i)%dadt
!         da2dt2(i) = compounds(i)%da2dt2
!
!         b(i) = compounds(i)%b
!      end do
!      ! =======================================================================
!
!      ! =======================================================================
!      !  Calculate the kij and lij matrices
!      ! -----------------------------------------------------------------------
!      select type(self)
!         class is (kij_exp_t)
!            call self%get_kij(T)
!            dkijdt = self%dkijdt
!            dkij2dt2 = self%dkij2dt2
!
!         class default
!             dkijdt = 0*self%kij
!             dkij2dt2 = 0*self%kij
!      end select
!
!      kij = self%kij
!      lij = self%lij
!      ! =======================================================================
!
!      ! =======================================================================
!      ! Calculate each aij
!      ! -----------------------------------------------------------------------
!      allocate(aij(nc, nc))
!      allocate(daijdt(nc, nc))
!      allocate(daij2dt2(nc, nc))
!      allocate(bij(nc, nc))
!
!      do i = 1, nc
!      end do
!
!      self%aij = aij
!      self%daijdt = daijdt
!      self%daij2dt2 = daij2dt2
!      self%bij = bij
!      ! =======================================================================
!
!      ! =======================================================================
!      !  Calculate mixture attractive and repulsive parameter, 
!      !  and their derivatives
!      ! -----------------------------------------------------------------------
!      select type(compounds)
!      type is (rkpr)
!          allocate(aux(nc))
!          do i = 1, nc
!             aux(i) = compounds(i)%del1
!          end do
!
!          call self%d_mix(concentrations)
!          call self%b_mix(concentrations)
!          call self%d1_mix(concentrations, d1)
!
!      class default
!          call self%d_mix(concentrations)
!          call self%b_mix(concentrations)
!      end select
!      ! =======================================================================
!   end subroutine quadratic_mix
!
!   subroutine attractive_mix(self, concentrations)
!      !! Attractive term
!      class(quadratic) :: self
!      real(wp), allocatable, intent(in) :: concentrations(:) !! Vector of concentrations of the mixture
!
!      integer :: i, j, nc
!      real(wp) :: D, dDdt, dD2dT2, aux, aux2
!      real(wp), allocatable :: dDdni(:), dD2dnit2(:), dD2dnij2(:, :)
!      real(wp), allocatable :: aij(:, :), daijdt(:, :), daij2dt2(:, :)
!      real(wp), allocatable :: n(:)
!
!      n = concentrations
!      nc = size(n)
!
!      D = 0.0_wp
!      dDdt = 0.0_wp
!      dD2dt2 = 0.0_wp
!
!      allocate(dDdni(nc))
!      allocate(dD2dnit2(nc))
!      allocate(dD2dnij2(nc, nc))
!
!      aij = self%aij
!      daijdt = self%daijdt
!      daij2dt2 = self%daij2dt2
!
!
!      do i = 1, nc
!         aux = 0.0_wp
!         aux2 = 0.0_wp
!         dDdni(i) = 0.0_wp
!         dD2dniT2(i) = 0.0_wp
!         do j = 1, nc
!            dDdni(i) = dDdni(i) + 2*n(j)*aij(i, j)
!            dD2dniT2(i) = dD2dniT2(i) + 2*n(j)*daijdT(i, j)
!            dD2dnij2(i, j) = 2*aij(i, j)
!
!            aux = aux + n(j)*aij(i, j)
!            aux2 = aux2 + n(j)*daij2dT2(i, j)
!         end do
!         D = D + n(i)*aux
!         dDdT = dDdT + n(i)*dD2dniT2(i)/2
!         dD2dT2 = dD2dT2 + n(i)*aux2
!      end do
!
!      self%D = D 
!      self%dDdt = dDdt 
!      self%dD2dt2 = dD2dt2 
!               
!      self%dDdni = dDdni
!      self%dD2dnit2 = dD2dnit2
!      self%dD2dnij2 = dD2dnij2
!   end subroutine attractive_mix
!
!   subroutine repulsive_mix(self, concentrations)
!       !! Repulsive parameter of the mixture and it's compositional derivatives.
!       class(quadratic) :: self
!       real(wp), allocatable :: concentrations(:)
!
!       real(wp) :: B_mix ! Mixture repulsive parameter
!       real(wp), allocatable :: dBdni(:) ! Repulsive parameter derivatives wrt number of moles
!       real(wp), allocatable :: dB2dnij2(:, :) ! Repulsive parameter second derivatives wrt number of moles
!
!       real(wp), allocatable :: aux(:)
!       real(wp) :: totn
!       real(wp), allocatable :: n(:)
!       real(wp), allocatable :: bij(:, :)
!
!       integer :: i, j, nc
!
!       n = concentrations
!       nc = size(n)
!
!       allocate(aux(nc))
!       allocate(dBdni(nc))
!       allocate(dB2dnij2(nc, nc))
!
!       bij = self%bij
!
!       totn = sum(concentrations)
!       B_mix = 0.0_wp
!       aux = 0.0_wp
!
!       do i = 1, nc
!          do j = 1, nc
!             aux(i) = aux(i) + n(j)*bij(i, j)
!          end do
!          B_mix = B_mix + n(i)*aux(i)
!       end do
!
!       B_mix = B_mix/totn
!
!       do i = 1, nc
!          dBdni(i) = (2*aux(i) - B_mix)/totn
!          do j = 1, i
!             dB2dnij2(i, j) = (2*bij(i, j) - dBdni(i) - dBdni(j))/totn
!             dB2dnij2(j, i) = dB2dnij2(i, j)
!          end do
!       end do
!
!       self%B = B_mix
!       self%dBdni = dBdni
!       self%dB2dnij2 = dB2dnij2
!   end subroutine repulsive_mix
!
!   subroutine delta1_mix(self, concentrations, d1)
!      !! Delta 1 parameter and compositional derivatives
!      class(quadratic) :: self
!      real(wp), allocatable, intent(in) :: concentrations(:) !! Array of concentrations
!      real(wp), allocatable, intent(in) :: d1(:) !! Array of delta_1 parameters
!
!      real(wp) :: totn ! Total number of moles
!      real(wp) :: d1_mix ! Mixture's d1
!      real(wp), allocatable :: n(:) ! Concentrations copy
!      real(wp), allocatable :: dD1dni(:)
!      real(wp), allocatable :: dD12dnij2(:, :)
!      
!      integer :: i, j, nc
!
!      n = concentrations
!      nc = size(n)
!      totn = sum(n)
!
!      d1_mix = sum(n*d1)/totn
!      dD1dni = (d1 - d1_mix)/totn
!
!      allocate(dD12dnij2(nc, nc))
!
!      do i = 1, nc
!         do j = 1, nc
!            dD12dnij2(i, j) = (2.0_wp*d1_mix - d1(i) - d1(j))/totn**2.0_wp
!         end do
!      end do
!
!      self%D1 = d1_mix
!      self%dD1dni = dD1dni
!      self%dD12dnij2 = dD12dnij2
!   end subroutine delta1_mix
end module mixing_rules
