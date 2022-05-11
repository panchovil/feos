module mixing_rules
    !! Module that contains the available mixing rules to be used.
   use constants
   implicit none

contains

   subroutine quadratic(nc, a, b, kij, dadt, da2dt2, dkijdt, dkij2dt2, lij, aij, daijdt, daij2dt2, bij)
      !! Classic quadratic mixing rules.
      integer, intent(in) :: nc !! Number of components
      real(wp), intent(in) :: a(nc) !! Atractive parameter at working temperature
      real(wp), intent(in) :: b(nc) !! Repulsive parameter
      real(wp), intent(in) :: kij(nc, nc) !! Kij matrix
      real(wp), intent(in) :: dadt(nc) !! First derivative with T
      real(wp), intent(in) :: da2dt2(nc) !! Second derivative with T
      real(wp), intent(in) :: dkijdt(nc, nc) !! Kij matrix first derivative
      real(wp), intent(in) :: dkij2dt2(nc, nc) !! Kij matrix second derivative
      real(wp), intent(in) :: lij(nc, nc) !! Lij matrix

      real(wp), intent(out) :: aij(nc, nc) !! Binary atractive parameters matrix
      real(wp), intent(out) :: daijdt(nc, nc) !! First derivative with T
      real(wp), intent(out) :: daij2dt2(nc, nc) !! Second derivative with T
      real(wp), intent(out) :: bij(nc, nc) !! Repulse parameter matrix

      integer :: i, j

      do i = 1, nc
         aij(i, i) = a(i)
         daijdT(i, i) = dadT(i)
         daij2dT2(i, i) = da2dT2(i)
         do j = 1, nc
            aij(j, i) = sqrt(a(i)*a(j))*(1 - kij(j, i))
            aij(i, j) = aij(j, i)

            daijdt(j, i) = (1 - Kij(j, i))*(sqrt(a(i)/a(j))*dadT(j) + sqrt(a(j)/a(i))*dadT(i))/2 &
                           - dkijdt(j, i)*sqrt(a(j)*a(i))
            daijdt(i, j) = daijdT(j, i)

            daij2dt2(j, i) = (1 - Kij(j, i))*(dadt(j)*dadt(i)/sqrt(a(i)*a(j)) &
                                              + sqrt(a(i)/a(j))*(da2dt2(j) - dadt(j)**2/(2*a(j))) &
                                              + sqrt(a(j)/a(i))*(da2dt2(i) - dadt(i)**2/(2*a(i))))/2 &
                             - dkijdt(j, i)*(a(j)*dadt(i) + a(i)*dadt(j))/sqrt(a(j)*a(i)) &
                             - dkij2dt2(j, i)*sqrt(a(j)*a(i))
            daij2dT2(i, j) = daij2dT2(j, i)

            bij(i, j) = (1 - lij(i, j))*(b(i) + b(j))/2
            bij(j, i) = bij(i, j)
         end do
      end do
   end subroutine quadratic

   subroutine other(nc, a, b, kij, dadt, da2dt2, dkijdt, dkij2dt2, lij, aij, daijdt, daij2dt2, bij)
        !! What's this??
      integer, intent(in) :: nc !! Number of components
      real(wp), intent(in) :: a(nc) !! Atractive parameter at working temperature
      real(wp), intent(in) :: b(nc) !! Repulsive parameter
      real(wp), intent(in) :: kij(nc, nc) !! Kij matrix
      real(wp), intent(in) :: dadt(nc) !! First derivative with T
      real(wp), intent(in) :: da2dt2(nc) !! Second derivative with T
      real(wp), intent(in) :: dkijdt(nc, nc) !! Kij matrix first derivative
      real(wp), intent(in) :: dkij2dt2(nc, nc) !! Kij matrix second derivative
      real(wp), intent(in) :: lij(nc, nc) !! Lij matrix

      real(wp), intent(out) :: aij(nc, nc) !! Binary atractive parameters matrix
      real(wp), intent(out) :: daijdt(nc, nc) !! First derivative with T
      real(wp), intent(out) :: daij2dt2(nc, nc) !! Second derivative with T
      real(wp), intent(out) :: bij(nc, nc) !! Repulse parameter matrix

      integer :: i, j

      real(wp) :: barrgij

      call quadratic(nc, a, b, kij, dadt, da2dt2, dkijdt, dkij2dt2, lij, aij, daijdt, daij2dt2, bij)
      do i = 1, nc - 1
         do j = i + 1, nc
            barrgij = bij(i, j)/sqrt(b(i)*b(j))
            aij(i, j) = barrgij*aij(i, j)
            aij(j, i) = aij(i, j)
            daijdT(i, j) = barrgij*daijdT(i, j)
            daijdT(j, i) = daijdT(i, j)
            daij2dt2(i, j) = barrgij*daij2dt2(i, j)
            daij2dt2(j, i) = daij2dt2(i, j)
         end do
      end do
   end subroutine other
end module mixing_rules
