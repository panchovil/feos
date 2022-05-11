module mixing_rules
    !! Module that contains the available mixing rules to be used.
    use constants
    implicit none

    contains

    subroutine quadratic(a, b, kij, dadt, da2dt2, dkijdt, dkij2dt2, lij, aij, daijdt, daij2dt2, bij, n)
        !! Classic quadratic mixing rules.
        real(wp), intent(in) :: a(n) !! Atractive parameter at working temperature
        real(wp), intent(in) :: b(n) !! Repulsive parameter
        real(wp), intent(in) :: kij(n, n) !! Kij matrix
        real(wp), intent(in) :: dadt(n) !! First derivative with T
        real(wp), intent(in) :: da2dt2(n) !! Second derivative with T
        real(wp), intent(in) :: dkijdt(n, n) !! Kij matrix first derivative
        real(wp), intent(in) :: dkij2dt2(n, n) !! Kij matrix second derivative
        real(wp), intent(in) :: lij(n, n) !! Lij matrix
        integer, intent(in) :: n !! Number of components

        real(wp), intent(out) :: aij(n, n) !! Binary atractive parameters matrix
        real(wp), intent(out) :: daijdt(n, n) !! First derivative with T
        real(wp), intent(out) :: daij2dt2(n, n) !! Second derivative with T
        real(wp), intent(out) :: bij(n, n) !! Repulse parameter matrix

        integer :: i, j

        do i = 1, n
            aij(i, i) = a(i)
            daijdT(i, i) = dadT(i)
            daij2dT2(i, i) = da2dT2(i)
            do j = 1, n
                aij(j, i) = sqrt(a(i)*a(j))*(1 - kij(j, i))
                aij(i, j) = aij(j, i)

                daijdt(j, i) = (1 - Kij(j, i))*(sqrt(a(i)/a(j))*dadT(j) + sqrt(a(j)/a(i))*dadT(i))/2 &
                                - dkijdt(j, i)*sqrt(a(j)*a(i))
                daijdt(i, j) = daijdT(j, i)

                daij2dt2(j, i) = (1 - Kij(j, i))*(dadt(j)*dadt(i)/sqrt(a(i)*a(j)) &
                                  + sqrt(a(i)/a(j))*(da2dt2(j) - dadt(j)**2/(2*a(j))) &
                                  + sqrt(a(j)/a(i))*(da2dt2(i) - dadt(i)**2/(2*a(i))))/2 &
                                  - dkijdt(j, i) * (a(j)*dadt(i) + a(i)*dadt(j))/sqrt(a(j)*a(i)) &
                                  - dkij2dt2(j, i) * sqrt(a(j)*a(i))
                daij2dT2(i, j) = daij2dT2(j, i)

                bij(i, j) = (1 - lij(i, j))*(b(i) + b(j))/2
                bij(j, i) = bij(i, j)
            end do
        end do
    end subroutine quadratic

    subroutine other(a, b, kij, dadt, da2dt2, dkijdt, dkij2dt2, lij, aij, daijdt, daij2dt2, bij, n)
        !! What's this??
        real(wp), intent(in) :: a(n) !! Atractive parameter at working temperature
        real(wp), intent(in) :: b(n) !! Repulsive parameter
        real(wp), intent(in) :: kij(n, n) !! Kij matrix
        real(wp), intent(in) :: dadt(n) !! First derivative with T
        real(wp), intent(in) :: da2dt2(n) !! Second derivative with T
        real(wp), intent(in) :: dkijdt(n, n) !! Kij matrix first derivative
        real(wp), intent(in) :: dkij2dt2(n, n) !! Kij matrix second derivative
        real(wp), intent(in) :: lij(n, n) !! Lij matrix
        integer, intent(in) :: n !! Number of components

        real(wp), intent(out) :: aij(n, n) !! Binary atractive parameters matrix
        real(wp), intent(out) :: daijdt(n, n) !! First derivative with T
        real(wp), intent(out) :: daij2dt2(n, n) !! Second derivative with T
        real(wp), intent(out) :: bij(n, n) !! Repulse parameter matrix
        
        integer :: i, j

        real(wp) :: barrgij

        call quadratic(a, b, kij, dadt, da2dt2, dkijdt, dkij2dt2, lij, aij, daijdt, daij2dt2, bij, n)
        do i = 1, n - 1
            do j = i + 1, n
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
