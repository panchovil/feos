module thermo
   !! Module that encompass the calculations of the residual Helmholtz energy
   !! and related properties like fugacity coefficents.
    use constants
    implicit none

    contains

    subroutine aTder(T, a, dadT, dadT2)
        !! Calculate the atractive parameter at T temperature.
        !! the subroutine will read the mixture's model and based on that
        !! will use the corresponding rule.
        use mixture, only: nmodel, ac, k, tc, nc
        real(wp), intent(in) :: T  !! Temperature

        real(wp), intent(out) :: a(nc)  !! Atractive parameter at T
        real(wp), intent(out) :: dadT(nc) !! First deritvative with T
        real(wp), intent(out) :: dadT2(nc)  !! Second derivative with T

        real(wp) :: Tr(nc) ! Reduced temperature

        Tr = T/Tc
        if (nmodel <= 2) then
            ! SRK and PR
            a = ac*(1 + k*(1 - sqrt(Tr)))**2
            dadT = ac*k*(k - (k + 1)/sqrt(Tr))/Tc
            dadT2 = ac*k*(k + 1)/(2*Tc**2*Tr**1.5)
        else if (nmodel == 3) then
            ! RKPR EOS
            a = ac*(3/(2 + Tr))**k
            dadT = -k*a/Tc/(2 + Tr)
            dadT2 = -(k + 1)*dadT/Tc/(2 + Tr)
        end if
    end subroutine aTder

    subroutine aijTder(NTD, T, aij, daijdT, daij2dT2)
        !! Calculate the binary atractive term matrix
        use mixture, only: nc, ncomb, ntdep, kij, kij_inf, kij0, T_star, lij
        implicit none

        integer, intent(in) :: NTD !! Calculate (or not) temperature derivatives (1 for yes)
        real(wp), intent(in) :: T !! Temperature

        real(wp), intent(out) :: aij(nc, nc) !! Atractive binary terms matrix
        real(wp), intent(out) :: daijdT(nc, nc) !! Atractive binary terms matrix first derivative with temperature
        real(wp), intent(out) :: daij2dT2(nc, nc) !! Atractive binary terms matrix second derivative with temperature

        real(wp) :: a(nc), dadT(nc), da2dT2(nc)
        real(wp) :: aux(nc, nc), ratK(nc, nc) ! ????
        real(wp) :: Tr(nc)
        real(wp) :: barrgij(nc, nc)

        integer :: i, j

        ! =============================================================================
        !  Setting up the Kij matrix
        ! -----------------------------------------------------------------------------
        if (ntdep >= 1) then
            kij = 0
            do i = 1, nc
                kij(:i - 1, i) = kij_inf(:i - 1, i) + kij0(:i - 1, i)*exp(-T/T_star(:i - 1, i))
            end do
        else
            kij = kij0
        end if
        ! =============================================================================

        ! =============================================================================
        !  Set up the aij and derivatives matrices
        ! -----------------------------------------------------------------------------
        call aTder(T, a, dadT, da2dT2)
        call classic_vdw(a, b, dadt, da2dt2, kij, lij, aij, daijdt, daij2dt2, bij, nc)
        do i = 1, nc
            aij(i, i) = a(i)
            daijdT(i, i) = dadT(i)
            daij2dT2(i, i) = da2dT2(i)

            if (i > 1) then
                do j = 1, i - 1
                    aij(j, i) = sqrt(a(i)*a(j))*(1 - Kij(j, i))
                    aij(i, j) = aij(j, i)

                    if (NTD == 1) then
                        daijdT(j, i) = (1 - Kij(j, i))*(sqrt(a(i)/a(j))*dadT(j) &
                                       + sqrt(a(j)/a(i))*dadT(i))/2
                        daijdT(i, j) = daijdT(j, i)
                        daij2dT2(j, i) = (1 - Kij(j, i))*(dadT(j)*dadT(i)/sqrt(a(i)*a(j)) &
                                         + sqrt(a(i)/a(j))*(da2dT2(j) - dadT(j)**2/(2*a(j))) &
                                         + sqrt(a(j)/a(i))*(da2dT2(i) - dadT(i)**2/(2*a(i))))/2
                        daij2dT2(i, j) = daij2dT2(j, i)
                    end if
                end do
            end if
        end do
        ! =============================================================================

        !! TODO: What combining rule is this?
        !if (ncomb == 1) then
        !    do i = 1, nc - 1
        !        do j = i + 1, nc
        !            barrgij = bij(i, j)/sqrt(b(i)*b(j))

        !            aij(i, j) = barrgij*aij(i, j)
        !            aij(j, i) = aij(i, j)

        !            daijdT(i, j) = barrgij*daijdT(i, j)
        !            daijdT(j, i) = daijdT(i, j)

        !            daij2dT2(i, j) = barrgij*daijdT2(i, j)
        !            daij2dT2(j, i) = daij2dT2(i, j)
        !        end do
        !    end do
        !end if

        !if (ntdep .ge. 1 .and. NTD .EQ. 1) then
        !    do i = 1, nc
        !        aux(:i - 1, i) = daijdT(:i - 1, i)
        !        ratK(:i - 1, i) = Kij(:i - 1, i)/(1 - Kij(:i - 1, i))/Tstar(:i - 1, i)
        !        daijdT(:i - 1, i) = aux(:i - 1, i) + aij(:i - 1, i)*ratK(:i - 1, i)
        !        daijdT(i, :i - 1) = daijdT(:i - 1, i)
        !        daij2dT2(:i - 1, i) = daij2dT2(:i - 1, i) + (2*aux(:i - 1, i) - aij(:i - 1, i)/Tstar(:i - 1, i))*ratK(:i - 1, i)
        !        daij2dT2(i, :i - 1) = daij2dT2(:i - 1, i)
        !    end do
        !end if
    end subroutine aijTder
end module thermo

module mixing_rules
    use constants
    use thermo
    implicit none

    contains

    subroutine classic_vdw(a, b, dadt, da2dt2, kij, lij, aij, daijdt, daij2dt2, bij, n)
        real(wp), intent(in) :: a(n)
        real(wp), intent(in) :: dadt(n)
        real(wp), intent(in) :: da2dt2(n)
        real(wp), intent(in) :: b(n)
        real(wp), intent(in) :: kij(n, n)
        real(wp), intent(in) :: lij(n, n)

        real(wp), intent(out) :: aij(n, n)
        real(wp), intent(out) :: daijdt(n, n)
        real(wp), intent(out) :: daij2dt2(n, n)
        real(wp), intent(out) :: bij(n, n)

        integer, intent(in) :: n

        integer :: i, j

        do i = 1, n
            aij(i, i) = a(i)
            daijdT(i, i) = dadT(i)
            daij2dT2(i, i) = da2dT2(i)
            do j = 1, n
                aij(j, i) = sqrt(a(i)*a(j))*(1 - kij(j, i))
                aij(i, j) = aij(j, i)

                daijdT(j, i) = (1 - Kij(j, i))*(sqrt(a(i)/a(j))*dadT(j) &
                                + sqrt(a(j)/a(i))*dadT(i))/2
                daijdT(i, j) = daijdT(j, i)

                daij2dT2(j, i) = (1 - Kij(j, i))*(dadT(j)*dadT(i)/sqrt(a(i)*a(j)) &
                                  + sqrt(a(i)/a(j))*(da2dT2(j) - dadT(j)**2/(2*a(j))) &
                                  + sqrt(a(j)/a(i))*(da2dT2(i) - dadT(i)**2/(2*a(i))))/2
                daij2dT2(i, j) = daij2dT2(j, i)

                bij(i, j) = (1 - lij(i, j))*(b(i) + b(j))/2
                bij(j, i) = bij(i, j)
            end do
        end do
    end subroutine classic_vdw
end module mixing_rules
